import type { Analysis, AngleRecord, AtomAnalysis, Molecule } from './types';
import type { Vec3 } from './vec';
import { v3, add, sub, mul, norm, len, bondAngle, dot } from './vec';
import { element } from './periodic';
import { neighborLists } from './smiles';
import { findRings } from './rings';
import { applyResonance } from './resonance';
import { buildGeometry, type AtomEnvironment } from './build3d';
import {
  arrangeGroups, axeNotation, bondWeight, geometryFor, hybridFor,
  lonePairCount, WEIGHT_LONE_PAIR,
} from './vsepr';
import { angleLabel, molecularMass, plainFormula, prettyFormula } from './formula';

/** Масштаб перевода качественной оценки дипольного момента в дебаи. */
const DIPOLE_SCALE = 1.25;
/** Вклад одной неподелённой пары в дипольный момент. */
const DIPOLE_LONE_PAIR = 0.50;

/**
 * Полный разбор молекулы: электронное строение каждого атома,
 * трёхмерная структура, валентные углы, полярность.
 */
export function analyze(mol: Molecule): Analysis {
  // Усредняем кратности равноценных связей: нитрат-ион должен получиться
  // правильным треугольником, а не «двойная плюс две одинарные».
  applyResonance(mol);

  const nb = neighborLists(mol);
  const n = mol.atoms.length;

  // --- сумма кратностей связей у каждого атома ------------------------------
  const orderSum = new Array(n).fill(0);
  const bondOrderTo: Map<string, number> = new Map();
  for (const b of mol.bonds) {
    orderSum[b.a] += b.order;
    orderSum[b.b] += b.order;
    bondOrderTo.set(`${b.a}:${b.b}`, b.order);
    bondOrderTo.set(`${b.b}:${b.a}`, b.order);
  }

  // --- окружение каждого атома по теории Гиллеспи ---------------------------
  const env: AtomEnvironment[] = [];
  for (let id = 0; id < n; id++) {
    const atom = mol.atoms[id];
    const neighbors = nb[id];
    const bondOrders = neighbors.map((k) => bondOrderTo.get(`${id}:${k}`) ?? 1);
    const lonePairs = lonePairCount(atom.el, atom.charge, orderSum[id]);
    const steric = neighbors.length + lonePairs;

    const weights = [
      ...bondOrders.map(bondWeight),
      ...new Array(lonePairs).fill(WEIGHT_LONE_PAIR),
    ];
    env.push({ neighbors, bondOrders, lonePairs, steric, localDirs: arrangeGroups(weights) });
  }

  const rings = findRings(mol);
  const build = buildGeometry(mol, env, rings);

  // --- разбор каждого атома -------------------------------------------------
  const ringAtoms = new Set<number>();
  rings.forEach((r) => r.atoms.forEach((a) => ringAtoms.add(a)));

  const atoms: AtomAnalysis[] = [];
  for (let id = 0; id < n; id++) {
    const atom = mol.atoms[id];
    const e = env[id];
    const sigma = e.neighbors.length;
    const geom = geometryFor(sigma, e.lonePairs);

    const orbitalDirs = e.localDirs.map((d) => normWorld(build.frames[id], d));

    atoms.push({
      id,
      el: atom.el,
      charge: atom.charge,
      sigma,
      lonePairs: e.lonePairs,
      steric: e.steric,
      axe: axeNotation(sigma, e.lonePairs),
      hybrid: hybridFor(e.steric),
      electronGeom: geom.electronGeom,
      molecularGeom: geom.molecularGeom,
      idealAngle: geom.idealAngle,
      predictedAngle: predictedAngleFor(e),
      isCentral: sigma >= 2,
      lonePairDirs: build.lonePairDirs[id],
      orbitalDirs,
      warning: warningFor(atom.el, e.steric, e.lonePairs, ringAtoms.has(id)),
    });
  }

  // --- валентные углы -------------------------------------------------------
  const expTable = normalizeExperimental(mol.meta?.exp);
  const angles: AngleRecord[] = [];
  for (let id = 0; id < n; id++) {
    const e = env[id];
    if (e.neighbors.length < 2) continue;
    for (let a = 0; a < e.neighbors.length; a++) {
      for (let b = a + 1; b < e.neighbors.length; b++) {
        const i = e.neighbors[a];
        const j = e.neighbors[b];
        const label = angleLabel(mol.atoms[i].el, mol.atoms[id].el, mol.atoms[j].el);
        const localAngle = angleBetween(e.localDirs[a], e.localDirs[b]);
        angles.push({
          i,
          center: id,
          j,
          label,
          actual: bondAngle(mol.atoms[i].pos, mol.atoms[id].pos, mol.atoms[j].pos),
          predicted: localAngle,
          experimental: null,
          inRing: ringAtoms.has(id) && ringAtoms.has(i) && ringAtoms.has(j),
        });
      }
    }
  }
  attachExperimental(angles, expTable);

  // --- дипольный момент (качественная оценка) -------------------------------
  let dipoleVec = v3();
  for (const b of mol.bonds) {
    const enA = element(mol.atoms[b.a].el).en;
    const enB = element(mol.atoms[b.b].el).en;
    const direction = norm(sub(mol.atoms[b.b].pos, mol.atoms[b.a].pos));
    dipoleVec = add(dipoleVec, mul(direction, enB - enA));
  }
  for (let id = 0; id < n; id++) {
    for (const d of build.lonePairDirs[id]) {
      dipoleVec = add(dipoleVec, mul(d, DIPOLE_LONE_PAIR));
    }
  }
  dipoleVec = mul(dipoleVec, DIPOLE_SCALE);
  const dipoleValue = len(dipoleVec);

  // --- сводка по гибридизациям ----------------------------------------------
  const hybridSummary: Record<string, number> = {};
  for (const a of atoms) {
    // концевые атомы в сводку не попадают: интересна гибридизация центров
    if (a.sigma < 2 || a.hybrid === '—') continue;
    hybridSummary[a.hybrid] = (hybridSummary[a.hybrid] ?? 0) + 1;
  }

  return {
    molecule: mol,
    formula: plainFormula(mol),
    formulaHtml: prettyFormula(mol),
    mass: molecularMass(mol),
    atoms,
    angles,
    rings,
    dipoleVec,
    dipoleValue,
    polar: dipoleValue > 0.06,
    hybridSummary,
  };
}

function normWorld(frame: number[], local: Vec3): Vec3 {
  return norm(v3(
    frame[0] * local.x + frame[1] * local.y + frame[2] * local.z,
    frame[3] * local.x + frame[4] * local.y + frame[5] * local.z,
    frame[6] * local.x + frame[7] * local.y + frame[8] * local.z,
  ));
}

function angleBetween(a: Vec3, b: Vec3): number {
  return (Math.acos(Math.max(-1, Math.min(1, dot(norm(a), norm(b))))) * 180) / Math.PI;
}

/** Наиболее характерный предсказанный угол между связями данного атома. */
function predictedAngleFor(e: AtomEnvironment): number | null {
  if (e.neighbors.length < 2) return null;
  const values: number[] = [];
  for (let a = 0; a < e.neighbors.length; a++) {
    for (let b = a + 1; b < e.neighbors.length; b++) {
      values.push(angleBetween(e.localDirs[a], e.localDirs[b]));
    }
  }
  values.sort((x, y) => x - y);
  // берём наименьший — именно его обычно и приводят в справочниках
  return values[0];
}

/** Приведение экспериментальных подписей к тому же виду, что и вычисленные. */
function normalizeExperimental(
  exp: Record<string, number | number[]> | undefined,
): Map<string, number[]> {
  const table = new Map<string, number[]>();
  if (!exp) return table;
  for (const [rawKey, value] of Object.entries(exp)) {
    const parts = rawKey.split('-');
    const key = parts.length === 3 ? angleLabel(parts[0], parts[1], parts[2]) : rawKey;
    table.set(key, Array.isArray(value) ? value.slice() : [value]);
  }
  return table;
}

/**
 * Привязка экспериментальных значений к конкретным углам.
 *
 * Подпись «F-S-F» в гексафториде серы относится сразу к пятнадцати углам:
 * двенадцати по 90° и трём по 180°. Поэтому каждое справочное значение
 * привязывается к тому ещё не занятому углу, предсказание для которого
 * к нему ближе всего.
 */
function attachExperimental(angles: AngleRecord[], table: Map<string, number[]>): void {
  for (const [label, values] of table) {
    const candidates = angles.filter((a) => a.label === label);
    if (candidates.length === 0) continue;
    for (const value of values) {
      let best: AngleRecord | null = null;
      let bestDiff = Infinity;
      for (const c of candidates) {
        if (c.experimental !== null) continue;
        const diff = Math.abs((c.predicted ?? 0) - value);
        if (diff < bestDiff) { bestDiff = diff; best = c; }
      }
      if (best) best.experimental = value;
    }
  }
}

/**
 * Предупреждения о случаях, где модель ОЭПВО заведомо расходится с опытом.
 * Это не ошибка программы, а известная граница применимости теории.
 */
function warningFor(el: string, steric: number, lonePairs: number, inRing: boolean): string | undefined {
  const period = element(el).period;

  if (period >= 3 && lonePairs > 0 && steric === 4 && el !== 'Cl' && el !== 'Br' && el !== 'I') {
    return 'У элемента 3-го периода и ниже s- и p-орбитали сильно различаются по энергии, гибридизация выражена слабо. Реальный угол заметно ближе к 90°, чем предсказывает модель.';
  }
  if (inRing && steric === 4) {
    return 'Атом входит в цикл: валентный угол задан геометрией кольца и может сильно отличаться от предсказанного. Разница и есть угловое напряжение цикла.';
  }
  if (steric >= 5 && lonePairs > 0) {
    return 'Для гипервалентных частиц модель верно определяет ФОРМУ, но тонкое искажение углов неподелёнными парами она передаёт лишь качественно.';
  }
  return undefined;
}
