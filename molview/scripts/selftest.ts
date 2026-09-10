/**
 * Самопроверка химического движка на всей базе молекул.
 * Запуск: npm run test
 *
 * Проверяется:
 *   • каждая молекула базы разбирается и строится без ошибок;
 *   • длины связей соответствуют сумме ковалентных радиусов;
 *   • атомы не налезают друг на друга;
 *   • предсказанные валентные углы сравниваются с экспериментальными;
 *   • контрольные молекулы имеют ожидаемую форму и гибридизацию.
 */
import { database, moleculeFromEntry } from '../src/chem/database';
import { analyze } from '../src/chem/analyze';
import { dist } from '../src/chem/vec';
import { bondLength } from '../src/chem/periodic';
import type { Analysis } from '../src/chem/types';

let failures = 0;
const problems: string[] = [];

function fail(message: string): void {
  failures++;
  problems.push(message);
}

const pad = (s: string, n: number): string => (s.length >= n ? s.slice(0, n) : s + ' '.repeat(n - s.length));
const padLeft = (s: string, n: number): string => (s.length >= n ? s : ' '.repeat(n - s.length) + s);

// ---------------------------------------------------------------------------
// 1. Полный прогон базы
// ---------------------------------------------------------------------------
console.log('\n=== РАЗБОР БАЗЫ ===\n');
console.log(pad('Молекула', 24) + pad('Формула', 12) + pad('Гибридизации', 22) + pad('Циклы', 7) + 'μ, Д');
console.log('-'.repeat(78));

const results: { name: string; analysis: Analysis }[] = [];

for (const entry of database()) {
  let analysis: Analysis;
  try {
    analysis = analyze(moleculeFromEntry(entry));
  } catch (err) {
    fail(`${entry.name}: разбор завершился ошибкой — ${(err as Error).message}`);
    continue;
  }
  results.push({ name: entry.name, analysis });

  const hyb = Object.entries(analysis.hybridSummary)
    .sort((a, b) => b[1] - a[1])
    .map(([h, c]) => `${h}×${c}`)
    .join(' ');

  console.log(
    pad(entry.name, 24) +
    pad(analysis.formulaHtml, 12) +
    pad(hyb || '—', 22) +
    pad(String(analysis.rings.length), 7) +
    analysis.dipoleValue.toFixed(2),
  );

  // --- геометрические проверки ---
  const mol = analysis.molecule;
  for (const b of mol.bonds) {
    const expected = bondLength(mol.atoms[b.a].el, mol.atoms[b.b].el, b.order);
    const actual = dist(mol.atoms[b.a].pos, mol.atoms[b.b].pos);
    if (Math.abs(actual - expected) > 0.12) {
      fail(`${entry.name}: связь ${mol.atoms[b.a].el}-${mol.atoms[b.b].el} = ${actual.toFixed(2)} Å, ожидалось ${expected.toFixed(2)} Å`);
    }
  }
  for (let i = 0; i < mol.atoms.length; i++) {
    for (let j = i + 1; j < mol.atoms.length; j++) {
      const d = dist(mol.atoms[i].pos, mol.atoms[j].pos);
      if (d < 0.85) {
        fail(`${entry.name}: атомы ${mol.atoms[i].el}${i} и ${mol.atoms[j].el}${j} слиплись (${d.toFixed(2)} Å)`);
      }
    }
  }
  for (const a of mol.atoms) {
    if (!Number.isFinite(a.pos.x + a.pos.y + a.pos.z)) {
      fail(`${entry.name}: у атома ${a.el}${a.id} нечисловые координаты`);
    }
  }
}

// ---------------------------------------------------------------------------
// 2. Сравнение с экспериментом
// ---------------------------------------------------------------------------
console.log('\n\n=== ПРЕДСКАЗАНИЕ ПРОТИВ ЭКСПЕРИМЕНТА ===\n');
console.log(
  pad('Молекула', 24) + pad('Угол', 11) + padLeft('ОЭПВО', 8) +
  padLeft('модель', 8) + padLeft('опыт', 8) + padLeft('Δ', 7) + '  примечание',
);
console.log('-'.repeat(84));

let compared = 0;
let sumAbsError = 0;
let ringCases = 0;
let heavyCases = 0;
const worst = { name: '', label: '', delta: 0 };

for (const { name, analysis } of results) {
  for (const angle of analysis.angles) {
    if (angle.experimental === null || angle.predicted === null) continue;

    const delta = angle.predicted - angle.experimental;
    const centerAtom = analysis.atoms[angle.center];
    const heavy = !angle.inRing && centerAtom.warning !== undefined;

    if (angle.inRing) ringCases++;
    else if (heavy) heavyCases++;
    else {
      compared++;
      sumAbsError += Math.abs(delta);
      if (Math.abs(delta) > Math.abs(worst.delta)) { worst.name = name; worst.label = angle.label; worst.delta = delta; }
    }

    let note = '';
    if (angle.inRing) note = 'угол задан циклом';
    else if (heavy) note = 'известная граница модели';
    else if (Math.abs(delta) > 6) note = '← расхождение';

    console.log(
      pad(name, 24) + pad(angle.label, 11) +
      padLeft(angle.predicted.toFixed(1), 8) +
      padLeft(angle.actual.toFixed(1), 8) +
      padLeft(angle.experimental.toFixed(1), 8) +
      padLeft((delta >= 0 ? '+' : '') + delta.toFixed(1), 7) + '  ' + note,
    );
  }
}

console.log('-'.repeat(84));
console.log(`Сравнений в области применимости модели: ${compared}, средняя ошибка ${(sumAbsError / compared).toFixed(2)}°`);
console.log(`Наибольшее расхождение: ${worst.name}, ${worst.label}, ${worst.delta.toFixed(1)}°`);
console.log(`Отдельно учтены: углы в циклах — ${ringCases}, известные границы модели — ${heavyCases}`);

// ---------------------------------------------------------------------------
// 3. Контрольные молекулы: форма и гибридизация
// ---------------------------------------------------------------------------
console.log('\n\n=== КОНТРОЛЬНЫЕ ПРОВЕРКИ ===\n');

interface Expectation {
  id: string;
  /** Индекс центрального атома (в порядке разбора SMILES). */
  center: number;
  hybrid: string;
  molecularGeom: string;
  angle?: [number, number];
  rings?: number;
  planar?: boolean;
}

const EXPECTATIONS: Expectation[] = [
  { id: 'water', center: 0, hybrid: 'sp³', molecularGeom: 'угловая', angle: [103, 106] },
  { id: 'ammonia', center: 0, hybrid: 'sp³', molecularGeom: 'тригонально-пирамидальная', angle: [105, 108] },
  { id: 'methane', center: 0, hybrid: 'sp³', molecularGeom: 'тетраэдрическая', angle: [109, 110] },
  { id: 'carbon-dioxide', center: 1, hybrid: 'sp', molecularGeom: 'линейная', angle: [179, 181] },
  { id: 'boron-trifluoride', center: 1, hybrid: 'sp²', molecularGeom: 'плоская треугольная', angle: [119, 121] },
  { id: 'ethene', center: 0, hybrid: 'sp²', molecularGeom: 'плоская треугольная' },
  { id: 'ethyne', center: 0, hybrid: 'sp', molecularGeom: 'линейная', angle: [179, 181] },
  { id: 'allene', center: 1, hybrid: 'sp', molecularGeom: 'линейная', angle: [179, 181] },
  { id: 'formaldehyde', center: 0, hybrid: 'sp²', molecularGeom: 'плоская треугольная', angle: [115, 119] },
  { id: 'sulfur-dioxide', center: 1, hybrid: 'sp²', molecularGeom: 'угловая', angle: [115, 120] },
  { id: 'phosphorus-pentachloride', center: 1, hybrid: 'sp³d', molecularGeom: 'тригонально-бипирамидальная', angle: [89, 91] },
  { id: 'sulfur-tetrafluoride', center: 1, hybrid: 'sp³d', molecularGeom: 'качели (дисфеноид)' },
  { id: 'chlorine-trifluoride', center: 1, hybrid: 'sp³d', molecularGeom: 'Т-образная', angle: [89, 91] },
  { id: 'xenon-difluoride', center: 1, hybrid: 'sp³d', molecularGeom: 'линейная', angle: [179, 181] },
  { id: 'xenon-tetrafluoride', center: 1, hybrid: 'sp³d²', molecularGeom: 'квадратная', angle: [89, 91] },
  { id: 'sulfur-hexafluoride', center: 1, hybrid: 'sp³d²', molecularGeom: 'октаэдрическая', angle: [89, 91] },
  { id: 'iodine-pentafluoride', center: 1, hybrid: 'sp³d²', molecularGeom: 'квадратно-пирамидальная' },
  { id: 'benzene', center: 0, hybrid: 'sp²', molecularGeom: 'плоская треугольная', rings: 1, planar: true },
  { id: 'naphthalene', center: 0, hybrid: 'sp²', molecularGeom: 'плоская треугольная', rings: 2, planar: true },
  { id: 'cyclohexane', center: 0, hybrid: 'sp³', molecularGeom: 'тетраэдрическая', rings: 1, planar: false },
  { id: 'cyclopropane', center: 0, hybrid: 'sp³', molecularGeom: 'тетраэдрическая', rings: 1 },
  { id: 'pyridine', center: 3, hybrid: 'sp²', molecularGeom: 'угловая', rings: 1, planar: true },
  { id: 'nitrate', center: 1, hybrid: 'sp²', molecularGeom: 'плоская треугольная', angle: [119, 121] },
  { id: 'sulfate', center: 1, hybrid: 'sp³', molecularGeom: 'тетраэдрическая', angle: [109, 110] },
  { id: 'ammonium', center: 0, hybrid: 'sp³', molecularGeom: 'тетраэдрическая', angle: [109, 110] },
  { id: 'dimethyl-sulfoxide', center: 1, hybrid: 'sp³', molecularGeom: 'тригонально-пирамидальная' },
  { id: 'glucose', center: 2, hybrid: 'sp³', molecularGeom: 'тетраэдрическая', rings: 1 },
];

for (const exp of EXPECTATIONS) {
  const entry = database().find((e) => e.id === exp.id);
  if (!entry) { fail(`Контроль: молекула «${exp.id}» отсутствует в базе`); continue; }

  const analysis = analyze(moleculeFromEntry(entry));
  const atom = analysis.atoms[exp.center];
  const notes: string[] = [];

  if (atom.hybrid !== exp.hybrid) {
    fail(`${entry.name}: гибридизация ${atom.hybrid}, ожидалась ${exp.hybrid}`);
    notes.push(`гибридизация ${atom.hybrid} ≠ ${exp.hybrid}`);
  }
  if (atom.molecularGeom !== exp.molecularGeom) {
    fail(`${entry.name}: форма «${atom.molecularGeom}», ожидалась «${exp.molecularGeom}»`);
    notes.push(`форма «${atom.molecularGeom}»`);
  }
  if (exp.angle) {
    const centered = analysis.angles.filter((a) => a.center === exp.center);
    const minAngle = Math.min(...centered.map((a) => a.actual));
    if (minAngle < exp.angle[0] || minAngle > exp.angle[1]) {
      fail(`${entry.name}: минимальный угол ${minAngle.toFixed(1)}°, ожидался в диапазоне ${exp.angle[0]}–${exp.angle[1]}°`);
      notes.push(`угол ${minAngle.toFixed(1)}°`);
    }
  }
  if (exp.rings !== undefined && analysis.rings.length !== exp.rings) {
    fail(`${entry.name}: найдено циклов ${analysis.rings.length}, ожидалось ${exp.rings}`);
    notes.push(`циклов ${analysis.rings.length}`);
  }
  if (exp.planar !== undefined) {
    const ringIds = new Set(analysis.rings.flatMap((r) => r.atoms));
    const heavy = analysis.molecule.atoms.filter((a) => ringIds.has(a.id));
    const spread = maxPlaneDeviation(heavy.map((a) => a.pos));
    const isPlanar = spread < 0.12;
    if (isPlanar !== exp.planar) {
      fail(`${entry.name}: кольцо ${isPlanar ? 'плоское' : 'неплоское'} (разброс ${spread.toFixed(2)} Å), ожидалось ${exp.planar ? 'плоское' : 'неплоское'}`);
      notes.push(`разброс ${spread.toFixed(2)} Å`);
    }
  }

  const status = notes.length === 0 ? 'ок' : 'ОШИБКА: ' + notes.join(', ');
  console.log(pad(entry.name, 26) + pad(atom.axe, 8) + pad(atom.hybrid, 7) + pad(atom.molecularGeom, 30) + status);
}

/** Максимальное отклонение точек от наилучшей плоскости. */
function maxPlaneDeviation(points: { x: number; y: number; z: number }[]): number {
  if (points.length < 4) return 0;
  const c = points.reduce((s, p) => ({ x: s.x + p.x, y: s.y + p.y, z: s.z + p.z }), { x: 0, y: 0, z: 0 });
  const n = points.length;
  const center = { x: c.x / n, y: c.y / n, z: c.z / n };

  let xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
  for (const p of points) {
    const d = { x: p.x - center.x, y: p.y - center.y, z: p.z - center.z };
    xx += d.x * d.x; xy += d.x * d.y; xz += d.x * d.z;
    yy += d.y * d.y; yz += d.y * d.z; zz += d.z * d.z;
  }
  // нормаль — собственный вектор с наименьшим собственным значением; ищем грубым перебором
  let best = Infinity;
  let bestNormal = { x: 0, y: 0, z: 1 };
  for (let a = 0; a < 60; a++) {
    for (let b = 0; b < 60; b++) {
      const theta = (Math.PI * a) / 60;
      const phi = (2 * Math.PI * b) / 60;
      const nx = Math.sin(theta) * Math.cos(phi);
      const ny = Math.sin(theta) * Math.sin(phi);
      const nz = Math.cos(theta);
      const value =
        xx * nx * nx + yy * ny * ny + zz * nz * nz +
        2 * (xy * nx * ny + xz * nx * nz + yz * ny * nz);
      if (value < best) { best = value; bestNormal = { x: nx, y: ny, z: nz }; }
    }
  }
  let maxDev = 0;
  for (const p of points) {
    const d = (p.x - center.x) * bestNormal.x + (p.y - center.y) * bestNormal.y + (p.z - center.z) * bestNormal.z;
    maxDev = Math.max(maxDev, Math.abs(d));
  }
  return maxDev;
}

// ---------------------------------------------------------------------------
console.log('\n' + '='.repeat(78));
if (failures === 0) {
  console.log(`ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ. Молекул в базе: ${database().length}.`);
} else {
  console.log(`ОБНАРУЖЕНО ПРОБЛЕМ: ${failures}\n`);
  problems.forEach((p) => console.log('  • ' + p));
  process.exitCode = 1;
}
