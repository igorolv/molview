import type { Molecule, Ring } from './types';
import type { Vec3, Mat3 } from './vec';
import {
  v3, add, sub, mul, dot, cross, norm, len, dist, clone,
  applyMat, rotationAxis, rotationFromTo, alignPair, centroid, jacobiEigen, IDENTITY,
} from './vec';
import { bondLength } from './periodic';
import { fuseRingSystems } from './rings';

/**
 * Построение трёхмерной структуры молекулы.
 *
 * Порядок работы:
 *   1. Циклы укладываются по геометрическим шаблонам (плоский многоугольник
 *      для ароматических, «кресло» для насыщенного шестичленного и т. д.).
 *   2. Остальные атомы наращиваются обходом графа в ширину: для каждого
 *      атома известны направления его электронных групп в СОБСТВЕННОЙ системе
 *      координат, остаётся найти поворот, совмещающий их с уже размещёнными
 *      соседями.
 *   3. При наращивании подбирается поворот вокруг связи (конформация) так,
 *      чтобы атомы не налезали друг на друга — получается заторможенная
 *      конформация этана и анти-конформация бутана.
 *   4. Лёгкая релаксация разводит оставшиеся сближения.
 *   5. Молекула разворачивается по главным осям — к зрителю самой широкой
 *      стороной.
 */

/** Данные об окружении атома, нужные построителю. */
export interface AtomEnvironment {
  neighbors: number[];
  bondOrders: number[];
  lonePairs: number;
  steric: number;
  /** Направления групп в локальной системе: сначала связи (в порядке neighbors), затем пары. */
  localDirs: Vec3[];
}

export interface BuildResult {
  /** Поворот из локальной системы атома в мировую. */
  frames: Mat3[];
  /** Мировые направления неподелённых пар для каждого атома. */
  lonePairDirs: Vec3[][];
}

export function buildGeometry(mol: Molecule, env: AtomEnvironment[], rings: Ring[]): BuildResult {
  const n = mol.atoms.length;
  const frames: Mat3[] = new Array(n).fill(null).map(() => IDENTITY.slice() as Mat3);
  const frameKnown = new Array<boolean>(n).fill(false);
  for (const atom of mol.atoms) { atom.placed = false; atom.pos = v3(); }

  // Атомы, уложенные по геометрическому шаблону кольца. Только их имеет смысл
  // держать неподвижными при релаксации: их геометрия заведомо правильная.
  // Мостиковые и прочие сложные системы шаблоном не покрываются, их достраивает
  // обход графа — и такие атомы релаксация обязана поправлять.
  const templatePlaced = new Set<number>();
  placeRings(mol, env, rings, templatePlaced);

  // --- наращивание остальной молекулы ---------------------------------------
  const queue: number[] = mol.atoms.filter((a) => a.placed).map((a) => a.id);

  if (queue.length === 0) {
    const root = pickRoot(mol, env);
    mol.atoms[root].pos = v3();
    mol.atoms[root].placed = true;
    queue.push(root);
  }

  let guard = 0;
  while (guard++ < n * 8) {
    if (queue.length === 0) {
      // отдельная несвязанная часть молекулы — отодвигаем её в сторону
      const next = mol.atoms.find((a) => !a.placed);
      if (!next) break;
      const shift = spanOf(mol) + 3;
      next.pos = v3(shift, 0, 0);
      next.placed = true;
      queue.push(next.id);
    }

    const id = queue.shift()!;
    const e = env[id];
    const pending = e.neighbors.filter((k) => !mol.atoms[k].placed);
    if (pending.length === 0 && frameKnown[id]) continue;

    const rot = resolveFrame(mol, env, id);
    frames[id] = rot;
    frameKnown[id] = true;

    for (const nbId of pending) {
      const slot = e.neighbors.indexOf(nbId);
      const dir = norm(applyMat(rot, e.localDirs[slot]));
      const L = bondLength(mol.atoms[id].el, mol.atoms[nbId].el, e.bondOrders[slot]);
      mol.atoms[nbId].pos = add(mol.atoms[id].pos, mul(dir, L));
      mol.atoms[nbId].placed = true;
      queue.push(nbId);
    }
  }

  // атомы, до которых обход не добрался (защита от зависания)
  for (const atom of mol.atoms) {
    if (!atom.placed) { atom.pos = v3(spanOf(mol) + 2, 0, 0); atom.placed = true; }
  }

  relax(mol, env, templatePlaced);

  // --- направления неподелённых пар в мировых координатах --------------------
  const lonePairDirs: Vec3[][] = [];
  for (let id = 0; id < n; id++) {
    if (!frameKnown[id]) {
      frames[id] = resolveFrame(mol, env, id);
      frameKnown[id] = true;
    }
    const e = env[id];
    const dirs: Vec3[] = [];
    for (let k = e.neighbors.length; k < e.localDirs.length; k++) {
      dirs.push(norm(applyMat(frames[id], e.localDirs[k])));
    }
    lonePairDirs.push(dirs);
  }

  orientByPrincipalAxes(mol);

  return { frames, lonePairDirs };
}

/** Корневой атом: самый «нагруженный» тяжёлый атом — от него удобнее расти. */
function pickRoot(mol: Molecule, env: AtomEnvironment[]): number {
  let best = 0;
  let bestScore = -1;
  for (const atom of mol.atoms) {
    const score = (atom.el === 'H' ? 0 : 100) + env[atom.id].neighbors.length;
    if (score > bestScore) { bestScore = score; best = atom.id; }
  }
  return best;
}

function spanOf(mol: Molecule): number {
  let max = 0;
  for (const a of mol.atoms) if (a.placed) max = Math.max(max, len(a.pos));
  return max;
}

// ---------------------------------------------------------------------------
// Поворот локальной системы атома в мировую
// ---------------------------------------------------------------------------

function resolveFrame(mol: Molecule, env: AtomEnvironment[], id: number): Mat3 {
  const e = env[id];
  const me = mol.atoms[id];
  const placed: { slot: number; local: Vec3; world: Vec3 }[] = [];

  for (let slot = 0; slot < e.neighbors.length; slot++) {
    const nbId = e.neighbors[slot];
    if (!mol.atoms[nbId].placed) continue;
    const world = sub(mol.atoms[nbId].pos, me.pos);
    if (len(world) < 1e-6) continue;
    placed.push({ slot, local: e.localDirs[slot], world: norm(world) });
  }

  if (placed.length === 0) return IDENTITY.slice() as Mat3;

  if (placed.length === 1) {
    const base = rotationFromTo(placed[0].local, placed[0].world);
    return chooseTwist(mol, env, id, base, placed[0].world);
  }

  // выбираем пару направлений с наибольшим углом между ними — так устойчивее
  let bi = 0, bj = 1, bestSep = -2;
  for (let i = 0; i < placed.length; i++) {
    for (let j = i + 1; j < placed.length; j++) {
      const sep = -dot(placed[i].world, placed[j].world);
      if (sep > bestSep) { bestSep = sep; bi = i; bj = j; }
    }
  }

  // Совмещаем по биссектрисе и нормали: это симметрично и не «перекашивает»
  // окружение, когда реальный угол в цикле отличается от идеального.
  let localSum = v3(), worldSum = v3();
  for (const p of placed) { localSum = add(localSum, p.local); worldSum = add(worldSum, p.world); }

  const localNormal = cross(placed[bi].local, placed[bj].local);
  const worldNormal = cross(placed[bi].world, placed[bj].world);

  if (len(localSum) > 0.15 && len(worldSum) > 0.15 && len(localNormal) > 0.05 && len(worldNormal) > 0.05) {
    return alignPair(localSum, localNormal, worldSum, worldNormal);
  }
  return alignPair(placed[bi].local, placed[bj].local, placed[bi].world, placed[bj].world);
}

/**
 * Подбор поворота вокруг единственной известной связи — выбор конформации.
 * Перебираем 36 положений и берём то, при котором новые атомы дальше всего
 * от уже построенных. Так сам собой получается заторможенный этан.
 */
function chooseTwist(
  mol: Molecule,
  env: AtomEnvironment[],
  id: number,
  base: Mat3,
  axis: Vec3,
): Mat3 {
  const e = env[id];
  const me = mol.atoms[id];
  const pending = e.neighbors
    .map((nbId, slot) => ({ nbId, slot }))
    .filter(({ nbId }) => !mol.atoms[nbId].placed);
  if (pending.length === 0) return base;

  const others = mol.atoms.filter((a) => a.placed && a.id !== id);
  if (others.length <= 1) return base;

  let bestRot = base;
  let bestScore = Infinity;

  for (let step = 0; step < 36; step++) {
    const angle = (step * Math.PI * 2) / 36;
    const twist = rotationAxis(axis, angle);
    const rot = multiply(twist, base);

    let score = 0;
    for (const { nbId, slot } of pending) {
      const dir = norm(applyMat(rot, e.localDirs[slot]));
      const L = bondLength(me.el, mol.atoms[nbId].el, e.bondOrders[slot]);
      const p = add(me.pos, mul(dir, L));
      for (const other of others) {
        const d = Math.max(0.4, dist(p, other.pos));
        score += 1 / (d * d * d * d);
      }
    }
    if (score < bestScore) { bestScore = score; bestRot = rot; }
  }
  return bestRot;
}

function multiply(a: Mat3, b: Mat3): Mat3 {
  const r = new Array(9).fill(0) as number[];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let s = 0;
      for (let k = 0; k < 3; k++) s += a[i * 3 + k] * b[k * 3 + j];
      r[i * 3 + j] = s;
    }
  }
  return r as Mat3;
}

// ---------------------------------------------------------------------------
// Укладка циклов по шаблонам
// ---------------------------------------------------------------------------

function placeRings(
  mol: Molecule, env: AtomEnvironment[], rings: Ring[], templatePlaced: Set<number>,
): void {
  if (rings.length === 0) return;

  for (const group of fuseRingSystems(rings)) {
    const allPlanar = group.every((r) =>
      r.atoms.every((a) => env[a].steric <= 3 || mol.atoms[a].aromatic));

    // Мостиковые системы (адамантан, норборнан): есть атом, входящий сразу
    // в три цикла — ни плоским многоугольником, ни «креслом» такое не уложить.
    // Их полностью строит обход графа, а расхождения правит релаксация.
    if (hasBridgehead(group)) continue;

    if (allPlanar) {
      placePlanarSystem(mol, group, templatePlaced);
    } else {
      // насыщенная система: аккуратно кладём только первое кольцо,
      // остальное достроит обход графа
      placeSaturatedRing(mol, group[0], templatePlaced);
    }
  }
}

/** Есть ли атом, принадлежащий трём и более циклам группы. */
function hasBridgehead(group: Ring[]): boolean {
  const count = new Map<number, number>();
  for (const ring of group) {
    for (const a of ring.atoms) count.set(a, (count.get(a) ?? 0) + 1);
  }
  return [...count.values()].some((n) => n >= 3);
}

/** Средняя длина связи в кольце. */
function ringBondLength(mol: Molecule, ring: Ring): number {
  let sum = 0;
  for (let i = 0; i < ring.atoms.length; i++) {
    const a = ring.atoms[i];
    const b = ring.atoms[(i + 1) % ring.atoms.length];
    const bond = mol.bonds.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a));
    sum += bondLength(mol.atoms[a].el, mol.atoms[b].el, bond?.order ?? 1);
  }
  return sum / ring.atoms.length;
}

function polygonRadius(sides: number, side: number): number {
  return side / (2 * Math.sin(Math.PI / sides));
}

/** Плоская (ароматическая или сопряжённая) система колец — всё в плоскости z = 0. */
function placePlanarSystem(mol: Molecule, group: Ring[], templatePlaced: Set<number>): void {
  const first = group[0];
  const L = ringBondLength(mol, first);
  const R = polygonRadius(first.atoms.length, L);

  first.atoms.forEach((id, k) => {
    const a = (2 * Math.PI * k) / first.atoms.length + Math.PI / 2;
    mol.atoms[id].pos = v3(R * Math.cos(a), R * Math.sin(a), 0);
    mol.atoms[id].placed = true;
    templatePlaced.add(id);
  });

  const remaining = group.slice(1);
  let progress = true;
  while (remaining.length > 0 && progress) {
    progress = false;
    for (let idx = 0; idx < remaining.length; idx++) {
      const ring = remaining[idx];
      if (tryAttachPlanarRing(mol, ring, templatePlaced)) {
        remaining.splice(idx, 1);
        progress = true;
        break;
      }
    }
  }
  // непристроенные кольца достроит общий обход
}

function tryAttachPlanarRing(mol: Molecule, ring: Ring, templatePlaced: Set<number>): boolean {
  const size = ring.atoms.length;
  const placedIdx: number[] = [];
  ring.atoms.forEach((id, k) => { if (mol.atoms[id].placed) placedIdx.push(k); });
  if (placedIdx.length < 2) return false;

  // ищем два соседних в кольце уже размещённых атома — общую связь
  let anchorA = -1, anchorB = -1;
  for (const k of placedIdx) {
    const next = (k + 1) % size;
    if (mol.atoms[ring.atoms[next]].placed) { anchorA = k; anchorB = next; break; }
  }
  if (anchorA < 0) return false;

  const pA = mol.atoms[ring.atoms[anchorA]].pos;
  const pB = mol.atoms[ring.atoms[anchorB]].pos;
  const L = ringBondLength(mol, ring);
  const R = polygonRadius(size, L);
  const apothem = Math.sqrt(Math.max(0, R * R - (L / 2) * (L / 2)));

  const mid = mul(add(pA, pB), 0.5);
  const along = norm(sub(pB, pA));
  // нормаль в плоскости z = 0
  let outward = norm(v3(-along.y, along.x, 0));

  // направление «наружу» от уже построенной части
  const built = mol.atoms.filter((a) => a.placed).map((a) => a.pos);
  const bulk = centroid(built);
  if (dot(outward, sub(mid, bulk)) < 0) outward = mul(outward, -1);

  const center = add(mid, mul(outward, apothem));

  const angleOf = (p: Vec3): number => Math.atan2(p.y - center.y, p.x - center.x);
  const aAngle = angleOf(pA);
  const bAngle = angleOf(pB);
  const stepSize = (2 * Math.PI) / size;
  let delta = bAngle - aAngle;
  while (delta > Math.PI) delta -= 2 * Math.PI;
  while (delta < -Math.PI) delta += 2 * Math.PI;
  const direction = delta > 0 ? 1 : -1;

  for (let step = 1; step < size; step++) {
    const ringIndex = (anchorB + step) % size;
    const id = ring.atoms[ringIndex];
    if (mol.atoms[id].placed) continue;
    const a = bAngle + direction * stepSize * step;
    mol.atoms[id].pos = v3(center.x + R * Math.cos(a), center.y + R * Math.sin(a), 0);
    mol.atoms[id].placed = true;
    templatePlaced.add(id);
  }
  return true;
}

/** Насыщенный цикл: «кресло» для шестичленного, «конверт» для пятичленного. */
function placeSaturatedRing(mol: Molecule, ring: Ring, templatePlaced: Set<number>): void {
  const size = ring.atoms.length;
  const L = ringBondLength(mol, ring);

  if (size === 6) {
    // Классическое кресло: радиус 0,948·L, отклонение по оси ±0,162·L.
    // При L = 1,54 Å это даёт длину связи 1,54 Å и валентный угол 111,5°.
    const R = 0.948 * L;
    const h = 0.162 * L;
    ring.atoms.forEach((id, k) => {
      const a = (2 * Math.PI * k) / 6;
      mol.atoms[id].pos = v3(R * Math.cos(a), R * Math.sin(a), k % 2 === 0 ? h : -h);
      mol.atoms[id].placed = true;
      templatePlaced.add(id);
    });
  } else if (size === 5) {
    // Конверт: четыре атома в плоскости, пятый приподнят.
    const R = polygonRadius(5, L);
    ring.atoms.forEach((id, k) => {
      const a = (2 * Math.PI * k) / 5;
      mol.atoms[id].pos = v3(R * Math.cos(a), R * Math.sin(a), k === 4 ? 0.28 * L : 0);
      mol.atoms[id].placed = true;
      templatePlaced.add(id);
    });
  } else {
    const R = polygonRadius(size, L);
    ring.atoms.forEach((id, k) => {
      const a = (2 * Math.PI * k) / size;
      // лёгкая складка для четырёхчленного цикла
      const z = size === 4 ? (k % 2 === 0 ? 0.05 * L : -0.05 * L) : 0;
      mol.atoms[id].pos = v3(R * Math.cos(a), R * Math.sin(a), z);
      mol.atoms[id].placed = true;
      templatePlaced.add(id);
    });
  }
}

// ---------------------------------------------------------------------------
// Релаксация: разводим атомы, которые оказались слишком близко
// ---------------------------------------------------------------------------

function relax(mol: Molecule, env: AtomEnvironment[], frozen: Set<number>): void {
  const n = mol.atoms.length;
  if (n < 4) return;

  const targetBond = new Map<string, number>();
  for (const b of mol.bonds) {
    targetBond.set(key(b.a, b.b), bondLength(mol.atoms[b.a].el, mol.atoms[b.b].el, b.order));
  }

  // целевые расстояния 1–3 (через один атом) — они и держат валентные углы
  const target13: { i: number; j: number; d: number }[] = [];
  for (let id = 0; id < n; id++) {
    const e = env[id];
    for (let a = 0; a < e.neighbors.length; a++) {
      for (let b = a + 1; b < e.neighbors.length; b++) {
        const i = e.neighbors[a], j = e.neighbors[b];
        // угол внутри уложенного по шаблону кольца задан геометрией самого кольца
        if (frozen.has(id) && frozen.has(i) && frozen.has(j)) continue;
        const angle = Math.acos(Math.max(-1, Math.min(1, dot(e.localDirs[a], e.localDirs[b]))));
        const la = bondLength(mol.atoms[id].el, mol.atoms[i].el, e.bondOrders[a]);
        const lb = bondLength(mol.atoms[id].el, mol.atoms[j].el, e.bondOrders[b]);
        target13.push({ i, j, d: Math.sqrt(la * la + lb * lb - 2 * la * lb * Math.cos(angle)) });
      }
    }
  }

  const bondedWithin3 = buildProximity(mol, 3);
  const movable = mol.atoms.map((a) => !frozen.has(a.id));

  for (let iter = 0; iter < 600; iter++) {
    const shift: Vec3[] = mol.atoms.map(() => v3());

    const applyPair = (i: number, j: number, targetDist: number, k: number): void => {
      const d = sub(mol.atoms[j].pos, mol.atoms[i].pos);
      const l = len(d);
      if (l < 1e-6) return;
      const correction = mul(d, ((l - targetDist) / l) * k * 0.5);
      shift[i] = add(shift[i], correction);
      shift[j] = sub(shift[j], correction);
    };

    for (const b of mol.bonds) applyPair(b.a, b.b, targetBond.get(key(b.a, b.b))!, 0.85);
    for (const t of target13) applyPair(t.i, t.j, t.d, 0.35);

    // мягкое отталкивание несвязанных атомов
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (bondedWithin3.has(key(i, j))) continue;
        const d = dist(mol.atoms[i].pos, mol.atoms[j].pos);
        const minDist = 2.4;
        if (d < minDist && d > 1e-6) applyPair(i, j, minDist, 0.25);
      }
    }

    for (let i = 0; i < n; i++) {
      if (!movable[i]) continue;
      mol.atoms[i].pos = add(mol.atoms[i].pos, mul(shift[i], 0.6));
    }
  }
}

const key = (a: number, b: number): string => (a < b ? `${a}-${b}` : `${b}-${a}`);

function buildProximity(mol: Molecule, maxBonds: number): Set<string> {
  const n = mol.atoms.length;
  const nb: number[][] = mol.atoms.map(() => []);
  for (const b of mol.bonds) { nb[b.a].push(b.b); nb[b.b].push(b.a); }

  const near = new Set<string>();
  for (let start = 0; start < n; start++) {
    const depth = new Map<number, number>([[start, 0]]);
    const queue = [start];
    while (queue.length > 0) {
      const cur = queue.shift()!;
      const d = depth.get(cur)!;
      if (d >= maxBonds) continue;
      for (const next of nb[cur]) {
        if (depth.has(next)) continue;
        depth.set(next, d + 1);
        near.add(key(start, next));
        queue.push(next);
      }
    }
  }
  return near;
}

// ---------------------------------------------------------------------------
// Разворот молекулы по главным осям
// ---------------------------------------------------------------------------

function orientByPrincipalAxes(mol: Molecule): void {
  const heavy = mol.atoms.filter((a) => a.el !== 'H');
  const reference = heavy.length >= 2 ? heavy : mol.atoms;
  const center = centroid(reference.map((a) => a.pos));

  for (const a of mol.atoms) a.pos = sub(a.pos, center);

  if (mol.atoms.length < 3) return;

  // матрица ковариации положений
  let xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
  for (const a of reference) {
    const p = a.pos;
    xx += p.x * p.x; xy += p.x * p.y; xz += p.x * p.z;
    yy += p.y * p.y; yz += p.y * p.z; zz += p.z * p.z;
  }
  const { values, vectors } = jacobiEigen([xx, xy, xz, xy, yy, yz, xz, yz, zz]);

  const order = [0, 1, 2].sort((a, b) => values[b] - values[a]);
  let ex = norm(vectors[order[0]]);
  let ey = norm(vectors[order[1]]);
  let ez = cross(ex, ey);
  ey = cross(ez, ex);

  for (const a of mol.atoms) {
    const p = a.pos;
    a.pos = v3(dot(p, ex), dot(p, ey), dot(p, ez));
  }

  // небольшой наклон, чтобы объём читался лучше, чем строго анфас
  const tilt = rotationAxis(v3(1, 0, 0), -0.32);
  const spin = rotationAxis(v3(0, 1, 0), 0.42);
  for (const a of mol.atoms) a.pos = applyMat(spin, applyMat(tilt, a.pos));
}

export function cloneMolecule(mol: Molecule): Molecule {
  return {
    atoms: mol.atoms.map((a) => ({ ...a, pos: clone(a.pos) })),
    bonds: mol.bonds.map((b) => ({ ...b })),
    meta: mol.meta,
  };
}
