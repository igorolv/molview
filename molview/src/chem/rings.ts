import type { Molecule, Ring } from './types';
import { neighborLists } from './smiles';

/**
 * Поиск наименьшего набора наименьших циклов (SSSR).
 *
 * Приём простой и надёжный: для каждой связи мысленно её удаляем и ищем
 * кратчайший путь между её концами. Путь плюс сама связь и есть цикл,
 * причём наименьший из проходящих через эту связь. Дальше отбираем
 * независимые циклы — столько, сколько даёт цикломатическое число графа.
 */
export function findRings(mol: Molecule): Ring[] {
  const nb = neighborLists(mol);
  const heavy = mol.atoms.filter((a) => a.el !== 'H').map((a) => a.id);
  const isHeavy = new Set(heavy);

  const candidates: number[][] = [];

  for (const bond of mol.bonds) {
    if (!isHeavy.has(bond.a) || !isHeavy.has(bond.b)) continue;
    const path = shortestPath(nb, isHeavy, bond.a, bond.b, bond);
    if (path && path.length >= 3) candidates.push(path);
  }

  // уникальные циклы по набору атомов, от меньших к большим
  const seen = new Set<string>();
  const unique: number[][] = [];
  for (const c of candidates.sort((x, y) => x.length - y.length)) {
    const key = [...c].sort((a, b) => a - b).join(',');
    if (seen.has(key)) continue;
    seen.add(key);
    unique.push(c);
  }

  // цикломатическое число: сколько независимых циклов есть в графе
  const heavyBonds = mol.bonds.filter((b) => isHeavy.has(b.a) && isHeavy.has(b.b));
  const components = countComponents(nb, heavy, isHeavy);
  const cycleCount = heavyBonds.length - heavy.length + components;

  const chosen: number[][] = [];
  const coveredBonds = new Set<string>();
  for (const c of unique) {
    if (chosen.length >= cycleCount) break;
    const bondKeys = ringBondKeys(c);
    if (bondKeys.some((k) => !coveredBonds.has(k)) || chosen.length === 0) {
      chosen.push(c);
      bondKeys.forEach((k) => coveredBonds.add(k));
    }
  }

  return chosen.map((atoms) => describeRing(mol, atoms));
}

function ringBondKeys(ring: number[]): string[] {
  const keys: string[] = [];
  for (let i = 0; i < ring.length; i++) {
    const a = ring[i];
    const b = ring[(i + 1) % ring.length];
    keys.push(a < b ? `${a}-${b}` : `${b}-${a}`);
  }
  return keys;
}

function shortestPath(
  nb: number[][],
  isHeavy: Set<number>,
  from: number,
  to: number,
  forbidden: { a: number; b: number },
): number[] | null {
  const prev = new Map<number, number>();
  const queue = [from];
  const visited = new Set<number>([from]);

  while (queue.length > 0) {
    const cur = queue.shift()!;
    if (cur === to) break;
    for (const next of nb[cur]) {
      if (!isHeavy.has(next)) continue;
      const isForbidden =
        (cur === forbidden.a && next === forbidden.b) || (cur === forbidden.b && next === forbidden.a);
      if (isForbidden) continue;
      if (visited.has(next)) continue;
      visited.add(next);
      prev.set(next, cur);
      queue.push(next);
    }
  }

  if (!visited.has(to)) return null;
  const path: number[] = [];
  let cur: number | undefined = to;
  while (cur !== undefined) {
    path.push(cur);
    if (cur === from) break;
    cur = prev.get(cur);
  }
  return path.length >= 3 ? path.reverse() : null;
}

function countComponents(nb: number[][], heavy: number[], isHeavy: Set<number>): number {
  const visited = new Set<number>();
  let count = 0;
  for (const start of heavy) {
    if (visited.has(start)) continue;
    count++;
    const stack = [start];
    visited.add(start);
    while (stack.length > 0) {
      const cur = stack.pop()!;
      for (const next of nb[cur]) {
        if (!isHeavy.has(next) || visited.has(next)) continue;
        visited.add(next);
        stack.push(next);
      }
    }
  }
  return count;
}

function describeRing(mol: Molecule, atoms: number[]): Ring {
  let aromatic = true;
  for (let i = 0; i < atoms.length; i++) {
    const a = atoms[i];
    const b = atoms[(i + 1) % atoms.length];
    const bond = mol.bonds.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a));
    if (!bond || !bond.aromatic) { aromatic = false; break; }
  }
  return { atoms, aromatic, planar: aromatic };
}

/** Группировка циклов в конденсированные системы (имеющие общие атомы). */
export function fuseRingSystems(rings: Ring[]): Ring[][] {
  const groups: Ring[][] = [];
  const assigned = new Array(rings.length).fill(false);

  for (let i = 0; i < rings.length; i++) {
    if (assigned[i]) continue;
    const group = [rings[i]];
    assigned[i] = true;
    let grew = true;
    while (grew) {
      grew = false;
      const inGroup = new Set<number>();
      group.forEach((r) => r.atoms.forEach((a) => inGroup.add(a)));
      for (let j = 0; j < rings.length; j++) {
        if (assigned[j]) continue;
        if (rings[j].atoms.some((a) => inGroup.has(a))) {
          group.push(rings[j]);
          assigned[j] = true;
          grew = true;
        }
      }
    }
    groups.push(group);
  }
  return groups;
}
