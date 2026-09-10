import type { Molecule } from './types';
import { neighborLists } from './smiles';

/**
 * ДЕЛОКАЛИЗАЦИЯ (РЕЗОНАНС)
 *
 * Формула Льюиса для нитрат-иона рисуется с одной двойной и двумя одинарными
 * связями, но на самом деле все три связи N-O одинаковы: длина 1,24 Å, порядок 1⅓.
 * Реальная частица — не одна из резонансных структур, а их наложение.
 *
 * Программа обнаруживает такие случаи сама: если центральный атом связан
 * с несколькими КОНЦЕВЫМИ атомами одного элемента, но связи получились разной
 * кратности, значит выбор «где рисовать двойную связь» произволен, и связи
 * надо усреднить. Заряд усредняется вместе с кратностью.
 *
 * Правило срабатывает для нитрат-, нитрит-, карбонат-, сульфат-, фосфат- и
 * ацетат-ионов, для нитрогруппы, для озона, для азотной кислоты — и НЕ
 * срабатывает для уксусной кислоты и сложных эфиров, где второй кислород
 * несёт водород или углеродный заместитель и потому неравноценен первому.
 */
export function applyResonance(mol: Molecule): void {
  const nb = neighborLists(mol);

  for (const center of mol.atoms) {
    if (nb[center.id].length < 2) continue;

    // концевые соседи, сгруппированные по элементу
    const groups = new Map<string, number[]>();
    for (const nbId of nb[center.id]) {
      const neighbor = mol.atoms[nbId];
      if (neighbor.el === 'H') continue;
      if (nb[nbId].length !== 1) continue; // не концевой — не равноценен
      if (!groups.has(neighbor.el)) groups.set(neighbor.el, []);
      groups.get(neighbor.el)!.push(nbId);
    }

    for (const [, members] of groups) {
      if (members.length < 2) continue;

      const bonds = members.map((id) => findBond(mol, center.id, id)!);
      if (bonds.some((b) => b.aromatic)) continue;

      const orders = bonds.map((b) => b.order);
      const charges = members.map((id) => mol.atoms[id].charge);
      const sameOrder = orders.every((o) => Math.abs(o - orders[0]) < 1e-9);
      const sameCharge = charges.every((c) => c === charges[0]);
      if (sameOrder && sameCharge) continue; // равноценны и без усреднения

      const avgOrder = orders.reduce((s, o) => s + o, 0) / orders.length;
      const avgCharge = charges.reduce((s, c) => s + c, 0) / charges.length;

      bonds.forEach((b) => { b.order = avgOrder; b.delocalized = true; });
      members.forEach((id) => { mol.atoms[id].charge = avgCharge; });
    }
  }
}

function findBond(mol: Molecule, a: number, b: number) {
  return mol.bonds.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a));
}

/** Красивая запись дробного порядка связи: 1⅓, 1½, 1⅔. */
export function formatBondOrder(order: number): string {
  const table: [number, string][] = [
    [1, '1'], [4 / 3, '1⅓'], [1.25, '1¼'], [1.5, '1½'], [5 / 3, '1⅔'],
    [2, '2'], [2.5, '2½'], [3, '3'],
  ];
  for (const [value, label] of table) {
    if (Math.abs(order - value) < 0.02) return label;
  }
  return order.toFixed(2).replace(/0$/, '');
}
