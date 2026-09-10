import elementsData from '../data/elements.json';

export interface ElementInfo {
  z: number;
  name: string;
  group: number;
  period: number;
  /** Число валентных электронов. */
  ve: number;
  /** Электроотрицательность по Полингу. */
  en: number;
  /** Ковалентный радиус одинарной связи, Å. */
  r: number;
  vdw: number;
  color: string;
  mass: number;
}

const TABLE = elementsData.elements as Record<string, ElementInfo>;

const UNKNOWN: ElementInfo = {
  z: 0, name: 'неизвестный', group: 14, period: 3,
  ve: 4, en: 2.2, r: 1.0, vdw: 1.7, color: '#9aa0aa', mass: 12,
};

export function element(symbol: string): ElementInfo {
  return TABLE[symbol] ?? UNKNOWN;
}

export function isKnownElement(symbol: string): boolean {
  return symbol in TABLE;
}

/**
 * Стандартные валентности для расчёта числа неявных атомов водорода.
 * Значения соответствуют «органическому подмножеству» нотации SMILES.
 */
const VALENCES: Record<string, number[]> = {
  H: [1], B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6],
  F: [1], Cl: [1], Br: [1], I: [1],
  Si: [4], Ge: [4], Sn: [4], As: [3, 5], Se: [2, 4, 6], Te: [2, 4, 6], Sb: [3, 5],
  Be: [2], Mg: [2], Al: [3], Li: [1], Na: [1], K: [1], Ca: [2],
  Xe: [0], Kr: [0], Ar: [0], Ne: [0], He: [0],
};

/** Элементы, у которых положительный заряд УМЕНЬШАЕТ валентность (левая часть таблицы). */
const ELECTRON_DEFICIENT = new Set(['B', 'C', 'Si', 'Ge', 'Sn', 'Al', 'Be', 'Mg']);

/**
 * Число неявных атомов водорода: берём наименьшую стандартную валентность,
 * которой хватает на уже имеющиеся связи, и добираем остаток водородом.
 *
 * Для АРОМАТИЧЕСКОГО атома действует особое правило. Кратности связей в кольце
 * не определены до расстановки формы Кекуле, поэтому считаем так: каждая связь
 * кольца даёт единицу, плюс единица за участие в π-системе. И используется
 * только НИЗШАЯ валентность: у нейтрального азота не может быть четырёх связей,
 * поэтому в кофеине азот с метильной группой водорода не получает.
 *
 *   бензол   c:  σ=2, π=1 → 4 − 3 = 1 атом H
 *   толуол   c:  σ=3, π=1 → 4 − 4 = 0
 *   пиридин  n:  σ=2, π=1 → 3 − 3 = 0
 *   кофеин   n:  σ=3, π=1 → 3 − 4 < 0 → 0
 *   фуран    o:  σ=2, π=1 → 2 − 3 < 0 → 0
 */
export function implicitHydrogens(
  symbol: string, bondSum: number, charge: number, aromatic = false,
): number {
  const valences = VALENCES[symbol];
  if (!valences || valences.length === 0 || valences[0] === 0) return 0;

  const shift = ELECTRON_DEFICIENT.has(symbol) ? -Math.abs(charge) : charge;

  if (aromatic) {
    return Math.max(0, Math.round(valences[0] + shift - bondSum));
  }

  let target = valences[valences.length - 1] + shift;
  for (const v of valences) {
    const adjusted = v + shift;
    if (adjusted >= bondSum - 1e-6) { target = adjusted; break; }
  }
  return Math.max(0, Math.round(target - bondSum));
}

/**
 * Длина связи как сумма ковалентных радиусов с поправкой на кратность.
 * Кратные связи короче одинарных примерно на 0,1 и 0,17 Å.
 */
export function bondLength(elA: string, elB: string, order: number): number {
  const base = element(elA).r + element(elB).r;
  return base - bondShortening(order);
}

/**
 * Укорочение связи при увеличении кратности, Å. Зависимость непрерывная,
 * поэтому корректно обрабатываются и дробные порядки: ароматическая связь
 * 1,5 и делокализованная 1⅓ в нитрат-ионе.
 */
export function bondShortening(order: number): number {
  const o = Math.max(1, Math.min(3, order));
  return o <= 2 ? 0.19 * (o - 1) : 0.19 + 0.13 * (o - 2);
}
