import type { Molecule } from './types';
import { element, isKnownElement } from './periodic';

/**
 * Брутто-формула в системе Хилла: сначала углерод, затем водород,
 * затем остальные элементы по алфавиту. Если углерода нет — всё по алфавиту.
 */
export function composition(mol: Molecule): Map<string, number> {
  const counts = new Map<string, number>();
  for (const atom of mol.atoms) {
    counts.set(atom.el, (counts.get(atom.el) ?? 0) + 1);
  }
  return counts;
}

export function totalCharge(mol: Molecule): number {
  return mol.atoms.reduce((s, a) => s + a.charge, 0);
}

export function molecularMass(mol: Molecule): number {
  return mol.atoms.reduce((s, a) => s + element(a.el).mass, 0);
}

function hillOrder(counts: Map<string, number>): string[] {
  const symbols = [...counts.keys()];
  const rest = symbols.filter((s) => s !== 'C' && s !== 'H').sort((a, b) => a.localeCompare(b));
  if (counts.has('C')) {
    return ['C', ...(counts.has('H') ? ['H'] : []), ...rest];
  }
  return symbols.sort((a, b) => a.localeCompare(b));
}

const SUB = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'];
const SUP: Record<string, string> = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '+': '⁺', '-': '⁻' };

const toSub = (n: number): string => String(n).split('').map((d) => SUB[+d]).join('');
const toSup = (s: string): string => s.split('').map((c) => SUP[c] ?? c).join('');

/** Формула обычными символами: "C2H6O". Используется для поиска. */
export function plainFormula(mol: Molecule): string {
  const counts = composition(mol);
  return hillOrder(counts).map((s) => s + (counts.get(s)! > 1 ? counts.get(s) : '')).join('');
}

/** Формула с настоящими подстрочными цифрами и зарядом: "C₂H₆O", "NO₃⁻". */
export function prettyFormula(mol: Molecule): string {
  const counts = composition(mol);
  let s = hillOrder(counts).map((sym) => sym + (counts.get(sym)! > 1 ? toSub(counts.get(sym)!) : '')).join('');
  const q = totalCharge(mol);
  if (q !== 0) {
    const magnitude = Math.abs(q) > 1 ? String(Math.abs(q)) : '';
    s += toSup(magnitude + (q > 0 ? '+' : '-'));
  }
  return s;
}

/**
 * Разбор формулы, введённой пользователем: "C2H6O", "h2o", "CH3COOH".
 * Возвращает нормализованный ключ состава либо null, если это не формула.
 */
export function parseFormulaQuery(input: string): string | null {
  const text = input.replace(/[\s()·.]/g, '');
  if (!text) return null;

  const counts = new Map<string, number>();
  let i = 0;
  while (i < text.length) {
    // элемент: заглавная + необязательная строчная
    let symbol: string | null = null;
    const two = text.slice(i, i + 2);
    if (two.length === 2 && /^[A-Za-z][a-z]$/.test(two)) {
      const candidate = two[0].toUpperCase() + two[1].toLowerCase();
      if (isKnownElement(candidate)) { symbol = candidate; i += 2; }
    }
    if (symbol === null) {
      const one = text[i].toUpperCase();
      if (!isKnownElement(one)) return null;
      symbol = one;
      i += 1;
    }
    let digits = '';
    while (i < text.length && /[0-9]/.test(text[i])) digits += text[i++];
    const n = digits === '' ? 1 : parseInt(digits, 10);
    counts.set(symbol, (counts.get(symbol) ?? 0) + n);
  }
  if (counts.size === 0) return null;
  return compositionKey(counts);
}

/** Канонический ключ состава для сравнения формул. */
export function compositionKey(counts: Map<string, number>): string {
  return [...counts.entries()]
    .filter(([, n]) => n > 0)
    .sort((a, b) => a[0].localeCompare(b[0]))
    .map(([s, n]) => `${s}${n}`)
    .join('');
}

export function moleculeCompositionKey(mol: Molecule): string {
  return compositionKey(composition(mol));
}

/** Нормализация подписи валентного угла: внешние атомы сортируются по алфавиту. */
export function angleLabel(outerA: string, center: string, outerB: string): string {
  const [x, y] = [outerA, outerB].sort((a, b) => a.localeCompare(b));
  return `${x}-${center}-${y}`;
}
