import type { Molecule } from './types';
import { element, isKnownElement } from './periodic';

/**
 * Состав молекулы: сколько атомов каждого элемента.
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

/**
 * Электроотрицательностный ряд элементов (IUPAC, Red Book, табл. VI).
 * В формуле неорганического вещества первым пишется элемент, стоящий в ряду
 * левее (более электроположительный), поэтому получается SO₂, PCl₃, XeF₂,
 * а кислород и фтор оказываются в конце. Водород стоит между азотом и серой:
 * отсюда NH₃, PH₃, SiH₄, но H₂O, H₂S, HCl.
 */
const EN_ORDER = [
  'Rn', 'Xe', 'Kr', 'Ar', 'Ne', 'He',
  'Fr', 'Cs', 'Rb', 'K', 'Na', 'Li',
  'Ra', 'Ba', 'Sr', 'Ca', 'Mg', 'Be',
  'Tl', 'In', 'Ga', 'Al', 'B',
  'Pb', 'Sn', 'Ge', 'Si', 'C',
  'Bi', 'Sb', 'As', 'P', 'N',
  'H',
  'Po', 'Te', 'Se', 'S',
  'At', 'I', 'Br', 'Cl',
  'O', 'F',
];

const enRank = new Map(EN_ORDER.map((s, i) => [s, i]));

/** Элементы вне ряда (металлы середины таблицы) — по возрастанию ЭО, как и весь ряд. */
const rankOf = (symbol: string): number => enRank.get(symbol) ?? -1;

/**
 * Вещества, которые по традиции записывают вопреки общему правилу:
 * гидроксид-ион OH⁻ (а не HO⁻), циановодород HCN (а не CHN по Хиллу).
 * Ключ — канонический ключ состава из compositionKey.
 */
const CONVENTIONAL_ORDER: Record<string, string[]> = {
  H1O1: ['O', 'H'],
  C1H1N1: ['H', 'C', 'N'],
};

/**
 * Порядок элементов в брутто-формуле.
 *
 * Для органики — система Хилла (C, H, остальное по алфавиту): C₂H₆O, C₆H₆.
 * Для неорганики — электроотрицательностный ряд: SO₂, PCl₃, NH₃, а не
 * алфавитные O₂S, Cl₃P, H₃N. У кислородсодержащих кислот водород по традиции
 * выносится вперёд: H₂SO₄, HNO₃, H₃PO₄.
 */
function formulaOrder(counts: Map<string, number>): string[] {
  const conventional = CONVENTIONAL_ORDER[compositionKey(counts)];
  if (conventional) return conventional;

  const symbols = [...counts.keys()];

  if (counts.has('C')) {
    const rest = symbols.filter((s) => s !== 'C' && s !== 'H').sort((a, b) => a.localeCompare(b));
    return ['C', ...(counts.has('H') ? ['H'] : []), ...rest];
  }

  const ordered = symbols.sort((a, b) => rankOf(a) - rankOf(b) || a.localeCompare(b));

  // кислота: есть и водород, и кислород, и ещё хотя бы один элемент
  if (counts.size > 2 && counts.has('H') && counts.has('O')) {
    return ['H', ...ordered.filter((s) => s !== 'H')];
  }
  return ordered;
}

const SUB = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'];
const SUP: Record<string, string> = { '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '+': '⁺', '-': '⁻' };

const toSub = (n: number): string => String(n).split('').map((d) => SUB[+d]).join('');
const toSup = (s: string): string => s.split('').map((c) => SUP[c] ?? c).join('');

/** Формула обычными символами: "C2H6O", "H2SO4". Используется для поиска. */
export function plainFormula(mol: Molecule): string {
  const counts = composition(mol);
  return formulaOrder(counts).map((s) => s + (counts.get(s)! > 1 ? counts.get(s) : '')).join('');
}

/** Формула с настоящими подстрочными цифрами и зарядом: "C₂H₆O", "NO₃⁻". */
export function prettyFormula(mol: Molecule): string {
  const counts = composition(mol);
  let s = formulaOrder(counts).map((sym) => sym + (counts.get(sym)! > 1 ? toSub(counts.get(sym)!) : '')).join('');
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
