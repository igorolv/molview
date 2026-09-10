import type { Atom, Bond, Molecule } from './types';
import { v3 } from './vec';
import { implicitHydrogens, isKnownElement } from './periodic';

/**
 * Разбор строки SMILES в граф молекулы.
 *
 * Поддерживается: органическое подмножество (B, C, N, O, P, S, F, Cl, Br, I),
 * атомы в квадратных скобках с зарядом и явным числом водородов,
 * кратные связи - = # :, ветвления, замыкания циклов (в том числе %nn),
 * ароматические атомы в нижнем регистре, разрыв '.'.
 *
 * Не поддерживается (сознательно, для школьного проекта не нужно):
 * стереохимия @/@@ и /\ — эти символы просто пропускаются.
 */

const ORGANIC_SUBSET = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']);
const AROMATIC_SYMBOLS = new Set(['b', 'c', 'n', 'o', 'p', 's']);

export class SmilesError extends Error {
  constructor(message: string, readonly position: number) {
    super(message);
    this.name = 'SmilesError';
  }
}

interface RawAtom {
  el: string;
  charge: number;
  aromatic: boolean;
  hCount: number;
  explicitH: boolean;
}

interface OpenRing {
  atom: number;
  order: number | null;
  aromatic: boolean;
  position: number;
}

export function parseSmiles(input: string): Molecule {
  const src = input.trim();
  if (!src) throw new SmilesError('Пустая строка', 0);

  const raw: RawAtom[] = [];
  const bonds: Bond[] = [];
  const rings = new Map<number, OpenRing>();
  const branchStack: number[] = [];

  let i = 0;
  let prev: number | null = null;
  let pendingOrder: number | null = null;
  let pendingAromatic = false;

  const fail: (msg: string) => never = (msg) => { throw new SmilesError(msg, i); };

  const linkTo = (target: number, order: number | null, aromatic: boolean, ringClosure: boolean) => {
    if (prev === null) return;
    const bothAromatic = raw[prev].aromatic && raw[target].aromatic;
    let finalOrder: number;
    let finalAromatic: boolean;
    if (order === null) {
      finalAromatic = bothAromatic;
      finalOrder = bothAromatic ? 1.5 : 1;
    } else {
      finalAromatic = aromatic;
      finalOrder = order;
    }
    bonds.push({ a: prev, b: target, order: finalOrder, aromatic: finalAromatic, ringClosure, delocalized: false });
  };

  const pushAtom = (a: RawAtom): number => {
    const id = raw.length;
    raw.push(a);
    if (prev !== null) {
      const bothAromatic = raw[prev].aromatic && a.aromatic;
      const order = pendingOrder ?? (bothAromatic ? 1.5 : 1);
      const aromatic = pendingOrder === null ? bothAromatic : pendingAromatic;
      bonds.push({ a: prev, b: id, order, aromatic, ringClosure: false, delocalized: false });
    }
    pendingOrder = null;
    pendingAromatic = false;
    prev = id;
    return id;
  };

  while (i < src.length) {
    const ch = src[i];

    // --- ветвления ---------------------------------------------------------
    if (ch === '(') {
      if (prev === null) fail('Ветвление «(» не может открывать формулу');
      branchStack.push(prev);
      i++;
      continue;
    }
    if (ch === ')') {
      const back = branchStack.pop();
      if (back === undefined) fail('Лишняя закрывающая скобка «)»');
      prev = back;
      i++;
      continue;
    }

    // --- символы связей ----------------------------------------------------
    if (ch === '-') { pendingOrder = 1; pendingAromatic = false; i++; continue; }
    if (ch === '=') { pendingOrder = 2; pendingAromatic = false; i++; continue; }
    if (ch === '#') { pendingOrder = 3; pendingAromatic = false; i++; continue; }
    if (ch === ':') { pendingOrder = 1.5; pendingAromatic = true; i++; continue; }
    if (ch === '/' || ch === '\\') { i++; continue; } // стереохимия двойной связи — игнорируем
    if (ch === '~') { pendingOrder = 1; i++; continue; }

    // --- разрыв цепи -------------------------------------------------------
    if (ch === '.') { prev = null; pendingOrder = null; i++; continue; }

    // --- замыкания циклов --------------------------------------------------
    if (ch === '%' || /[0-9]/.test(ch)) {
      let label: number;
      if (ch === '%') {
        const digits = src.slice(i + 1, i + 3);
        if (!/^[0-9]{2}$/.test(digits)) fail('После «%» должны идти две цифры номера цикла');
        label = parseInt(digits, 10);
        i += 3;
      } else {
        label = parseInt(ch, 10);
        i += 1;
      }
      if (prev === null) fail('Номер цикла не может стоять до первого атома');

      const open = rings.get(label);
      if (open === undefined) {
        rings.set(label, { atom: prev, order: pendingOrder, aromatic: pendingAromatic, position: i });
      } else {
        if (open.atom === prev) fail(`Цикл ${label} замыкается сам на себя`);
        const order = pendingOrder ?? open.order;
        const aromatic = pendingOrder !== null ? pendingAromatic : open.aromatic;
        linkTo(open.atom, order, aromatic, true);
        rings.delete(label);
      }
      pendingOrder = null;
      pendingAromatic = false;
      continue;
    }

    // --- атом в квадратных скобках ----------------------------------------
    if (ch === '[') {
      const close = src.indexOf(']', i);
      if (close < 0) fail('Не закрыта квадратная скобка «[»');
      const body = src.slice(i + 1, close);
      const atom = parseBracketAtom(body, i);
      i = close + 1;
      pushAtom(atom);
      continue;
    }

    // --- атом органического подмножества ----------------------------------
    const two = src.slice(i, i + 2);
    if (ORGANIC_SUBSET.has(two)) {
      i += 2;
      pushAtom({ el: two, charge: 0, aromatic: false, hCount: 0, explicitH: false });
      continue;
    }
    if (ORGANIC_SUBSET.has(ch)) {
      i += 1;
      pushAtom({ el: ch, charge: 0, aromatic: false, hCount: 0, explicitH: false });
      continue;
    }
    if (AROMATIC_SYMBOLS.has(ch)) {
      i += 1;
      pushAtom({ el: ch.toUpperCase(), charge: 0, aromatic: true, hCount: 0, explicitH: false });
      continue;
    }
    if (ch === '*') {
      i += 1;
      pushAtom({ el: 'C', charge: 0, aromatic: false, hCount: 0, explicitH: false });
      continue;
    }

    fail(`Непонятный символ «${ch}» в позиции ${i + 1}`);
  }

  if (branchStack.length > 0) fail('Не закрыта скобка «(»');
  if (rings.size > 0) {
    const label = [...rings.keys()][0];
    throw new SmilesError(`Цикл с номером ${label} открыт, но не замкнут`, src.length);
  }
  if (raw.length === 0) throw new SmilesError('Не найдено ни одного атома', 0);

  return buildMolecule(raw, bonds);
}

function parseBracketAtom(body: string, position: number): RawAtom {
  let p = 0;
  // изотоп — пропускаем
  while (p < body.length && /[0-9]/.test(body[p])) p++;

  let el = '';
  let aromatic = false;
  if (p < body.length && /[a-z]/.test(body[p]) && AROMATIC_SYMBOLS.has(body[p])) {
    aromatic = true;
    el = body[p].toUpperCase();
    p++;
  } else if (p < body.length && /[A-Z]/.test(body[p])) {
    el = body[p];
    p++;
    if (p < body.length && /[a-z]/.test(body[p]) && isKnownElement(el + body[p])) {
      el += body[p];
      p++;
    }
  } else {
    throw new SmilesError(`Не удалось определить элемент в «[${body}]»`, position);
  }

  // хиральность — пропускаем
  while (p < body.length && body[p] === '@') p++;

  let hCount = 0;
  if (p < body.length && body[p] === 'H') {
    p++;
    let digits = '';
    while (p < body.length && /[0-9]/.test(body[p])) digits += body[p++];
    hCount = digits === '' ? 1 : parseInt(digits, 10);
  }

  let charge = 0;
  while (p < body.length && (body[p] === '+' || body[p] === '-')) {
    const sign = body[p] === '+' ? 1 : -1;
    p++;
    let digits = '';
    while (p < body.length && /[0-9]/.test(body[p])) digits += body[p++];
    if (digits) {
      charge += sign * parseInt(digits, 10);
    } else {
      charge += sign;
      while (p < body.length && body[p] === (sign > 0 ? '+' : '-')) { charge += sign; p++; }
    }
  }

  // класс атома «:12» — пропускаем
  if (p < body.length && body[p] === ':') {
    p++;
    while (p < body.length && /[0-9]/.test(body[p])) p++;
  }

  if (p !== body.length) {
    throw new SmilesError(`Лишние символы в «[${body}]»`, position);
  }
  if (!isKnownElement(el)) {
    throw new SmilesError(`Элемент «${el}» отсутствует в справочнике программы`, position);
  }

  return { el, charge, aromatic, hCount, explicitH: true };
}

/** Развёртывание неявных водородов в явные атомы. */
function buildMolecule(raw: RawAtom[], bonds: Bond[]): Molecule {
  const bondSum = new Array(raw.length).fill(0);
  // для ароматических атомов: число связей плюс «лишняя» кратность
  // неароматических связей — см. пояснение в implicitHydrogens
  const aromaticSum = new Array(raw.length).fill(0);
  for (const b of bonds) {
    bondSum[b.a] += b.order;
    bondSum[b.b] += b.order;
    const extra = b.aromatic ? 0 : b.order - 1;
    aromaticSum[b.a] += 1 + extra;
    aromaticSum[b.b] += 1 + extra;
  }

  const atoms: Atom[] = raw.map((r, id) => ({
    id,
    el: r.el,
    charge: r.charge,
    aromatic: r.aromatic,
    hCount: r.explicitH
      ? r.hCount
      : r.aromatic
        ? implicitHydrogens(r.el, aromaticSum[id] + 1, r.charge, true)
        : implicitHydrogens(r.el, bondSum[id], r.charge),
    pos: v3(),
    placed: false,
    fromImplicitH: false,
  }));

  const allBonds: Bond[] = bonds.slice();
  const heavyCount = atoms.length;
  for (let id = 0; id < heavyCount; id++) {
    const n = atoms[id].hCount;
    for (let k = 0; k < n; k++) {
      const hid = atoms.length;
      atoms.push({
        id: hid, el: 'H', charge: 0, aromatic: false, hCount: 0,
        pos: v3(), placed: false, fromImplicitH: true,
      });
      allBonds.push({ a: id, b: hid, order: 1, aromatic: false, ringClosure: false, delocalized: false });
    }
  }

  return { atoms, bonds: allBonds };
}

/** Списки соседей для каждого атома. */
export function neighborLists(mol: Molecule): number[][] {
  const nb: number[][] = mol.atoms.map(() => []);
  for (const b of mol.bonds) {
    nb[b.a].push(b.b);
    nb[b.b].push(b.a);
  }
  return nb;
}

/** Сумма кратностей связей у атома. */
export function bondOrderSum(mol: Molecule, id: number): number {
  let s = 0;
  for (const b of mol.bonds) {
    if (b.a === id || b.b === id) s += b.order;
  }
  return s;
}

export function findBond(mol: Molecule, a: number, b: number): Bond | undefined {
  return mol.bonds.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a));
}
