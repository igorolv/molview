import moleculesData from '../data/molecules.json';
import type { Molecule, MoleculeMeta } from './types';
import { parseSmiles, SmilesError } from './smiles';
import { moleculeCompositionKey, parseFormulaQuery, plainFormula, prettyFormula } from './formula';

export interface DbEntry extends MoleculeMeta {
  compositionKey: string;
  formula: string;
  formulaPretty: string;
  searchText: string;
}

function buildIndex(): DbEntry[] {
  const raw = moleculesData.molecules as MoleculeMeta[];
  const entries: DbEntry[] = [];

  for (const meta of raw) {
    let mol: Molecule;
    try {
      mol = parseSmiles(meta.smiles);
    } catch (err) {
      // молекула с ошибкой в базе не должна ронять всю программу
      console.warn(`Молекула «${meta.name}» пропущена: ${(err as Error).message}`);
      continue;
    }
    const formula = plainFormula(mol);
    const searchText = [meta.name, ...(meta.syn ?? []), formula, meta.smiles]
      .join(' ')
      .toLowerCase()
      .replace(/ё/g, 'е');

    entries.push({
      ...meta,
      compositionKey: moleculeCompositionKey(mol),
      formula,
      formulaPretty: prettyFormula(mol),
      searchText,
    });
  }
  return entries;
}

let cachedIndex: DbEntry[] | null = null;

export function database(): DbEntry[] {
  if (cachedIndex === null) cachedIndex = buildIndex();
  return cachedIndex;
}

export function entryById(id: string): DbEntry | undefined {
  return database().find((e) => e.id === id);
}

/** Категории в порядке появления в базе. */
export function categories(): { name: string; entries: DbEntry[] }[] {
  const map = new Map<string, DbEntry[]>();
  for (const e of database()) {
    if (!map.has(e.cat)) map.set(e.cat, []);
    map.get(e.cat)!.push(e);
  }
  return [...map.entries()].map(([name, entries]) => ({ name, entries }));
}

export type Resolution =
  | { kind: 'entries'; entries: DbEntry[]; reason: 'name' | 'formula' }
  | { kind: 'smiles'; molecule: Molecule; smiles: string }
  | { kind: 'error'; message: string; hint?: string };

const normalize = (s: string): string => s.trim().toLowerCase().replace(/ё/g, 'е');

/**
 * Разбор пользовательского запроса.
 *
 * Порядок попыток: название → брутто-формула → SMILES.
 * Формула специально идёт раньше SMILES: на запрос «C2H6O» программа должна
 * показать ОБА изомера, а не молчаливо разобрать строку как SMILES.
 */
export function resolveQuery(input: string): Resolution {
  const query = input.trim();
  if (!query) return { kind: 'error', message: 'Введите формулу, название или SMILES' };

  const db = database();
  const needle = normalize(query);

  // 1. точное совпадение названия или синонима
  const exact = db.filter(
    (e) => normalize(e.name) === needle || (e.syn ?? []).some((s) => normalize(s) === needle),
  );
  if (exact.length > 0) return { kind: 'entries', entries: exact, reason: 'name' };

  // 2. брутто-формула — сюда попадают все изомеры
  const key = parseFormulaQuery(query);
  if (key) {
    const byFormula = db.filter((e) => e.compositionKey === key);
    if (byFormula.length > 0) return { kind: 'entries', entries: byFormula, reason: 'formula' };
  }

  // 3. вхождение в название
  const partial = db.filter((e) => e.searchText.includes(needle));
  if (partial.length > 0) return { kind: 'entries', entries: partial, reason: 'name' };

  // 4. SMILES
  try {
    const molecule = parseSmiles(query);
    return { kind: 'smiles', molecule, smiles: query };
  } catch (err) {
    if (!(err instanceof SmilesError)) {
      return { kind: 'error', message: (err as Error).message };
    }
    if (key) {
      return {
        kind: 'error',
        message: `Формулы ${query} нет в справочнике`,
        hint: 'Строение можно задать напрямую в виде SMILES — например, CCO для этанола '
          + 'или c1ccccc1 для бензола. Подсказка по синтаксису — кнопка «?» справа вверху.',
      };
    }
    // запрос без скобок и знаков связи скорее всего задумывался как формула
    const looksLikeFormula = /^[A-Za-z][A-Za-z0-9]*$/.test(query);
    return {
      kind: 'error',
      message: looksLikeFormula
        ? `Не похоже ни на формулу, ни на название: ${err.message.toLowerCase()}`
        : `Не удалось разобрать запрос: ${err.message}`,
      hint: looksLikeFormula
        ? 'Если это брутто-формула, проверьте символы элементов: они пишутся с заглавной буквы '
          + '(Cl, Br, Xe). Можно также ввести название по-русски или строение в виде SMILES.'
        : 'Введите брутто-формулу (H2O, C2H6O), название («вода», «бензол») '
          + 'или строение в виде SMILES (CCO, c1ccccc1).',
    };
  }
}

/** Подсказки автодополнения. */
export function suggestions(input: string, limit = 8): DbEntry[] {
  const needle = normalize(input);
  if (!needle) return [];
  const db = database();

  const starts = db.filter((e) => normalize(e.name).startsWith(needle));
  const contains = db.filter((e) => !starts.includes(e) && e.searchText.includes(needle));
  return [...starts, ...contains].slice(0, limit);
}

/** Молекулы с тем же составом — изомеры. */
export function isomersOf(entry: DbEntry): DbEntry[] {
  return database().filter((e) => e.compositionKey === entry.compositionKey && e.id !== entry.id);
}

export function moleculeFromEntry(entry: DbEntry): Molecule {
  const mol = parseSmiles(entry.smiles);
  mol.meta = entry;
  return mol;
}
