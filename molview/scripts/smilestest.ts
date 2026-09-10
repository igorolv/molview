/**
 * Проверка разбора и построения для молекул, которых НЕТ в базе.
 * Запуск: npm run test:smiles
 *
 * Именно эта проверка показывает, что программа — не листалка готовых картинок:
 * структура строится по строке SMILES для любой молекулы.
 */
import { parseSmiles, SmilesError } from '../src/chem/smiles';
import { analyze } from '../src/chem/analyze';
import { dist } from '../src/chem/vec';
import { bondLength } from '../src/chem/periodic';

interface Case { smiles: string; name: string; expect?: { rings?: number; formula?: string } }

const CASES: Case[] = [
  { smiles: 'CC(=O)Oc1ccccc1C(=O)O', name: 'Аспирин', expect: { rings: 1, formula: 'C9H8O4' } },
  { smiles: 'CC(=O)Nc1ccc(O)cc1', name: 'Парацетамол', expect: { rings: 1, formula: 'C8H9NO2' } },
  { smiles: 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', name: 'Кофеин', expect: { rings: 2, formula: 'C8H10N4O2' } },
  { smiles: 'c1ccc2cc3ccccc3cc2c1', name: 'Антрацен', expect: { rings: 3, formula: 'C14H10' } },
  { smiles: 'c1ccc2[nH]ccc2c1', name: 'Индол', expect: { rings: 2 } },
  { smiles: 'c1cncnc1', name: 'Пиримидин', expect: { rings: 1, formula: 'C4H4N2' } },
  { smiles: 'C1C2CC3CC1CC(C2)C3', name: 'Адамантан (мостиковый)', expect: { formula: 'C10H16' } },
  { smiles: 'C1CC2CCC1C2', name: 'Норборнан (мостиковый)' },
  { smiles: 'OCC(O)CO', name: 'Глицерин', expect: { rings: 0, formula: 'C3H8O3' } },
  { smiles: 'CC(N)C(=O)O', name: 'Аланин', expect: { formula: 'C3H7NO2' } },
  { smiles: 'CCOC(C)=O', name: 'Этилацетат', expect: { formula: 'C4H8O2' } },
  { smiles: 'OC(=O)CCC(=O)O', name: 'Янтарная кислота' },
  { smiles: 'c1ccccc1[N+](=O)[O-]', name: 'Нитробензол', expect: { rings: 1 } },
  { smiles: 'Cc1ccc(cc1)S(=O)(=O)O', name: 'п-Толуолсульфокислота', expect: { rings: 1 } },
  { smiles: 'F[Br](F)(F)(F)F', name: 'Пентафторид брома (sp³d²)' },
  { smiles: 'F[I](F)(F)(F)(F)(F)F', name: 'Гептафторид иода (sp³d³)' },
  { smiles: '[O-][Cl](=O)(=O)=O', name: 'Перхлорат-ион' },
  { smiles: 'ClS(Cl)=O', name: 'Хлористый тионил' },
  { smiles: 'ClP(Cl)(Cl)=O', name: 'Хлорокись фосфора' },
  { smiles: '[Na+].[Cl-]', name: 'Разорванная запись (две частицы)' },
  { smiles: 'C1=CC=CC=C1', name: 'Бензол в форме Кекуле', expect: { rings: 1, formula: 'C6H6' } },
  { smiles: 'CC(C)(C)C(C)(C)C', name: 'Гексаметилэтан (сильно загромождён)' },
  { smiles: 'OCC1OC(O)C(O)C(O)C1O', name: 'Глюкоза' },
  { smiles: 'C/C=C/C', name: 'Бутен-2 со стереометками' },
  { smiles: 'CC#CC', name: 'Бутин-2' },
];

const BAD: { smiles: string; why: string }[] = [
  { smiles: 'C1CC', why: 'незамкнутый цикл' },
  { smiles: 'CC(', why: 'незакрытая скобка' },
  { smiles: 'CC)', why: 'лишняя скобка' },
  { smiles: 'CQC', why: 'неизвестный символ' },
  { smiles: '[Zz]', why: 'неизвестный элемент' },
  { smiles: '', why: 'пустая строка' },
];

const pad = (s: string, n: number): string => (s.length >= n ? s.slice(0, n) : s + ' '.repeat(n - s.length));
let failures = 0;
const problems: string[] = [];
const fail = (m: string): void => { failures++; problems.push(m); };

console.log('\n=== ПОСТРОЕНИЕ ПО SMILES (молекул нет в базе) ===\n');
console.log(pad('Молекула', 32) + pad('Формула', 12) + pad('Атомов', 8) + pad('Циклов', 8) + pad('Мин. расст.', 13) + 'Макс. ошибка связи');
console.log('-'.repeat(96));

for (const test of CASES) {
  let analysis;
  try {
    analysis = analyze(parseSmiles(test.smiles));
  } catch (err) {
    fail(`${test.name}: ${(err as Error).message}`);
    console.log(pad(test.name, 32) + 'ОШИБКА: ' + (err as Error).message);
    continue;
  }

  const mol = analysis.molecule;
  let minDist = Infinity;
  for (let i = 0; i < mol.atoms.length; i++) {
    for (let j = i + 1; j < mol.atoms.length; j++) {
      minDist = Math.min(minDist, dist(mol.atoms[i].pos, mol.atoms[j].pos));
    }
  }
  let maxBondError = 0;
  for (const b of mol.bonds) {
    const expected = bondLength(mol.atoms[b.a].el, mol.atoms[b.b].el, b.order);
    maxBondError = Math.max(maxBondError, Math.abs(dist(mol.atoms[b.a].pos, mol.atoms[b.b].pos) - expected));
  }

  console.log(
    pad(test.name, 32) + pad(analysis.formula, 12) +
    pad(String(mol.atoms.length), 8) + pad(String(analysis.rings.length), 8) +
    pad(minDist.toFixed(2) + ' Å', 13) + maxBondError.toFixed(2) + ' Å',
  );

  if (minDist < 0.85) fail(`${test.name}: атомы слиплись (${minDist.toFixed(2)} Å)`);
  if (maxBondError > 0.25) fail(`${test.name}: длина связи отклонилась на ${maxBondError.toFixed(2)} Å`);
  if (mol.atoms.some((a) => !Number.isFinite(a.pos.x + a.pos.y + a.pos.z))) {
    fail(`${test.name}: нечисловые координаты`);
  }
  if (test.expect?.rings !== undefined && analysis.rings.length !== test.expect.rings) {
    fail(`${test.name}: циклов ${analysis.rings.length}, ожидалось ${test.expect.rings}`);
  }
  if (test.expect?.formula !== undefined && analysis.formula !== test.expect.formula) {
    fail(`${test.name}: формула ${analysis.formula}, ожидалась ${test.expect.formula}`);
  }
}

console.log('\n=== ОБРАБОТКА ОШИБОК ВВОДА ===\n');
for (const bad of BAD) {
  try {
    parseSmiles(bad.smiles);
    fail(`«${bad.smiles}» (${bad.why}) разобралось без ошибки`);
    console.log(pad(`«${bad.smiles}»`, 16) + 'ОШИБКА: разбор прошёл, хотя не должен был');
  } catch (err) {
    const ok = err instanceof SmilesError;
    if (!ok) fail(`«${bad.smiles}»: исключение не того типа`);
    console.log(pad(`«${bad.smiles}»`, 16) + pad(bad.why, 22) + '→ ' + (err as Error).message);
  }
}

console.log('\n' + '='.repeat(96));
if (failures === 0) {
  console.log(`ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ. Разобрано молекул вне базы: ${CASES.length}.`);
} else {
  console.log(`ОБНАРУЖЕНО ПРОБЛЕМ: ${failures}\n`);
  problems.forEach((p) => console.log('  • ' + p));
  process.exitCode = 1;
}
