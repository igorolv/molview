/**
 * Выгрузка эталонных результатов для сверки с портом на C++.
 * Запуск: npm run reference
 *
 * Сохраняются только величины, НЕ зависящие от поворота молекулы в пространстве:
 * гибридизации, стерические числа, длины связей, валентные углы. Координаты
 * специально не сохраняются — минимизация энергии может прийти к тому же самому
 * минимуму, повёрнутому как угодно, и сравнивать координаты бессмысленно.
 */
import { writeFileSync } from 'node:fs';
import { database, moleculeFromEntry } from '../src/chem/database';
import { analyze } from '../src/chem/analyze';
import { dist } from '../src/chem/vec';

const round = (x: number, digits = 2): number => Number(x.toFixed(digits));

const molecules = database().map((entry) => {
  const analysis = analyze(moleculeFromEntry(entry));
  const mol = analysis.molecule;

  return {
    id: entry.id,
    smiles: entry.smiles,
    formula: analysis.formula,
    mass: round(analysis.mass, 3),
    charge: mol.atoms.reduce((s, a) => s + a.charge, 0),
    rings: analysis.rings
      .map((r) => ({ size: r.atoms.length, aromatic: r.aromatic }))
      .sort((a, b) => a.size - b.size || Number(a.aromatic) - Number(b.aromatic)),
    polar: analysis.polar,

    atoms: analysis.atoms.map((a) => ({
      i: a.id,
      el: a.el,
      q: round(a.charge, 4),
      sigma: a.sigma,
      lp: a.lonePairs,
      sn: a.steric,
      axe: a.axe,
      hybrid: a.hybrid,
      eGeom: a.electronGeom,
      mGeom: a.molecularGeom,
    })),

    bonds: mol.bonds.map((b) => ({
      a: b.a,
      b: b.b,
      order: round(b.order, 4),
      aromatic: b.aromatic,
      delocalized: b.delocalized,
      length: round(dist(mol.atoms[b.a].pos, mol.atoms[b.b].pos), 3),
    })),

    // углы отсортированы, чтобы порядок обхода не влиял на сравнение
    angles: analysis.angles
      .map((r) => ({
        i: r.i, c: r.center, j: r.j,
        label: r.label,
        predicted: r.predicted === null ? null : round(r.predicted, 2),
        actual: round(r.actual, 2),
      }))
      .sort((x, y) => x.c - y.c || x.i - y.i || x.j - y.j),
  };
});

const out = {
  _comment: 'Эталон, выгруженный из версии на TypeScript. Порт на C++ должен '
    + 'воспроизводить эти значения. Допуск: длины связей 0,01 Å, углы 0,05°, '
    + 'всё остальное — точное совпадение.',
  generated: new Date().toISOString().slice(0, 10),
  count: molecules.length,
  molecules,
};

const file = 'reference.json';
writeFileSync(file, JSON.stringify(out, null, 1), 'utf8');

const atoms = molecules.reduce((s, m) => s + m.atoms.length, 0);
const angles = molecules.reduce((s, m) => s + m.angles.length, 0);
console.log(`Записано ${file}: ${molecules.length} молекул, ${atoms} атомов, ${angles} валентных углов.`);
