import type { Vec3 } from './vec';

/** Атом в графе молекулы. */
export interface Atom {
  id: number;
  el: string;
  /** Формальный заряд. */
  charge: number;
  /** Входит ли атом в ароматическую систему. */
  aromatic: boolean;
  /** Число неявных атомов водорода (до их развёртывания в явные). */
  hCount: number;
  /** Координаты, Å. Заполняются построителем геометрии. */
  pos: Vec3;
  placed: boolean;
  /** true — атом появился при развёртывании неявных водородов. */
  fromImplicitH: boolean;
}

/** Связь. order: 1 / 2 / 3, для ароматических — 1.5. */
export interface Bond {
  a: number;
  b: number;
  order: number;
  aromatic: boolean;
  /** Признак замыкающей цикл связи (для отладки построителя). */
  ringClosure: boolean;
  /** Кратность усреднена из-за делокализации (нитрат-ион, карбоксилат). */
  delocalized: boolean;
}

export interface Molecule {
  atoms: Atom[];
  bonds: Bond[];
  /** Заголовок из базы, если молекула найдена в справочнике. */
  meta?: MoleculeMeta;
}

/** Запись справочника molecules.json. */
export interface MoleculeMeta {
  id: string;
  name: string;
  smiles: string;
  cat: string;
  syn?: string[];
  /** Экспериментальные валентные углы, ключ вида "H-O-H". Можно указать несколько значений. */
  exp?: Record<string, number | number[]>;
  /** Экспериментальный дипольный момент, Д. */
  dipole?: number;
  note?: string;
  /** Брутто-формула, вычисляется при загрузке базы. */
  formula?: string;
}

/** Результат применения теории Гиллеспи к одному атому. */
export interface AtomAnalysis {
  id: number;
  el: string;
  charge: number;
  /** Число σ-связей = число соседей. */
  sigma: number;
  /** Число неподелённых электронных пар. */
  lonePairs: number;
  /** Стерическое число = σ + НЭП. */
  steric: number;
  /** Обозначение типа: AX3E1. */
  axe: string;
  hybrid: string;
  /** Геометрия электронных пар. */
  electronGeom: string;
  /** Геометрия молекулы (по положениям ядер). */
  molecularGeom: string;
  /** Идеальный угол для электронной геометрии. */
  idealAngle: number | null;
  /** Угол с учётом поправки на неподелённые пары. */
  predictedAngle: number | null;
  /** Является ли атом центральным (два и более соседа). */
  isCentral: boolean;
  /** Направления неподелённых пар (единичные векторы). */
  lonePairDirs: Vec3[];
  /** Направления гибридных орбиталей (все, включая связывающие). */
  orbitalDirs: Vec3[];
  /** Пояснение, если модель даёт нетипичный результат. */
  warning?: string;
}

/** Один валентный угол. */
export interface AngleRecord {
  i: number;
  center: number;
  j: number;
  label: string;
  /** Угол, измеренный по построенным координатам. */
  actual: number;
  /** Предсказание ОЭПВО. */
  predicted: number | null;
  /** Экспериментальное значение из справочника. */
  experimental: number | null;
  /** Оба атома входят в один цикл — угол задан геометрией кольца. */
  inRing: boolean;
}

export interface Ring {
  atoms: number[];
  aromatic: boolean;
  planar: boolean;
}

export interface Analysis {
  molecule: Molecule;
  formula: string;
  formulaHtml: string;
  mass: number;
  atoms: AtomAnalysis[];
  angles: AngleRecord[];
  rings: Ring[];
  /** Расчётный вектор дипольного момента (качественная оценка). */
  dipoleVec: Vec3;
  dipoleValue: number;
  polar: boolean;
  /** Сводка по гибридизациям, например { sp3: 2, sp2: 1 }. */
  hybridSummary: Record<string, number>;
}
