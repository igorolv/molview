import type { Vec3 } from './vec';
import { v3, norm, sub, mul, add, dot, len } from './vec';
import { element } from './periodic';

/**
 * ТЕОРИЯ ОТТАЛКИВАНИЯ ЭЛЕКТРОННЫХ ПАР ВАЛЕНТНОЙ ОБОЛОЧКИ (ОЭПВО / VSEPR)
 *
 * Главный принцип реализации: расположение электронных групп вокруг атома
 * не берётся из готовой таблицы, а НАХОДИТСЯ как минимум энергии их взаимного
 * отталкивания (модель точек на сфере, Р. Бартелл). Энергия:
 *
 *     E = Σ (w_i · w_j) / r_ij^6
 *
 * где w — «размер» электронной группы: неподелённая пара крупнее связывающей,
 * кратная связь крупнее одинарной. Из одной этой минимизации получаются все
 * правила учебника:
 *   • 2 группы → линия, 3 → треугольник, 4 → тетраэдр, 6 → октаэдр;
 *   • неподелённые пары сжимают углы (109,5° → 107° → 104,5°);
 *   • кратная связь «раздвигает» соседей (в H₂C=O угол H-C-H меньше 120°);
 *   • в тригональной бипирамиде пара уходит в ЭКВАТОРИАЛЬНУЮ позицию;
 *   • две пары в октаэдре становятся ТРАНС, давая плоский квадрат.
 *
 * Для стерического числа 5 и больше свободная минимизация даёт небольшие
 * искажения там, где симметрия требует точных 90° и 120°, поэтому там
 * используется точный шаблон геометрии, а расположение неподелённых пар
 * по его вершинам всё равно выбирается ПО МИНИМУМУ ЭНЕРГИИ — то есть
 * правила «пара в экватор» и «пары в транс» программа по-прежнему выводит,
 * а не берёт готовыми.
 */

// --- эмпирические «размеры» электронных групп -------------------------------
// Откалиброваны по справочным углам H₂O (104,5°), NH₃ (107,0°), H₂C=O (116,5°).
export const WEIGHT_LONE_PAIR = 1.45;
export const WEIGHT_SINGLE_BOND = 1.0;
/** Прибавка за каждую единицу кратности сверх одинарной связи. */
export const WEIGHT_PER_EXTRA_ORDER = 0.18;
/** Показатель степени в законе отталкивания E ~ 1/r^n. */
const REPULSION_EXPONENT = 6;

export function bondWeight(order: number): number {
  return WEIGHT_SINGLE_BOND + WEIGHT_PER_EXTRA_ORDER * Math.max(0, order - 1);
}

// ---------------------------------------------------------------------------
// Электронное строение: неподелённые пары
// ---------------------------------------------------------------------------

/**
 * Неподелённые пары по схеме Льюиса:
 *   НЭП = (валентные электроны − формальный заряд − Σ кратностей связей) / 2
 *
 * У ароматических атомов сумма кратностей дробная (1,5 + 1,5 = 3), и деление
 * даёт «половинку» — она означает пару, ушедшую в π-систему кольца. Такую
 * половинку отбрасываем: у азота в пирроле остаётся 0 пар в плоскости кольца,
 * у кислорода во фуране — 1, у азота в пиридине — 1. Это в точности
 * соответствует химии этих соединений.
 */
export function lonePairCount(el: string, charge: number, bondOrderSum: number): number {
  const ve = element(el).ve;
  const raw = (ve - charge - bondOrderSum) / 2;
  return Math.max(0, Math.round(raw - 1e-3));
}

export interface GeometryInfo {
  electronGeom: string;
  molecularGeom: string;
  /** Идеальный угол электронной геометрии, «из учебника». */
  idealAngle: number | null;
  /** Из чего складывается окружение атома. Про сам атом, без примеров. */
  hint: string;
  /**
   * Вещества с таким же окружением центрального атома. Это именно примеры
   * «из учебника», к разбираемой молекуле они отношения не имеют, поэтому
   * в интерфейсе выводятся отдельной строкой и с явной пометкой.
   */
  examples: string;
}

const FALLBACK: GeometryInfo = {
  electronGeom: 'не определена',
  molecularGeom: 'не определена',
  idealAngle: null,
  hint: 'Стерическое число выходит за рамки модели ОЭПВО.',
  examples: '',
};

/** Ключ: «число σ-связей : число неподелённых пар». */
const GEOMETRY: Record<string, GeometryInfo> = {
  '1:0': { electronGeom: 'одна связь', molecularGeom: 'концевой атом', idealAngle: null, hint: 'Единственный сосед — валентного угла нет.', examples: 'водород в H₂O и CH₄' },
  '1:1': { electronGeom: 'линейная', molecularGeom: 'концевой атом', idealAngle: null, hint: 'Одна связь и одна неподелённая пара.', examples: 'азот в HC≡N, углерод в CO' },
  '1:2': { electronGeom: 'тригональная', molecularGeom: 'концевой атом', idealAngle: null, hint: 'Одна связь и две неподелённые пары.', examples: 'концевой кислород двойной связи — в C=O, SO₂, NO₃⁻' },
  '1:3': { electronGeom: 'тетраэдрическая', molecularGeom: 'концевой атом', idealAngle: null, hint: 'Одна связь и три неподелённые пары.', examples: 'галогены в HCl и CCl₄, кислород в OH⁻ и NO₃⁻' },

  '2:0': { electronGeom: 'линейная', molecularGeom: 'линейная', idealAngle: 180, hint: 'Две связи без пар — они расходятся на 180°.', examples: 'CO₂, BeCl₂, ацетилен' },
  '2:1': { electronGeom: 'тригональная', molecularGeom: 'угловая', idealAngle: 120, hint: 'Две связи и одна пара — частица угловая.', examples: 'SO₂, озон, нитрит-ион' },
  '2:2': { electronGeom: 'тетраэдрическая', molecularGeom: 'угловая', idealAngle: 109.5, hint: 'Две связи и две пары — частица угловая.', examples: 'вода, сероводород, спирты и эфиры' },
  '2:3': { electronGeom: 'тригонально-бипирамидальная', molecularGeom: 'линейная', idealAngle: 180, hint: 'Две связи и три пары: пары в экваторе, связи на оси.', examples: 'XeF₂' },

  '3:0': { electronGeom: 'тригональная', molecularGeom: 'плоская треугольная', idealAngle: 120, hint: 'Три связи без пар — плоский треугольник.', examples: 'BF₃, SO₃, нитрат-ион, любой sp²-углерод' },
  '3:1': { electronGeom: 'тетраэдрическая', molecularGeom: 'тригонально-пирамидальная', idealAngle: 109.5, hint: 'Три связи и одна пара — пирамида.', examples: 'аммиак, PCl₃, ион гидроксония, сера в ДМСО' },
  '3:2': { electronGeom: 'тригонально-бипирамидальная', molecularGeom: 'Т-образная', idealAngle: 90, hint: 'Три связи и две пары: обе пары в экваторе, форма «Т».', examples: 'ClF₃' },

  '4:0': { electronGeom: 'тетраэдрическая', molecularGeom: 'тетраэдрическая', idealAngle: 109.5, hint: 'Четыре связи без пар — правильный тетраэдр.', examples: 'метан, ион аммония, сульфат-ион' },
  '4:1': { electronGeom: 'тригонально-бипирамидальная', molecularGeom: 'качели (дисфеноид)', idealAngle: 120, hint: 'Четыре связи и одна пара: пара в экваторе, форма качелей.', examples: 'SF₄' },
  '4:2': { electronGeom: 'октаэдрическая', molecularGeom: 'квадратная', idealAngle: 90, hint: 'Четыре связи и две пары: пары транс, связи в квадрате.', examples: 'XeF₄' },

  '5:0': { electronGeom: 'тригонально-бипирамидальная', molecularGeom: 'тригонально-бипирамидальная', idealAngle: 120, hint: 'Пять связей: 3 экваториальные и 2 аксиальные.', examples: 'PCl₅' },
  '5:1': { electronGeom: 'октаэдрическая', molecularGeom: 'квадратно-пирамидальная', idealAngle: 90, hint: 'Пять связей и одна пара: квадратная пирамида.', examples: 'IF₅, BrF₅' },

  '6:0': { electronGeom: 'октаэдрическая', molecularGeom: 'октаэдрическая', idealAngle: 90, hint: 'Шесть равноценных связей — правильный октаэдр.', examples: 'SF₆' },

  '7:0': { electronGeom: 'пентагонально-бипирамидальная', molecularGeom: 'пентагонально-бипирамидальная', idealAngle: 72, hint: 'Семь связей — редкая геометрия.', examples: 'IF₇' },
};

export function geometryFor(sigma: number, lonePairs: number): GeometryInfo {
  return GEOMETRY[`${sigma}:${lonePairs}`] ?? FALLBACK;
}

export function hybridFor(steric: number): string {
  switch (steric) {
    case 1: return '—';
    case 2: return 'sp';
    case 3: return 'sp²';
    case 4: return 'sp³';
    case 5: return 'sp³d';
    case 6: return 'sp³d²';
    case 7: return 'sp³d³';
    default: return '—';
  }
}

/** Пояснение, из чего складывается гибридизация. */
export function hybridExplanation(steric: number): string {
  switch (steric) {
    case 2: return 'одна s- и одна p-орбиталь дают две sp-орбитали под 180°';
    case 3: return 'одна s- и две p-орбитали дают три sp²-орбитали под 120° в одной плоскости';
    case 4: return 'одна s- и три p-орбитали дают четыре sp³-орбитали, направленные к вершинам тетраэдра';
    case 5: return 'к s- и p-орбиталям добавляется одна d-орбиталь: пять sp³d-орбиталей образуют тригональную бипирамиду';
    case 6: return 'к s- и p-орбиталям добавляются две d-орбитали: шесть sp³d²-орбиталей образуют октаэдр';
    case 7: return 'семь sp³d³-орбиталей образуют пентагональную бипирамиду';
    default: return '';
  }
}

const SUB = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'];
const subscript = (n: number): string => String(n).split('').map((d) => SUB[+d]).join('');

/** Обозначение типа молекулы по Гиллеспи: AX₄E₂. */
export function axeNotation(sigma: number, lonePairs: number): string {
  let s = 'A';
  if (sigma > 0) s += 'X' + (sigma > 1 ? subscript(sigma) : '');
  if (lonePairs > 0) s += 'E' + (lonePairs > 1 ? subscript(lonePairs) : '');
  return s;
}

// ---------------------------------------------------------------------------
// Расстановка электронных групп вокруг атома
// ---------------------------------------------------------------------------

function repulsionEnergy(points: Vec3[], weights: number[]): number {
  let e = 0;
  for (let i = 0; i < points.length; i++) {
    for (let j = i + 1; j < points.length; j++) {
      const d = Math.max(1e-6, len(sub(points[i], points[j])));
      e += (weights[i] * weights[j]) / Math.pow(d, REPULSION_EXPONENT);
    }
  }
  return e;
}

/** Детерминированный генератор — одна молекула всегда строится одинаково. */
function makeRandom(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

function relaxOnSphere(points: Vec3[], weights: number[], steps: number, step0: number, decay: number): void {
  const n = points.length;
  let step = step0;
  for (let iter = 0; iter < steps; iter++) {
    const forces: Vec3[] = points.map(() => v3());
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const diff = sub(points[i], points[j]);
        const d = Math.max(1e-6, len(diff));
        const magnitude =
          (REPULSION_EXPONENT * weights[i] * weights[j]) / Math.pow(d, REPULSION_EXPONENT + 2);
        const f = mul(diff, magnitude);
        forces[i] = add(forces[i], f);
        forces[j] = sub(forces[j], f);
      }
    }
    for (let i = 0; i < n; i++) {
      // сохраняем только касательную составляющую: точка обязана остаться на сфере
      const p = points[i];
      const f = forces[i];
      const tangential = sub(f, mul(p, dot(f, p)));
      points[i] = norm(add(p, mul(tangential, step)));
    }
    step *= decay;
  }
}

/** Точные шаблоны геометрии для стерических чисел 5–7. */
function idealTemplate(n: number): Vec3[] {
  const s32 = Math.sqrt(3) / 2;
  switch (n) {
    case 5: // тригональная бипирамида: 3 экваториальные + 2 аксиальные
      return [
        v3(1, 0, 0), v3(-0.5, s32, 0), v3(-0.5, -s32, 0),
        v3(0, 0, 1), v3(0, 0, -1),
      ];
    case 6: // октаэдр
      return [
        v3(0, 0, 1), v3(0, 0, -1),
        v3(1, 0, 0), v3(-1, 0, 0),
        v3(0, 1, 0), v3(0, -1, 0),
      ];
    case 7: { // пентагональная бипирамида
      const eq: Vec3[] = [];
      for (let k = 0; k < 5; k++) {
        const a = (2 * Math.PI * k) / 5;
        eq.push(v3(Math.cos(a), Math.sin(a), 0));
      }
      return [...eq, v3(0, 0, 1), v3(0, 0, -1)];
    }
    default:
      return [];
  }
}

/** Все перестановки индексов 0..n-1 (n ≤ 7, поэтому это дёшево). */
function permutations(n: number): number[][] {
  const result: number[][] = [];
  const current: number[] = [];
  const used = new Array(n).fill(false);
  const walk = (): void => {
    if (current.length === n) { result.push(current.slice()); return; }
    for (let i = 0; i < n; i++) {
      if (used[i]) continue;
      used[i] = true;
      current.push(i);
      walk();
      current.pop();
      used[i] = false;
    }
  };
  walk();
  return result;
}

const permCache = new Map<number, number[][]>();
function cachedPermutations(n: number): number[][] {
  let p = permCache.get(n);
  if (!p) { p = permutations(n); permCache.set(n, p); }
  return p;
}

const directionCache = new Map<string, Vec3[]>();

/**
 * Направления электронных групп вокруг атома.
 * weights[k] — «размер» k-й группы. Порядок групп в ответе совпадает со входом.
 */
export function arrangeGroups(weights: number[]): Vec3[] {
  const n = weights.length;
  if (n === 0) return [];
  if (n === 1) return [v3(0, 0, 1)];

  const key = weights.map((w) => w.toFixed(4)).join('|');
  const cached = directionCache.get(key);
  if (cached) return cached.map((p) => ({ ...p }));

  let result: Vec3[];
  if (n <= 4) {
    result = minimizeFreely(weights);
  } else if (n <= 7) {
    result = fitToTemplate(weights, idealTemplate(n));
  } else {
    result = minimizeFreely(weights);
  }

  canonicalize(result);
  directionCache.set(key, result.map((p) => ({ ...p })));
  return result.map((p) => ({ ...p }));
}

/** Свободная минимизация энергии на сфере (стерическое число ≤ 4). */
function minimizeFreely(weights: number[]): Vec3[] {
  const n = weights.length;
  let best: Vec3[] | null = null;
  let bestEnergy = Infinity;

  for (let attempt = 0; attempt < 24; attempt++) {
    const rnd = makeRandom(0x9e3779b9 + attempt * 7919 + n * 104729);
    const points: Vec3[] = [];
    for (let i = 0; i < n; i++) {
      const u = rnd() * 2 - 1;
      const phi = rnd() * Math.PI * 2;
      const r = Math.sqrt(Math.max(0, 1 - u * u));
      points.push(v3(r * Math.cos(phi), r * Math.sin(phi), u));
    }
    relaxOnSphere(points, weights, 3000, 0.3, 0.9985);
    relaxOnSphere(points, weights, 600, 0.02, 0.999);
    const e = repulsionEnergy(points, weights);
    if (e < bestEnergy - 1e-10) { bestEnergy = e; best = points; }
  }
  return best!;
}

/**
 * Стерическое число 5–7: берём точный шаблон геометрии, но РАСПРЕДЕЛЕНИЕ групп
 * по его вершинам выбираем по минимуму энергии отталкивания. Именно так
 * программа выводит, что неподелённая пара идёт в экваториальную позицию
 * тригональной бипирамиды, а две пары в октаэдре становятся друг напротив друга.
 */
function fitToTemplate(weights: number[], template: Vec3[]): Vec3[] {
  const n = weights.length;
  let bestPerm: number[] | null = null;
  let bestEnergy = Infinity;

  for (const perm of cachedPermutations(n)) {
    const points = perm.map((slot) => template[slot]);
    const e = repulsionEnergy(points, weights);
    if (e < bestEnergy - 1e-12) { bestEnergy = e; bestPerm = perm; }
  }
  return bestPerm!.map((slot) => ({ ...template[slot] }));
}

/**
 * Приведение к воспроизводимой ориентации: первая группа смотрит вдоль +Z,
 * вторая ложится в плоскость XZ. Иначе одна и та же молекула при каждом
 * запуске выглядела бы повёрнутой по-разному.
 */
function canonicalize(points: Vec3[]): void {
  if (points.length < 2) return;

  const first = norm(points[0]);
  const target = v3(0, 0, 1);
  const axis = v3(
    first.y * target.z - first.z * target.y,
    first.z * target.x - first.x * target.z,
    first.x * target.y - first.y * target.x,
  );
  const axisLen = len(axis);
  const cosA = dot(first, target);

  if (axisLen > 1e-9) {
    const a = mul(axis, 1 / axisLen);
    const angle = Math.acos(Math.max(-1, Math.min(1, cosA)));
    const c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
    for (let i = 0; i < points.length; i++) {
      const p = points[i];
      points[i] = v3(
        (t * a.x * a.x + c) * p.x + (t * a.x * a.y - s * a.z) * p.y + (t * a.x * a.z + s * a.y) * p.z,
        (t * a.x * a.y + s * a.z) * p.x + (t * a.y * a.y + c) * p.y + (t * a.y * a.z - s * a.x) * p.z,
        (t * a.x * a.z - s * a.y) * p.x + (t * a.y * a.z + s * a.x) * p.y + (t * a.z * a.z + c) * p.z,
      );
    }
  } else if (cosA < 0) {
    for (let i = 0; i < points.length; i++) points[i] = mul(points[i], -1);
  }
  points[0] = v3(0, 0, 1);

  const second = points[1];
  const planar = Math.hypot(second.x, second.y);
  if (planar > 1e-9) {
    const c = second.x / planar;
    const s = second.y / planar;
    for (let i = 0; i < points.length; i++) {
      const p = points[i];
      points[i] = v3(c * p.x + s * p.y, -s * p.x + c * p.y, p.z);
    }
  }
}
