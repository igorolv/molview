/**
 * Векторная алгебра в трёхмерном пространстве.
 * Модуль полностью независим от браузера — при переносе программы
 * на C++ или Python переписывается один в один.
 */

export interface Vec3 {
  x: number;
  y: number;
  z: number;
}

export const v3 = (x = 0, y = 0, z = 0): Vec3 => ({ x, y, z });

export const add = (a: Vec3, b: Vec3): Vec3 => ({ x: a.x + b.x, y: a.y + b.y, z: a.z + b.z });
export const sub = (a: Vec3, b: Vec3): Vec3 => ({ x: a.x - b.x, y: a.y - b.y, z: a.z - b.z });
export const mul = (a: Vec3, k: number): Vec3 => ({ x: a.x * k, y: a.y * k, z: a.z * k });
export const dot = (a: Vec3, b: Vec3): number => a.x * b.x + a.y * b.y + a.z * b.z;

export const cross = (a: Vec3, b: Vec3): Vec3 => ({
  x: a.y * b.z - a.z * b.y,
  y: a.z * b.x - a.x * b.z,
  z: a.x * b.y - a.y * b.x,
});

export const len = (a: Vec3): number => Math.sqrt(dot(a, a));
export const dist = (a: Vec3, b: Vec3): number => len(sub(a, b));
export const clone = (a: Vec3): Vec3 => ({ x: a.x, y: a.y, z: a.z });

export function norm(a: Vec3): Vec3 {
  const l = len(a);
  return l < 1e-12 ? v3(0, 0, 1) : mul(a, 1 / l);
}

/** Любой единичный вектор, перпендикулярный данному. */
export function anyPerp(a: Vec3): Vec3 {
  const n = norm(a);
  const helper = Math.abs(n.x) < 0.9 ? v3(1, 0, 0) : v3(0, 1, 0);
  return norm(cross(n, helper));
}

/** Угол между векторами в градусах. */
export function angleDeg(a: Vec3, b: Vec3): number {
  const c = dot(norm(a), norm(b));
  return (Math.acos(Math.max(-1, Math.min(1, c))) * 180) / Math.PI;
}

/** Валентный угол i–c–j (в градусах) по трём точкам. */
export function bondAngle(pi: Vec3, pc: Vec3, pj: Vec3): number {
  return angleDeg(sub(pi, pc), sub(pj, pc));
}

// ---------------------------------------------------------------------------
// Матрицы 3×3 (по строкам: m[строка][столбец])
// ---------------------------------------------------------------------------

export type Mat3 = [number, number, number, number, number, number, number, number, number];

export const IDENTITY: Mat3 = [1, 0, 0, 0, 1, 0, 0, 0, 1];

export function applyMat(m: Mat3, v: Vec3): Vec3 {
  return {
    x: m[0] * v.x + m[1] * v.y + m[2] * v.z,
    y: m[3] * v.x + m[4] * v.y + m[5] * v.z,
    z: m[6] * v.x + m[7] * v.y + m[8] * v.z,
  };
}

export function matMul(a: Mat3, b: Mat3): Mat3 {
  const r = new Array(9).fill(0) as number[];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let s = 0;
      for (let k = 0; k < 3; k++) s += a[i * 3 + k] * b[k * 3 + j];
      r[i * 3 + j] = s;
    }
  }
  return r as Mat3;
}

export function transpose(m: Mat3): Mat3 {
  return [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];
}

/** Матрица поворота вокруг произвольной оси на угол (радианы). Формула Родрига. */
export function rotationAxis(axis: Vec3, angle: number): Mat3 {
  const a = norm(axis);
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  const t = 1 - c;
  const { x, y, z } = a;
  return [
    t * x * x + c, t * x * y - s * z, t * x * z + s * y,
    t * x * y + s * z, t * y * y + c, t * y * z - s * x,
    t * x * z - s * y, t * y * z + s * x, t * z * z + c,
  ];
}

/** Минимальный поворот, переводящий единичный вектор `from` в `to`. */
export function rotationFromTo(from: Vec3, to: Vec3): Mat3 {
  const a = norm(from);
  const b = norm(to);
  const d = dot(a, b);
  if (d > 0.999999) return IDENTITY;
  if (d < -0.999999) return rotationAxis(anyPerp(a), Math.PI);
  const axis = cross(a, b);
  return rotationAxis(axis, Math.acos(Math.max(-1, Math.min(1, d))));
}

/**
 * Ортонормированная система координат, построенная по двум векторам
 * (первый задаёт ось X, второй — плоскость XY). Возвращается матрица,
 * столбцы которой — базисные векторы.
 */
export function frameFrom(a: Vec3, b: Vec3): Mat3 {
  const e1 = norm(a);
  let e2 = sub(b, mul(e1, dot(b, e1)));
  if (len(e2) < 1e-8) e2 = anyPerp(e1);
  e2 = norm(e2);
  const e3 = cross(e1, e2);
  // столбцы: e1, e2, e3
  return [e1.x, e2.x, e3.x, e1.y, e2.y, e3.y, e1.z, e2.z, e3.z];
}

/**
 * Поворот, наилучшим образом совмещающий пару направлений (s1,s2) с (t1,t2).
 * Используется при достройке окружения атома, у которого уже размещены соседи.
 */
export function alignPair(s1: Vec3, s2: Vec3, t1: Vec3, t2: Vec3): Mat3 {
  const fs = frameFrom(s1, s2);
  const ft = frameFrom(t1, t2);
  return matMul(ft, transpose(fs));
}

/** Центр масс (геометрический центр) набора точек. */
export function centroid(points: Vec3[]): Vec3 {
  if (points.length === 0) return v3();
  let s = v3();
  for (const p of points) s = add(s, p);
  return mul(s, 1 / points.length);
}

/**
 * Собственные векторы симметричной матрицы 3×3 методом Якоби.
 * Нужны для разворота молекулы по главным осям инерции — чтобы
 * при первом показе она была повёрнута к зрителю самой «широкой» стороной.
 */
export function jacobiEigen(inputMatrix: Mat3): { values: number[]; vectors: Vec3[] } {
  const a = inputMatrix.slice() as Mat3;
  let v: Mat3 = IDENTITY.slice() as Mat3;

  for (let sweep = 0; sweep < 64; sweep++) {
    let off = 0;
    for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++) off += a[i * 3 + j] ** 2;
    if (off < 1e-18) break;

    for (let p = 0; p < 3; p++) {
      for (let q = p + 1; q < 3; q++) {
        const apq = a[p * 3 + q];
        if (Math.abs(apq) < 1e-18) continue;
        const theta = (a[q * 3 + q] - a[p * 3 + p]) / (2 * apq);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const s = t * c;

        for (let k = 0; k < 3; k++) {
          const akp = a[k * 3 + p];
          const akq = a[k * 3 + q];
          a[k * 3 + p] = c * akp - s * akq;
          a[k * 3 + q] = s * akp + c * akq;
        }
        for (let k = 0; k < 3; k++) {
          const apk = a[p * 3 + k];
          const aqk = a[q * 3 + k];
          a[p * 3 + k] = c * apk - s * aqk;
          a[q * 3 + k] = s * apk + c * aqk;
        }
        for (let k = 0; k < 3; k++) {
          const vkp = v[k * 3 + p];
          const vkq = v[k * 3 + q];
          v[k * 3 + p] = c * vkp - s * vkq;
          v[k * 3 + q] = s * vkp + c * vkq;
        }
      }
    }
  }

  const values = [a[0], a[4], a[8]];
  const vectors: Vec3[] = [
    v3(v[0], v[3], v[6]),
    v3(v[1], v[4], v[7]),
    v3(v[2], v[5], v[8]),
  ];
  return { values, vectors };
}
