#pragma once

// Векторная алгебра в трёхмерном пространстве.
// Модуль ничего не знает ни про Windows, ни про химию — это чистая математика.
// Порядок арифметических действий сохранён как в оригинале на TypeScript:
// у вырожденных минимумов энергии разница в последнем бите способна выбрать
// другую конфигурацию, и тогда сверка с эталоном разойдётся.

#include <array>
#include <cmath>
#include <vector>

namespace chem {

/** Пи. Своя константа: M_PI при -std=c++17 (строгий ANSI) может отсутствовать. */
constexpr double PI = 3.14159265358979323846;

/**
 * Округление ровно как Math.round в JavaScript: половинка уходит ВВЕРХ,
 * а не «от нуля», как у std::round. В расчёте неподелённых пар на этом
 * держится отбрасывание «половинки» ароматической пары.
 */
inline double jsRound(double x) { return std::floor(x + 0.5); }

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Vec3 v3(double x = 0.0, double y = 0.0, double z = 0.0) { return Vec3{x, y, z}; }

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(const Vec3& a, double k) { return Vec3{a.x * k, a.y * k, a.z * k}; }
inline Vec3 operator*(double k, const Vec3& a) { return a * k; }
inline Vec3 operator-(const Vec3& a) { return Vec3{-a.x, -a.y, -a.z}; }
inline Vec3& operator+=(Vec3& a, const Vec3& b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
inline Vec3& operator-=(Vec3& a, const Vec3& b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }

// Имена add/sub/mul оставлены рядом с операторами: так проще сверять
// спорные места с оригиналом строка в строку.
inline Vec3 add(const Vec3& a, const Vec3& b) { return a + b; }
inline Vec3 sub(const Vec3& a, const Vec3& b) { return a - b; }
inline Vec3 mul(const Vec3& a, double k) { return a * k; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

inline double len(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline double dist(const Vec3& a, const Vec3& b) { return len(sub(a, b)); }

/** Нормировка. Нулевой вектор превращается в (0,0,1) — как в оригинале. */
Vec3 norm(const Vec3& a);

/** Любой единичный вектор, перпендикулярный данному. */
Vec3 anyPerp(const Vec3& a);

/** Угол между векторами в градусах. */
double angleDeg(const Vec3& a, const Vec3& b);

/** Валентный угол i–c–j (в градусах) по трём точкам. */
double bondAngle(const Vec3& pi, const Vec3& pc, const Vec3& pj);

// ---------------------------------------------------------------------------
// Матрицы 3×3 (по строкам: m[строка * 3 + столбец])
// ---------------------------------------------------------------------------

using Mat3 = std::array<double, 9>;

extern const Mat3 IDENTITY;

Vec3 applyMat(const Mat3& m, const Vec3& v);
Mat3 matMul(const Mat3& a, const Mat3& b);
Mat3 transpose(const Mat3& m);

/** Матрица поворота вокруг произвольной оси на угол (радианы). Формула Родрига. */
Mat3 rotationAxis(const Vec3& axis, double angle);

/** Минимальный поворот, переводящий единичный вектор from в to. */
Mat3 rotationFromTo(const Vec3& from, const Vec3& to);

/**
 * Ортонормированная система координат по двум векторам: первый задаёт ось X,
 * второй — плоскость XY. Столбцы матрицы — базисные векторы.
 */
Mat3 frameFrom(const Vec3& a, const Vec3& b);

/** Поворот, наилучшим образом совмещающий пару направлений (s1,s2) с (t1,t2). */
Mat3 alignPair(const Vec3& s1, const Vec3& s2, const Vec3& t1, const Vec3& t2);

/** Геометрический центр набора точек. */
Vec3 centroid(const std::vector<Vec3>& points);

struct Eigen3 {
    std::array<double, 3> values;
    std::array<Vec3, 3> vectors;
};

/**
 * Собственные векторы симметричной матрицы 3×3 методом Якоби.
 * Нужны, чтобы развернуть молекулу к зрителю самой широкой стороной.
 */
Eigen3 jacobiEigen(const Mat3& inputMatrix);

}  // namespace chem
