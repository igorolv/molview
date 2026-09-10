#include "chem/vec.h"

#include <algorithm>

namespace chem {

const Mat3 IDENTITY = {1, 0, 0, 0, 1, 0, 0, 0, 1};

namespace {
/** Ограничение косинуса: из-за ошибок округления он выходит за [-1, 1]. */
double clampCos(double c) { return std::max(-1.0, std::min(1.0, c)); }
}  // namespace

Vec3 norm(const Vec3& a) {
    const double l = len(a);
    return l < 1e-12 ? v3(0, 0, 1) : mul(a, 1.0 / l);
}

Vec3 anyPerp(const Vec3& a) {
    const Vec3 n = norm(a);
    const Vec3 helper = std::fabs(n.x) < 0.9 ? v3(1, 0, 0) : v3(0, 1, 0);
    return norm(cross(n, helper));
}

double angleDeg(const Vec3& a, const Vec3& b) {
    const double c = dot(norm(a), norm(b));
    return std::acos(clampCos(c)) * 180.0 / PI;
}

double bondAngle(const Vec3& pi, const Vec3& pc, const Vec3& pj) {
    return angleDeg(sub(pi, pc), sub(pj, pc));
}

Vec3 applyMat(const Mat3& m, const Vec3& v) {
    return Vec3{
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[3] * v.x + m[4] * v.y + m[5] * v.z,
        m[6] * v.x + m[7] * v.y + m[8] * v.z,
    };
}

Mat3 matMul(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double s = 0;
            for (int k = 0; k < 3; k++) s += a[i * 3 + k] * b[k * 3 + j];
            r[i * 3 + j] = s;
        }
    }
    return r;
}

Mat3 transpose(const Mat3& m) {
    return Mat3{m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]};
}

Mat3 rotationAxis(const Vec3& axis, double angle) {
    const Vec3 a = norm(axis);
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double t = 1 - c;
    const double x = a.x, y = a.y, z = a.z;
    return Mat3{
        t * x * x + c,     t * x * y - s * z, t * x * z + s * y,
        t * x * y + s * z, t * y * y + c,     t * y * z - s * x,
        t * x * z - s * y, t * y * z + s * x, t * z * z + c,
    };
}

Mat3 rotationFromTo(const Vec3& from, const Vec3& to) {
    const Vec3 a = norm(from);
    const Vec3 b = norm(to);
    const double d = dot(a, b);
    if (d > 0.999999) return IDENTITY;
    if (d < -0.999999) return rotationAxis(anyPerp(a), PI);
    return rotationAxis(cross(a, b), std::acos(clampCos(d)));
}

Mat3 frameFrom(const Vec3& a, const Vec3& b) {
    const Vec3 e1 = norm(a);
    Vec3 e2 = sub(b, mul(e1, dot(b, e1)));
    if (len(e2) < 1e-8) e2 = anyPerp(e1);
    e2 = norm(e2);
    const Vec3 e3 = cross(e1, e2);
    // столбцы: e1, e2, e3
    return Mat3{e1.x, e2.x, e3.x, e1.y, e2.y, e3.y, e1.z, e2.z, e3.z};
}

Mat3 alignPair(const Vec3& s1, const Vec3& s2, const Vec3& t1, const Vec3& t2) {
    return matMul(frameFrom(t1, t2), transpose(frameFrom(s1, s2)));
}

Vec3 centroid(const std::vector<Vec3>& points) {
    if (points.empty()) return v3();
    Vec3 s = v3();
    for (const Vec3& p : points) s = add(s, p);
    return mul(s, 1.0 / static_cast<double>(points.size()));
}

Eigen3 jacobiEigen(const Mat3& inputMatrix) {
    Mat3 a = inputMatrix;
    Mat3 v = IDENTITY;

    for (int sweep = 0; sweep < 64; sweep++) {
        double off = 0;
        for (int i = 0; i < 3; i++)
            for (int j = i + 1; j < 3; j++) off += a[i * 3 + j] * a[i * 3 + j];
        if (off < 1e-18) break;

        for (int p = 0; p < 3; p++) {
            for (int q = p + 1; q < 3; q++) {
                const double apq = a[p * 3 + q];
                if (std::fabs(apq) < 1e-18) continue;
                const double theta = (a[q * 3 + q] - a[p * 3 + p]) / (2 * apq);
                // Math.sign(theta || 1): у нулевого theta знак берётся положительным
                const double sign = theta > 0 ? 1.0 : (theta < 0 ? -1.0 : 1.0);
                const double t = sign / (std::fabs(theta) + std::sqrt(theta * theta + 1));
                const double c = 1 / std::sqrt(t * t + 1);
                const double s = t * c;

                for (int k = 0; k < 3; k++) {
                    const double akp = a[k * 3 + p];
                    const double akq = a[k * 3 + q];
                    a[k * 3 + p] = c * akp - s * akq;
                    a[k * 3 + q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; k++) {
                    const double apk = a[p * 3 + k];
                    const double aqk = a[q * 3 + k];
                    a[p * 3 + k] = c * apk - s * aqk;
                    a[q * 3 + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; k++) {
                    const double vkp = v[k * 3 + p];
                    const double vkq = v[k * 3 + q];
                    v[k * 3 + p] = c * vkp - s * vkq;
                    v[k * 3 + q] = s * vkp + c * vkq;
                }
            }
        }
    }

    Eigen3 result;
    result.values = {a[0], a[4], a[8]};
    result.vectors = {
        v3(v[0], v[3], v[6]),
        v3(v[1], v[4], v[7]),
        v3(v[2], v[5], v[8]),
    };
    return result;
}

}  // namespace chem
