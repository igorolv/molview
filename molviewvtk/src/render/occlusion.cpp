#include "render/occlusion.h"

#include <cmath>

namespace render {

using chem::dot;
using chem::mul;
using chem::sub;
using chem::Vec3;

bool Occluder::visible(const Vec3& eye, const Vec3& point, int ignore) const {
    const Vec3 ray = sub(point, eye);
    const double rayLen2 = dot(ray, ray);
    if (rayLen2 < 1e-12) return true;

    // Параметр t вдоль луча: 0 — глаз, 1 — проверяемая точка. Отступы с обоих
    // концов убирают два случая: касание в самой точке и шар, оказавшийся
    // позади неё.
    const double eps = 1e-3;

    for (std::size_t i = 0; i < spheres.size(); i++) {
        if (static_cast<int>(i) == ignore) continue;
        const Sphere& s = spheres[i];

        const Vec3 toCenter = sub(s.center, eye);
        // t ближайшего к центру шара места на луче
        const double t = dot(toCenter, ray) / rayLen2;
        if (t <= eps || t >= 1 - eps) continue;

        // квадрат расстояния от центра шара до луча
        const Vec3 nearest = chem::add(eye, mul(ray, t));
        const Vec3 offset = sub(s.center, nearest);
        if (dot(offset, offset) < s.radius * s.radius) return false;
    }
    return true;
}

bool Occluder::rectVisible(const Vec3& eye, const Vec3& center, double halfWidth, double halfHeight,
                           const Vec3& right, const Vec3& up, int ignore) const {
    if (!visible(eye, center, ignore)) return false;
    for (int sx = -1; sx <= 1; sx += 2) {
        for (int sy = -1; sy <= 1; sy += 2) {
            const Vec3 corner = chem::add(center, chem::add(mul(right, sx * halfWidth), mul(up, sy * halfHeight)));
            if (!visible(eye, corner, ignore)) return false;
        }
    }
    return true;
}

}  // namespace render
