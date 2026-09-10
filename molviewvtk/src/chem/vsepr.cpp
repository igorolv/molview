#include "chem/vsepr.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <unordered_map>

#include "chem/periodic.h"

namespace chem {

namespace {

/** Показатель степени в законе отталкивания E ~ 1/r^n. */
constexpr int REPULSION_EXPONENT = 6;

double repulsionEnergy(const std::vector<Vec3>& points, const std::vector<double>& weights) {
    double e = 0;
    for (std::size_t i = 0; i < points.size(); i++) {
        for (std::size_t j = i + 1; j < points.size(); j++) {
            const double d = std::max(1e-6, len(sub(points[i], points[j])));
            e += (weights[i] * weights[j]) / std::pow(d, REPULSION_EXPONENT);
        }
    }
    return e;
}

/**
 * Детерминированный генератор — одна молекула всегда строится одинаково.
 * Линейный конгруэнтный, переполнение uint32_t повторяет «>>> 0» в JavaScript.
 */
class Random {
public:
    explicit Random(std::uint32_t seed) : state(seed) {}
    double next() {
        state = state * 1664525u + 1013904223u;
        return static_cast<double>(state) / 4294967296.0;
    }

private:
    std::uint32_t state;
};

void relaxOnSphere(std::vector<Vec3>& points, const std::vector<double>& weights,
                   int steps, double step0, double decay) {
    const std::size_t n = points.size();
    double step = step0;
    for (int iter = 0; iter < steps; iter++) {
        std::vector<Vec3> forces(n, v3());
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = i + 1; j < n; j++) {
                const Vec3 diff = sub(points[i], points[j]);
                const double d = std::max(1e-6, len(diff));
                const double magnitude =
                    (REPULSION_EXPONENT * weights[i] * weights[j]) / std::pow(d, REPULSION_EXPONENT + 2);
                const Vec3 f = mul(diff, magnitude);
                forces[i] = add(forces[i], f);
                forces[j] = sub(forces[j], f);
            }
        }
        for (std::size_t i = 0; i < n; i++) {
            // сохраняем только касательную составляющую: точка обязана остаться на сфере
            const Vec3 p = points[i];
            const Vec3 f = forces[i];
            const Vec3 tangential = sub(f, mul(p, dot(f, p)));
            points[i] = norm(add(p, mul(tangential, step)));
        }
        step *= decay;
    }
}

/** Точные шаблоны геометрии для стерических чисел 5–7. */
std::vector<Vec3> idealTemplate(int n) {
    const double s32 = std::sqrt(3.0) / 2.0;
    switch (n) {
        case 5:  // тригональная бипирамида: 3 экваториальные + 2 аксиальные
            return {
                v3(1, 0, 0), v3(-0.5, s32, 0), v3(-0.5, -s32, 0),
                v3(0, 0, 1), v3(0, 0, -1),
            };
        case 6:  // октаэдр
            return {
                v3(0, 0, 1), v3(0, 0, -1),
                v3(1, 0, 0), v3(-1, 0, 0),
                v3(0, 1, 0), v3(0, -1, 0),
            };
        case 7: {  // пентагональная бипирамида
            std::vector<Vec3> result;
            for (int k = 0; k < 5; k++) {
                const double a = (2 * PI * k) / 5;
                result.push_back(v3(std::cos(a), std::sin(a), 0));
            }
            result.push_back(v3(0, 0, 1));
            result.push_back(v3(0, 0, -1));
            return result;
        }
        default:
            return {};
    }
}

/** Все перестановки индексов 0..n-1 (n ≤ 7, поэтому это дёшево). */
std::vector<std::vector<int>> permutations(int n) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    std::vector<bool> used(static_cast<std::size_t>(n), false);

    // рекурсия обычным замыканием: порядок перебора обязан совпасть с оригиналом,
    // потому что при равной энергии выигрывает первая найденная перестановка
    std::function<void()> walk = [&]() {
        if (static_cast<int>(current.size()) == n) { result.push_back(current); return; }
        for (int i = 0; i < n; i++) {
            if (used[i]) continue;
            used[i] = true;
            current.push_back(i);
            walk();
            current.pop_back();
            used[i] = false;
        }
    };
    walk();
    return result;
}

const std::vector<std::vector<int>>& cachedPermutations(int n) {
    static std::map<int, std::vector<std::vector<int>>> cache;
    auto it = cache.find(n);
    if (it == cache.end()) it = cache.emplace(n, permutations(n)).first;
    return it->second;
}

/** Свободная минимизация энергии на сфере (стерическое число ≤ 4). */
std::vector<Vec3> minimizeFreely(const std::vector<double>& weights) {
    const int n = static_cast<int>(weights.size());
    std::vector<Vec3> best;
    double bestEnergy = std::numeric_limits<double>::infinity();

    for (int attempt = 0; attempt < 24; attempt++) {
        Random rnd(static_cast<std::uint32_t>(0x9e3779b9u + attempt * 7919 + n * 104729));
        std::vector<Vec3> points;
        for (int i = 0; i < n; i++) {
            const double u = rnd.next() * 2 - 1;
            const double phi = rnd.next() * PI * 2;
            const double r = std::sqrt(std::max(0.0, 1 - u * u));
            points.push_back(v3(r * std::cos(phi), r * std::sin(phi), u));
        }
        relaxOnSphere(points, weights, 3000, 0.3, 0.9985);
        relaxOnSphere(points, weights, 600, 0.02, 0.999);
        const double e = repulsionEnergy(points, weights);
        if (e < bestEnergy - 1e-10) { bestEnergy = e; best = points; }
    }
    return best;
}

/**
 * Стерическое число 5–7: берём точный шаблон геометрии, но РАСПРЕДЕЛЕНИЕ групп
 * по его вершинам выбираем по минимуму энергии отталкивания. Именно так
 * программа выводит, что неподелённая пара идёт в экваториальную позицию
 * тригональной бипирамиды, а две пары в октаэдре становятся друг напротив друга.
 */
std::vector<Vec3> fitToTemplate(const std::vector<double>& weights,
                                const std::vector<Vec3>& tmpl) {
    const int n = static_cast<int>(weights.size());
    const std::vector<int>* bestPerm = nullptr;
    double bestEnergy = std::numeric_limits<double>::infinity();

    for (const std::vector<int>& perm : cachedPermutations(n)) {
        std::vector<Vec3> points;
        points.reserve(perm.size());
        for (int slot : perm) points.push_back(tmpl[slot]);
        const double e = repulsionEnergy(points, weights);
        if (e < bestEnergy - 1e-12) { bestEnergy = e; bestPerm = &perm; }
    }

    std::vector<Vec3> result;
    for (int slot : *bestPerm) result.push_back(tmpl[slot]);
    return result;
}

/**
 * Приведение к воспроизводимой ориентации: первая группа смотрит вдоль +Z,
 * вторая ложится в плоскость XZ. Иначе одна и та же молекула при каждом
 * запуске выглядела бы повёрнутой по-разному.
 */
void canonicalize(std::vector<Vec3>& points) {
    if (points.size() < 2) return;

    const Vec3 first = norm(points[0]);
    const Vec3 target = v3(0, 0, 1);
    const Vec3 axis = v3(
        first.y * target.z - first.z * target.y,
        first.z * target.x - first.x * target.z,
        first.x * target.y - first.y * target.x);
    const double axisLen = len(axis);
    const double cosA = dot(first, target);

    if (axisLen > 1e-9) {
        const Vec3 a = mul(axis, 1.0 / axisLen);
        const double angle = std::acos(std::max(-1.0, std::min(1.0, cosA)));
        const double c = std::cos(angle), s = std::sin(angle), t = 1 - c;
        for (Vec3& p : points) {
            p = v3(
                (t * a.x * a.x + c) * p.x + (t * a.x * a.y - s * a.z) * p.y + (t * a.x * a.z + s * a.y) * p.z,
                (t * a.x * a.y + s * a.z) * p.x + (t * a.y * a.y + c) * p.y + (t * a.y * a.z - s * a.x) * p.z,
                (t * a.x * a.z - s * a.y) * p.x + (t * a.y * a.z + s * a.x) * p.y + (t * a.z * a.z + c) * p.z);
        }
    } else if (cosA < 0) {
        for (Vec3& p : points) p = mul(p, -1);
    }
    points[0] = v3(0, 0, 1);

    const Vec3 second = points[1];
    const double planar = std::hypot(second.x, second.y);
    if (planar > 1e-9) {
        const double c = second.x / planar;
        const double s = second.y / planar;
        for (Vec3& p : points) {
            p = v3(c * p.x + s * p.y, -s * p.x + c * p.y, p.z);
        }
    }
}

/** Ключ кэша: веса, отформатированные с четырьмя знаками, через разделитель. */
std::string cacheKey(const std::vector<double>& weights) {
    std::string key;
    char buffer[32];
    for (std::size_t i = 0; i < weights.size(); i++) {
        if (i > 0) key += '|';
        std::snprintf(buffer, sizeof(buffer), "%.4f", weights[i]);
        key += buffer;
    }
    return key;
}

}  // namespace

double bondWeight(double order) {
    return WEIGHT_SINGLE_BOND + WEIGHT_PER_EXTRA_ORDER * std::max(0.0, order - 1.0);
}

int lonePairCount(const std::string& el, double charge, double bondOrderSum) {
    const double ve = element(el).ve;
    const double raw = (ve - charge - bondOrderSum) / 2.0;
    return static_cast<int>(std::max(0.0, jsRound(raw - 1e-3)));
}

// ---------------------------------------------------------------------------
// Геометрия по числу связей и неподелённых пар
// ---------------------------------------------------------------------------

const GeometryInfo& geometryFor(int sigma, int lonePairs) {
    static const GeometryInfo fallback = {
        "не определена", "не определена", 0.0, false,
        "Стерическое число выходит за рамки модели ОЭПВО.", "",
    };

    // Ключ: «число σ-связей : число неподелённых пар».
    static const std::map<std::string, GeometryInfo> table = {
        {"1:0", {"одна связь", "концевой атом", 0, false, "Единственный сосед — валентного угла нет.", "водород в H₂O и CH₄"}},
        {"1:1", {"линейная", "концевой атом", 0, false, "Одна связь и одна неподелённая пара.", "азот в HC≡N, углерод в CO"}},
        {"1:2", {"тригональная", "концевой атом", 0, false, "Одна связь и две неподелённые пары.", "концевой кислород двойной связи — в C=O, SO₂, NO₃⁻"}},
        {"1:3", {"тетраэдрическая", "концевой атом", 0, false, "Одна связь и три неподелённые пары.", "галогены в HCl и CCl₄, кислород в OH⁻ и NO₃⁻"}},

        {"2:0", {"линейная", "линейная", 180, true, "Две связи без пар — они расходятся на 180°.", "CO₂, BeCl₂, ацетилен"}},
        {"2:1", {"тригональная", "угловая", 120, true, "Две связи и одна пара — частица угловая.", "SO₂, озон, нитрит-ион"}},
        {"2:2", {"тетраэдрическая", "угловая", 109.5, true, "Две связи и две пары — частица угловая.", "вода, сероводород, спирты и эфиры"}},
        {"2:3", {"тригонально-бипирамидальная", "линейная", 180, true, "Две связи и три пары: пары в экваторе, связи на оси.", "XeF₂"}},

        {"3:0", {"тригональная", "плоская треугольная", 120, true, "Три связи без пар — плоский треугольник.", "BF₃, SO₃, нитрат-ион, любой sp²-углерод"}},
        {"3:1", {"тетраэдрическая", "тригонально-пирамидальная", 109.5, true, "Три связи и одна пара — пирамида.", "аммиак, PCl₃, ион гидроксония, сера в ДМСО"}},
        {"3:2", {"тригонально-бипирамидальная", "Т-образная", 90, true, "Три связи и две пары: обе пары в экваторе, форма «Т».", "ClF₃"}},

        {"4:0", {"тетраэдрическая", "тетраэдрическая", 109.5, true, "Четыре связи без пар — правильный тетраэдр.", "метан, ион аммония, сульфат-ион"}},
        {"4:1", {"тригонально-бипирамидальная", "качели (дисфеноид)", 120, true, "Четыре связи и одна пара: пара в экваторе, форма качелей.", "SF₄"}},
        {"4:2", {"октаэдрическая", "квадратная", 90, true, "Четыре связи и две пары: пары транс, связи в квадрате.", "XeF₄"}},

        {"5:0", {"тригонально-бипирамидальная", "тригонально-бипирамидальная", 120, true, "Пять связей: 3 экваториальные и 2 аксиальные.", "PCl₅"}},
        {"5:1", {"октаэдрическая", "квадратно-пирамидальная", 90, true, "Пять связей и одна пара: квадратная пирамида.", "IF₅, BrF₅"}},

        {"6:0", {"октаэдрическая", "октаэдрическая", 90, true, "Шесть равноценных связей — правильный октаэдр.", "SF₆"}},

        {"7:0", {"пентагонально-бипирамидальная", "пентагонально-бипирамидальная", 72, true, "Семь связей — редкая геометрия.", "IF₇"}},
    };

    const auto it = table.find(std::to_string(sigma) + ":" + std::to_string(lonePairs));
    return it == table.end() ? fallback : it->second;
}

std::string hybridFor(int steric) {
    switch (steric) {
        case 1: return "—";
        case 2: return "sp";
        case 3: return "sp²";
        case 4: return "sp³";
        case 5: return "sp³d";
        case 6: return "sp³d²";
        case 7: return "sp³d³";
        default: return "—";
    }
}

std::string hybridExplanation(int steric) {
    switch (steric) {
        case 2: return "одна s- и одна p-орбиталь дают две sp-орбитали под 180°";
        case 3: return "одна s- и две p-орбитали дают три sp²-орбитали под 120° в одной плоскости";
        case 4: return "одна s- и три p-орбитали дают четыре sp³-орбитали, направленные к вершинам тетраэдра";
        case 5: return "к s- и p-орбиталям добавляется одна d-орбиталь: пять sp³d-орбиталей образуют тригональную бипирамиду";
        case 6: return "к s- и p-орбиталям добавляются две d-орбитали: шесть sp³d²-орбиталей образуют октаэдр";
        case 7: return "семь sp³d³-орбиталей образуют пентагональную бипирамиду";
        default: return "";
    }
}

std::string axeNotation(int sigma, int lonePairs) {
    static const char* const SUB[] = {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};
    const auto subscript = [](int n) {
        std::string digits = std::to_string(n);
        std::string out;
        for (char c : digits) out += SUB[c - '0'];
        return out;
    };

    std::string s = "A";
    if (sigma > 0) s += "X" + (sigma > 1 ? subscript(sigma) : std::string());
    if (lonePairs > 0) s += "E" + (lonePairs > 1 ? subscript(lonePairs) : std::string());
    return s;
}

std::vector<Vec3> arrangeGroups(const std::vector<double>& weights) {
    const int n = static_cast<int>(weights.size());
    if (n == 0) return {};
    if (n == 1) return {v3(0, 0, 1)};

    // Без кэша минимизация гонялась бы заново для каждого атома молекулы.
    static std::unordered_map<std::string, std::vector<Vec3>> directionCache;
    const std::string key = cacheKey(weights);
    const auto cached = directionCache.find(key);
    if (cached != directionCache.end()) return cached->second;

    std::vector<Vec3> result;
    if (n <= 4) {
        result = minimizeFreely(weights);
    } else if (n <= 7) {
        result = fitToTemplate(weights, idealTemplate(n));
    } else {
        result = minimizeFreely(weights);
    }

    canonicalize(result);
    directionCache.emplace(key, result);
    return result;
}

}  // namespace chem
