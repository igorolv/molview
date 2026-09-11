#include "chem/orbital.h"

#include <algorithm>
#include <cmath>

namespace chem {

namespace {

/**
 * Коэффициенты при s и p в гибриде. Сумма квадратов равна единице — это и есть
 * условие нормировки гибрида, потому что 2s и 2p ортонормированы между собой.
 */
struct Mix {
    double s;
    double p;
};

Mix mixFor(Orbital kind) {
    switch (kind) {
        case Orbital::S: return {1.0, 0.0};
        case Orbital::P: return {0.0, 1.0};
        case Orbital::Sp: return {1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0)};
        case Orbital::Sp2: return {1.0 / std::sqrt(3.0), std::sqrt(2.0 / 3.0)};
        case Orbital::Sp3: break;
    }
    return {0.5, std::sqrt(3.0) / 2.0};
}

/**
 * Нормировочный множитель 2p при заряде ядра zeff, Å^(−3/2).
 *
 * У водородоподобных 2s и 2p он ОДИН И ТОТ ЖЕ: в R₂₀Y₀₀ и R₂₁Y₁₀ после
 * сокращения √3/√6 = 1/√2 остаётся 1/(4√(2π)) и там, и там. Поэтому задание и
 * записывает обе функции «с точностью до нормировки» одинаково.
 *
 * Множитель 1/BOHR^(3/2) переводит его из атомных единиц в ангстремы: сетка
 * живёт в ангстремах, и без перевода ∫ψ²dV по ней дало бы не единицу, а куб
 * радиуса Бора.
 */
double normalizationP(double zeff) {
    return std::pow(zeff, 1.5) / (4.0 * std::sqrt(2.0 * PI)) / std::pow(BOHR, 1.5);
}

/**
 * Нормировочный множитель слейтеровской 2s, Å^(−3/2).
 *
 * У неё та же радиальная часть r·e^(−r/2), что у 2p, но нет угловой: где у p
 * стоит √(3/4π)·cosθ, у s стоит 1/√(4π). Отсюда ровно √3 разницы, и второй раз
 * интеграл считать не нужно.
 */
double normalizationSlaterS(double zeff) { return normalizationP(zeff) / std::sqrt(3.0); }

}  // namespace

Orbital orbitalForSteric(int steric) {
    if (steric <= 2) return Orbital::Sp;
    if (steric == 3) return Orbital::Sp2;
    return Orbital::Sp3;
}

Vec3 OrbitalGrid::point(int i, int j, int k) const {
    const double s = step();
    return Vec3{center.x - half + s * i, center.y - half + s * j, center.z - half + s * k};
}

double OrbitalGrid::peak() const {
    double best = 0;
    for (double v : values) best = std::max(best, std::abs(v));
    return best;
}

double OrbitalGrid::highest() const {
    double best = 0;
    for (double v : values) best = std::max(best, v);
    return best;
}

double OrbitalGrid::lowest() const {
    double best = 0;
    for (double v : values) best = std::min(best, v);
    return best;
}

Vec3 OrbitalGrid::axis() const {
    Vec3 moment;
    const double s = step();
    std::size_t n = 0;
    for (int k = 0; k < size; k++) {
        const double z = -half + s * k;
        for (int j = 0; j < size; j++) {
            const double y = -half + s * j;
            for (int i = 0; i < size; i++, n++) {
                moment += Vec3{-half + s * i, y, z} * values[n];
            }
        }
    }
    return norm(moment);
}

double overlap(const OrbitalGrid& a, const OrbitalGrid& b) {
    if (a.size != b.size || a.values.size() != b.values.size()) return 0;
    double sum = 0;
    for (std::size_t n = 0; n < a.values.size(); n++) sum += a.values[n] * b.values[n];
    return sum * a.cellVolume();
}

std::vector<OrbitalGrid> orbitalGrids(Orbital kind, const std::vector<Vec3>& dirs,
                                      const Vec3& center, const OrbitalOptions& options) {
    int size = std::max(3, options.size);
    if (size % 2 == 0) size++;
    const double half = options.half;
    const double zeff = options.zeff > 0 ? options.zeff : 1.0;

    const std::size_t total = static_cast<std::size_t>(size) * size * size;
    const double step = 2 * half / (size - 1);
    const Mix mix = mixFor(kind);
    const double scaleP = normalizationP(zeff);
    const double scaleS =
        options.radial == Radial::Slater ? normalizationSlaterS(zeff) : scaleP;

    // Радиальная часть от направления гибрида не зависит, а корень и экспонента
    // в каждом узле — это почти весь расчёт. Считаем её один раз на все
    // направления: sPart — готовая 2s, а pUnit домножается на (d·r), где r
    // берётся в АНГСТРЕМАХ, поэтому перевод длины уже внесён в множитель.
    std::vector<double> sPart(total);
    std::vector<double> pUnit(total);
    std::size_t n = 0;
    for (int k = 0; k < size; k++) {
        const double z = -half + step * k;
        for (int j = 0; j < size; j++) {
            const double y = -half + step * j;
            for (int i = 0; i < size; i++, n++) {
                const double x = -half + step * i;
                const double r = zeff * std::sqrt(x * x + y * y + z * z) / BOHR;
                const double decay = std::exp(-r / 2);
                sPart[n] = scaleS * decay * (options.radial == Radial::Slater ? r : 2 - r);
                pUnit[n] = scaleP * decay * zeff / BOHR;
            }
        }
    }

    std::vector<OrbitalGrid> out;
    out.reserve(dirs.size());
    for (const Vec3& raw : dirs) {
        const Vec3 d = norm(raw);
        OrbitalGrid grid;
        grid.size = size;
        grid.half = half;
        grid.center = center;
        grid.values.resize(total);
        n = 0;
        for (int k = 0; k < size; k++) {
            const double z = -half + step * k;
            for (int j = 0; j < size; j++) {
                const double y = -half + step * j;
                for (int i = 0; i < size; i++, n++) {
                    const double projection = d.x * (-half + step * i) + d.y * y + d.z * z;
                    grid.values[n] = mix.s * sPart[n] + mix.p * pUnit[n] * projection;
                }
            }
        }
        out.push_back(std::move(grid));
    }
    return out;
}

OrbitalGrid orbitalGrid(Orbital kind, const Vec3& dir, const Vec3& center,
                        const OrbitalOptions& options) {
    return orbitalGrids(kind, std::vector<Vec3>{dir}, center, options).front();
}

}  // namespace chem
