#include "chem/analyze.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "chem/build3d.h"
#include "chem/formula.h"
#include "chem/periodic.h"
#include "chem/resonance.h"
#include "chem/rings.h"
#include "chem/smiles.h"
#include "chem/vsepr.h"

namespace chem {

namespace {

/** Масштаб перевода качественной оценки дипольного момента в дебаи. */
constexpr double DIPOLE_SCALE = 1.25;
/** Вклад одной неподелённой пары в дипольный момент. */
constexpr double DIPOLE_LONE_PAIR = 0.50;

Vec3 normWorld(const Mat3& frame, const Vec3& local) {
    return norm(v3(
        frame[0] * local.x + frame[1] * local.y + frame[2] * local.z,
        frame[3] * local.x + frame[4] * local.y + frame[5] * local.z,
        frame[6] * local.x + frame[7] * local.y + frame[8] * local.z));
}

double angleBetween(const Vec3& a, const Vec3& b) {
    return std::acos(std::max(-1.0, std::min(1.0, dot(norm(a), norm(b))))) * 180.0 / PI;
}

/** Наиболее характерный предсказанный угол между связями данного атома. */
bool predictedAngleFor(const AtomEnvironment& e, double& out) {
    if (e.neighbors.size() < 2) return false;
    std::vector<double> values;
    for (std::size_t a = 0; a < e.neighbors.size(); a++) {
        for (std::size_t b = a + 1; b < e.neighbors.size(); b++) {
            values.push_back(angleBetween(e.localDirs[a], e.localDirs[b]));
        }
    }
    std::sort(values.begin(), values.end());
    // берём наименьший — именно его обычно и приводят в справочниках
    out = values[0];
    return true;
}

/**
 * Предупреждения о случаях, где модель ОЭПВО заведомо расходится с опытом.
 * Это не ошибка программы, а известная граница применимости теории.
 */
std::string warningFor(const std::string& el, int steric, int lonePairs, bool inRing) {
    const int period = element(el).period;

    if (period >= 3 && lonePairs > 0 && steric == 4 && el != "Cl" && el != "Br" && el != "I") {
        return "У элемента 3-го периода и ниже s- и p-орбитали сильно различаются по энергии, "
               "гибридизация выражена слабо. Реальный угол заметно ближе к 90°, чем предсказывает модель.";
    }
    if (inRing && steric == 4) {
        return "Атом входит в цикл: валентный угол задан геометрией кольца и может сильно "
               "отличаться от предсказанного. Разница и есть угловое напряжение цикла.";
    }
    if (steric >= 5 && lonePairs > 0) {
        return "Для гипервалентных частиц модель верно определяет ФОРМУ, но тонкое искажение "
               "углов неподелёнными парами она передаёт лишь качественно.";
    }
    return "";
}

/** Приведение экспериментальных подписей к тому же виду, что и вычисленные. */
std::vector<std::pair<std::string, std::vector<double>>> normalizeExperimental(
    const MoleculeMeta* meta) {
    std::vector<std::pair<std::string, std::vector<double>>> table;
    if (meta == nullptr) return table;

    for (const auto& item : meta->exp) {
        // подпись «H-O-H» разбирается на три символа и нормализуется
        std::vector<std::string> parts;
        std::string current;
        for (char c : item.first) {
            if (c == '-') { parts.push_back(current); current.clear(); }
            else current += c;
        }
        parts.push_back(current);

        const std::string key = parts.size() == 3
            ? angleLabel(parts[0], parts[1], parts[2])
            : item.first;
        table.emplace_back(key, item.second);
    }
    return table;
}

/**
 * Привязка экспериментальных значений к конкретным углам.
 *
 * Подпись «F-S-F» в гексафториде серы относится сразу к пятнадцати углам:
 * двенадцати по 90° и трём по 180°. Поэтому каждое справочное значение
 * привязывается к тому ещё не занятому углу, предсказание для которого
 * к нему ближе всего.
 */
void attachExperimental(std::vector<AngleRecord>& angles,
                        const std::vector<std::pair<std::string, std::vector<double>>>& table) {
    for (const auto& row : table) {
        std::vector<AngleRecord*> candidates;
        for (AngleRecord& a : angles) {
            if (a.label == row.first) candidates.push_back(&a);
        }
        if (candidates.empty()) continue;

        for (double value : row.second) {
            AngleRecord* best = nullptr;
            double bestDiff = std::numeric_limits<double>::infinity();
            for (AngleRecord* c : candidates) {
                if (c->hasExperimental) continue;
                const double diff = std::fabs((c->hasPredicted ? c->predicted : 0.0) - value);
                if (diff < bestDiff) { bestDiff = diff; best = c; }
            }
            if (best != nullptr) { best->experimental = value; best->hasExperimental = true; }
        }
    }
}

}  // namespace

Analysis analyze(Molecule mol) {
    // Усредняем кратности равноценных связей: нитрат-ион должен получиться
    // правильным треугольником, а не «двойная плюс две одинарные».
    applyResonance(mol);

    const std::vector<std::vector<int>> nb = neighborLists(mol);
    const std::size_t n = mol.atoms.size();

    // --- сумма кратностей связей у каждого атома ------------------------------
    std::vector<double> orderSum(n, 0.0);
    std::map<std::pair<int, int>, double> bondOrderTo;
    for (const Bond& b : mol.bonds) {
        orderSum[b.a] += b.order;
        orderSum[b.b] += b.order;
        bondOrderTo[{b.a, b.b}] = b.order;
        bondOrderTo[{b.b, b.a}] = b.order;
    }

    // --- окружение каждого атома по теории Гиллеспи ---------------------------
    std::vector<AtomEnvironment> env;
    env.reserve(n);
    for (std::size_t id = 0; id < n; id++) {
        const Atom& atom = mol.atoms[id];
        AtomEnvironment e;
        e.neighbors = nb[id];
        for (int k : e.neighbors) {
            const auto it = bondOrderTo.find({static_cast<int>(id), k});
            e.bondOrders.push_back(it == bondOrderTo.end() ? 1.0 : it->second);
        }
        e.lonePairs = lonePairCount(atom.el, atom.charge, orderSum[id]);
        e.steric = static_cast<int>(e.neighbors.size()) + e.lonePairs;

        std::vector<double> weights;
        for (double order : e.bondOrders) weights.push_back(bondWeight(order));
        for (int k = 0; k < e.lonePairs; k++) weights.push_back(WEIGHT_LONE_PAIR);

        e.localDirs = arrangeGroups(weights);
        env.push_back(std::move(e));
    }

    const std::vector<Ring> rings = findRings(mol);
    const BuildResult build = buildGeometry(mol, env, rings);

    // --- разбор каждого атома -------------------------------------------------
    std::set<int> ringAtoms;
    for (const Ring& r : rings) {
        for (int a : r.atoms) ringAtoms.insert(a);
    }

    std::vector<AtomAnalysis> atoms;
    for (std::size_t id = 0; id < n; id++) {
        const Atom& atom = mol.atoms[id];
        const AtomEnvironment& e = env[id];
        const int sigma = static_cast<int>(e.neighbors.size());
        const GeometryInfo& geom = geometryFor(sigma, e.lonePairs);

        AtomAnalysis a;
        a.id = static_cast<int>(id);
        a.el = atom.el;
        a.charge = atom.charge;
        a.sigma = sigma;
        a.lonePairs = e.lonePairs;
        a.steric = e.steric;
        a.axe = axeNotation(sigma, e.lonePairs);
        a.hybrid = hybridFor(e.steric);
        a.electronGeom = geom.electronGeom;
        a.molecularGeom = geom.molecularGeom;
        a.idealAngle = geom.idealAngle;
        a.hasIdealAngle = geom.hasIdealAngle;
        a.hasPredictedAngle = predictedAngleFor(e, a.predictedAngle);
        a.isCentral = sigma >= 2;
        a.lonePairDirs = build.lonePairDirs[id];
        for (const Vec3& d : e.localDirs) a.orbitalDirs.push_back(normWorld(build.frames[id], d));
        a.warning = warningFor(atom.el, e.steric, e.lonePairs, ringAtoms.count(static_cast<int>(id)) != 0);
        atoms.push_back(std::move(a));
    }

    // --- валентные углы -------------------------------------------------------
    std::vector<AngleRecord> angles;
    for (std::size_t id = 0; id < n; id++) {
        const AtomEnvironment& e = env[id];
        if (e.neighbors.size() < 2) continue;
        for (std::size_t a = 0; a < e.neighbors.size(); a++) {
            for (std::size_t b = a + 1; b < e.neighbors.size(); b++) {
                const int i = e.neighbors[a];
                const int j = e.neighbors[b];
                AngleRecord rec;
                rec.i = i;
                rec.center = static_cast<int>(id);
                rec.j = j;
                rec.label = angleLabel(mol.atoms[i].el, mol.atoms[id].el, mol.atoms[j].el);
                rec.actual = bondAngle(mol.atoms[i].pos, mol.atoms[id].pos, mol.atoms[j].pos);
                rec.predicted = angleBetween(e.localDirs[a], e.localDirs[b]);
                rec.hasPredicted = true;
                rec.inRing = ringAtoms.count(static_cast<int>(id)) != 0
                          && ringAtoms.count(i) != 0 && ringAtoms.count(j) != 0;
                angles.push_back(std::move(rec));
            }
        }
    }
    attachExperimental(angles, normalizeExperimental(mol.meta));

    // --- дипольный момент (качественная оценка) -------------------------------
    Vec3 dipoleVec = v3();
    for (const Bond& b : mol.bonds) {
        const double enA = element(mol.atoms[b.a].el).en;
        const double enB = element(mol.atoms[b.b].el).en;
        const Vec3 direction = norm(sub(mol.atoms[b.b].pos, mol.atoms[b.a].pos));
        dipoleVec = add(dipoleVec, mul(direction, enB - enA));
    }
    for (std::size_t id = 0; id < n; id++) {
        for (const Vec3& d : build.lonePairDirs[id]) {
            dipoleVec = add(dipoleVec, mul(d, DIPOLE_LONE_PAIR));
        }
    }
    dipoleVec = mul(dipoleVec, DIPOLE_SCALE);
    const double dipoleValue = len(dipoleVec);

    // --- сводка по гибридизациям ----------------------------------------------
    std::vector<std::pair<std::string, int>> hybridSummary;
    for (const AtomAnalysis& a : atoms) {
        // концевые атомы в сводку не попадают: интересна гибридизация центров
        if (a.sigma < 2 || a.hybrid == "—") continue;
        auto it = hybridSummary.begin();
        while (it != hybridSummary.end() && it->first != a.hybrid) ++it;
        if (it == hybridSummary.end()) hybridSummary.emplace_back(a.hybrid, 1);
        else it->second++;
    }

    Analysis result;
    result.formula = plainFormula(mol);
    result.formulaHtml = prettyFormula(mol);
    result.mass = molecularMass(mol);
    result.atoms = std::move(atoms);
    result.angles = std::move(angles);
    result.rings = rings;
    result.dipoleVec = dipoleVec;
    result.dipoleValue = dipoleValue;
    result.polar = dipoleValue > 0.06;
    result.hybridSummary = std::move(hybridSummary);
    result.molecule = std::move(mol);
    return result;
}

}  // namespace chem
