#include "chem/periodic.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "chem/json.h"
#include "chem/vec.h"

namespace chem {

namespace {

/**
 * Стандартные валентности для расчёта числа неявных атомов водорода.
 * Значения соответствуют «органическому подмножеству» нотации SMILES.
 */
const std::map<std::string, std::vector<int>>& valenceTable() {
    static const std::map<std::string, std::vector<int>> table = {
        {"H", {1}}, {"B", {3}}, {"C", {4}}, {"N", {3, 5}}, {"O", {2}},
        {"P", {3, 5}}, {"S", {2, 4, 6}},
        {"F", {1}}, {"Cl", {1}}, {"Br", {1}}, {"I", {1}},
        {"Si", {4}}, {"Ge", {4}}, {"Sn", {4}}, {"As", {3, 5}},
        {"Se", {2, 4, 6}}, {"Te", {2, 4, 6}}, {"Sb", {3, 5}},
        {"Be", {2}}, {"Mg", {2}}, {"Al", {3}}, {"Li", {1}}, {"Na", {1}},
        {"K", {1}}, {"Ca", {2}},
        {"Xe", {0}}, {"Kr", {0}}, {"Ar", {0}}, {"Ne", {0}}, {"He", {0}},
    };
    return table;
}

/** Элементы, у которых положительный заряд УМЕНЬШАЕТ валентность (левая часть таблицы). */
bool isElectronDeficient(const std::string& symbol) {
    static const std::set<std::string> set = {"B", "C", "Si", "Ge", "Sn", "Al", "Be", "Mg"};
    return set.count(symbol) != 0;
}

const ElementInfo UNKNOWN = {
    0, "неизвестный", 14, 3, 4, 2.2, 1.0, 1.7, "#9aa0aa", 12.0,
};

/** Справочник читается с диска один раз при первом обращении. */
const std::map<std::string, ElementInfo>& table() {
    static const std::map<std::string, ElementInfo> loaded = [] {
        std::map<std::string, ElementInfo> result;
        const Json root = Json::parseFile(findDataFile("elements.json"));
        for (const auto& field : root["elements"].fields()) {
            const Json& e = field.second;
            ElementInfo info;
            info.z = static_cast<int>(e["z"].number());
            info.name = e["name"].str();
            info.group = static_cast<int>(e["group"].number());
            info.period = static_cast<int>(e["period"].number());
            info.ve = static_cast<int>(e["ve"].number());
            info.en = e["en"].number();
            info.r = e["r"].number();
            info.vdw = e["vdw"].number();
            info.color = e["color"].str();
            info.mass = e["mass"].number();
            result.emplace(field.first, info);
        }
        return result;
    }();
    return loaded;
}

}  // namespace

const ElementInfo& element(const std::string& symbol) {
    const auto& t = table();
    const auto it = t.find(symbol);
    return it == t.end() ? UNKNOWN : it->second;
}

bool isKnownElement(const std::string& symbol) {
    return table().count(symbol) != 0;
}

int implicitHydrogens(const std::string& symbol, double bondSum, double charge, bool aromatic) {
    const auto it = valenceTable().find(symbol);
    if (it == valenceTable().end()) return 0;
    const std::vector<int>& valences = it->second;
    if (valences.empty() || valences[0] == 0) return 0;

    const double shift = isElectronDeficient(symbol) ? -std::fabs(charge) : charge;

    if (aromatic) {
        return static_cast<int>(std::max(0.0, jsRound(valences[0] + shift - bondSum)));
    }

    double target = valences.back() + shift;
    for (int v : valences) {
        const double adjusted = v + shift;
        if (adjusted >= bondSum - 1e-6) { target = adjusted; break; }
    }
    return static_cast<int>(std::max(0.0, jsRound(target - bondSum)));
}

double bondLength(const std::string& elA, const std::string& elB, double order) {
    const double base = element(elA).r + element(elB).r;
    return base - bondShortening(order);
}

double bondShortening(double order) {
    const double o = std::max(1.0, std::min(3.0, order));
    return o <= 2 ? 0.19 * (o - 1) : 0.19 + 0.13 * (o - 2);
}

}  // namespace chem
