/**
 * Самопроверка химического ядра на всей базе молекул.
 * Запуск: make test  (или build/selftest.exe из корня проекта)
 *
 * Проверяется:
 *   • каждая молекула базы разбирается и строится без ошибок;
 *   • длины связей соответствуют сумме ковалентных радиусов;
 *   • атомы не налезают друг на друга;
 *   • предсказанные валентные углы сравниваются с экспериментальными;
 *   • контрольные молекулы имеют ожидаемую форму и гибридизацию;
 *   • все вычисленные величины сходятся с эталоном data/reference.json,
 *     снятым с исходной версии на TypeScript.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "chem/analyze.h"
#include "chem/database.h"
#include "chem/formula.h"
#include "chem/json.h"
#include "chem/periodic.h"
#include "chem/types.h"
#include "chem/vec.h"
#include "testutil.h"

using namespace chem;
using testutil::fixed;
using testutil::pad;
using testutil::padLeft;
using testutil::repeat;

namespace {

testutil::Problems problems;

void say(const std::string& line) { std::printf("%s\n", line.c_str()); }

/** Максимальное отклонение точек от наилучшей плоскости. */
double maxPlaneDeviation(const std::vector<Vec3>& points) {
    if (points.size() < 4) return 0;
    const Vec3 center = centroid(points);

    double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
    for (const Vec3& p : points) {
        const Vec3 d = sub(p, center);
        xx += d.x * d.x; xy += d.x * d.y; xz += d.x * d.z;
        yy += d.y * d.y; yz += d.y * d.z; zz += d.z * d.z;
    }
    // нормаль — собственный вектор с наименьшим собственным значением; ищем грубым перебором
    double best = std::numeric_limits<double>::infinity();
    Vec3 bestNormal = v3(0, 0, 1);
    for (int a = 0; a < 60; a++) {
        for (int b = 0; b < 60; b++) {
            const double theta = (PI * a) / 60;
            const double phi = (2 * PI * b) / 60;
            const double nx = std::sin(theta) * std::cos(phi);
            const double ny = std::sin(theta) * std::sin(phi);
            const double nz = std::cos(theta);
            const double value = xx * nx * nx + yy * ny * ny + zz * nz * nz
                               + 2 * (xy * nx * ny + xz * nx * nz + yz * ny * nz);
            if (value < best) { best = value; bestNormal = v3(nx, ny, nz); }
        }
    }
    double maxDev = 0;
    for (const Vec3& p : points) {
        maxDev = std::max(maxDev, std::fabs(dot(sub(p, center), bestNormal)));
    }
    return maxDev;
}

// ---------------------------------------------------------------------------
// 1. Полный прогон базы
// ---------------------------------------------------------------------------

std::vector<std::pair<std::string, Analysis>> runDatabase() {
    say("");
    say("=== РАЗБОР БАЗЫ ===");
    say("");
    say(pad("Молекула", 24) + pad("Формула", 12) + pad("Гибридизации", 22) + pad("Циклы", 7) + "μ, Д");
    say(repeat("-", 78));

    std::vector<std::pair<std::string, Analysis>> results;

    for (const DbEntry& entry : database()) {
        Analysis analysis;
        try {
            analysis = analyze(moleculeFromEntry(entry));
        } catch (const std::exception& err) {
            problems.fail(entry.name + ": разбор завершился ошибкой — " + err.what());
            continue;
        }

        // в сводке сначала самые многочисленные гибридизации
        std::vector<std::pair<std::string, int>> summary = analysis.hybridSummary;
        std::stable_sort(summary.begin(), summary.end(),
                         [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
                             return a.second > b.second;
                         });
        std::string hyb;
        for (std::size_t k = 0; k < summary.size(); k++) {
            if (k > 0) hyb += " ";
            hyb += summary[k].first + "×" + std::to_string(summary[k].second);
        }

        say(pad(entry.name, 24) + pad(analysis.formulaHtml, 12) + pad(hyb.empty() ? "—" : hyb, 22)
            + pad(std::to_string(analysis.rings.size()), 7) + fixed(analysis.dipoleValue, 2));

        // --- геометрические проверки ---
        const Molecule& mol = analysis.molecule;
        for (const Bond& b : mol.bonds) {
            const double expected = bondLength(mol.atoms[b.a].el, mol.atoms[b.b].el, b.order);
            const double actual = dist(mol.atoms[b.a].pos, mol.atoms[b.b].pos);
            if (std::fabs(actual - expected) > 0.12) {
                problems.fail(entry.name + ": связь " + mol.atoms[b.a].el + "-" + mol.atoms[b.b].el
                              + " = " + fixed(actual, 2) + " Å, ожидалось " + fixed(expected, 2) + " Å");
            }
        }
        for (std::size_t i = 0; i < mol.atoms.size(); i++) {
            for (std::size_t j = i + 1; j < mol.atoms.size(); j++) {
                const double d = dist(mol.atoms[i].pos, mol.atoms[j].pos);
                if (d < 0.85) {
                    problems.fail(entry.name + ": атомы " + mol.atoms[i].el + std::to_string(i)
                                  + " и " + mol.atoms[j].el + std::to_string(j)
                                  + " слиплись (" + fixed(d, 2) + " Å)");
                }
            }
        }
        for (const Atom& a : mol.atoms) {
            if (!std::isfinite(a.pos.x + a.pos.y + a.pos.z)) {
                problems.fail(entry.name + ": у атома " + a.el + std::to_string(a.id)
                              + " нечисловые координаты");
            }
        }

        results.emplace_back(entry.name, std::move(analysis));
    }
    return results;
}

// ---------------------------------------------------------------------------
// 2. Сравнение с экспериментом
// ---------------------------------------------------------------------------

void compareWithExperiment(const std::vector<std::pair<std::string, Analysis>>& results) {
    say("");
    say("");
    say("=== ПРЕДСКАЗАНИЕ ПРОТИВ ЭКСПЕРИМЕНТА ===");
    say("");
    say(pad("Молекула", 24) + pad("Угол", 11) + padLeft("ОЭПВО", 8) + padLeft("модель", 8)
        + padLeft("опыт", 8) + padLeft("Δ", 7) + "  примечание");
    say(repeat("-", 84));

    int compared = 0;
    double sumAbsError = 0;
    int ringCases = 0;
    int heavyCases = 0;
    std::string worstName, worstLabel;
    double worstDelta = 0;

    for (const auto& item : results) {
        const std::string& name = item.first;
        const Analysis& analysis = item.second;
        for (const AngleRecord& angle : analysis.angles) {
            if (!angle.hasExperimental || !angle.hasPredicted) continue;

            const double delta = angle.predicted - angle.experimental;
            const AtomAnalysis& centerAtom = analysis.atoms[angle.center];
            const bool heavy = !angle.inRing && !centerAtom.warning.empty();

            if (angle.inRing) {
                ringCases++;
            } else if (heavy) {
                heavyCases++;
            } else {
                compared++;
                sumAbsError += std::fabs(delta);
                if (std::fabs(delta) > std::fabs(worstDelta)) {
                    worstName = name; worstLabel = angle.label; worstDelta = delta;
                }
            }

            std::string note;
            if (angle.inRing) note = "угол задан циклом";
            else if (heavy) note = "известная граница модели";
            else if (std::fabs(delta) > 6) note = "← расхождение";

            say(pad(name, 24) + pad(angle.label, 11)
                + padLeft(fixed(angle.predicted, 1), 8)
                + padLeft(fixed(angle.actual, 1), 8)
                + padLeft(fixed(angle.experimental, 1), 8)
                + padLeft((delta >= 0 ? "+" : "") + fixed(delta, 1), 7) + "  " + note);
        }
    }

    say(repeat("-", 84));
    say("Сравнений в области применимости модели: " + std::to_string(compared)
        + ", средняя ошибка " + fixed(sumAbsError / compared, 2) + "°");
    say("Наибольшее расхождение: " + worstName + ", " + worstLabel + ", " + fixed(worstDelta, 1) + "°");
    say("Отдельно учтены: углы в циклах — " + std::to_string(ringCases)
        + ", известные границы модели — " + std::to_string(heavyCases));
}

// ---------------------------------------------------------------------------
// 3. Контрольные молекулы: форма и гибридизация
// ---------------------------------------------------------------------------

struct Expectation {
    const char* id;
    /** Индекс центрального атома (в порядке разбора SMILES). */
    int center;
    const char* hybrid;
    const char* molecularGeom;
    bool hasAngle;
    double angleLo;
    double angleHi;
    int rings;       // -1 — не проверять
    int planar;      // -1 не проверять, 0 неплоское, 1 плоское
};

const Expectation EXPECTATIONS[] = {
    {"water", 0, "sp³", "угловая", true, 103, 106, -1, -1},
    {"ammonia", 0, "sp³", "тригонально-пирамидальная", true, 105, 108, -1, -1},
    {"methane", 0, "sp³", "тетраэдрическая", true, 109, 110, -1, -1},
    {"carbon-dioxide", 1, "sp", "линейная", true, 179, 181, -1, -1},
    {"boron-trifluoride", 1, "sp²", "плоская треугольная", true, 119, 121, -1, -1},
    {"ethene", 0, "sp²", "плоская треугольная", false, 0, 0, -1, -1},
    {"ethyne", 0, "sp", "линейная", true, 179, 181, -1, -1},
    {"allene", 1, "sp", "линейная", true, 179, 181, -1, -1},
    {"formaldehyde", 0, "sp²", "плоская треугольная", true, 115, 119, -1, -1},
    {"sulfur-dioxide", 1, "sp²", "угловая", true, 115, 120, -1, -1},
    {"phosphorus-pentachloride", 1, "sp³d", "тригонально-бипирамидальная", true, 89, 91, -1, -1},
    {"sulfur-tetrafluoride", 1, "sp³d", "качели (дисфеноид)", false, 0, 0, -1, -1},
    {"chlorine-trifluoride", 1, "sp³d", "Т-образная", true, 89, 91, -1, -1},
    {"xenon-difluoride", 1, "sp³d", "линейная", true, 179, 181, -1, -1},
    {"xenon-tetrafluoride", 1, "sp³d²", "квадратная", true, 89, 91, -1, -1},
    {"sulfur-hexafluoride", 1, "sp³d²", "октаэдрическая", true, 89, 91, -1, -1},
    {"iodine-pentafluoride", 1, "sp³d²", "квадратно-пирамидальная", false, 0, 0, -1, -1},
    {"benzene", 0, "sp²", "плоская треугольная", false, 0, 0, 1, 1},
    {"naphthalene", 0, "sp²", "плоская треугольная", false, 0, 0, 2, 1},
    {"cyclohexane", 0, "sp³", "тетраэдрическая", false, 0, 0, 1, 0},
    {"cyclopropane", 0, "sp³", "тетраэдрическая", false, 0, 0, 1, -1},
    {"pyridine", 3, "sp²", "угловая", false, 0, 0, 1, 1},
    {"nitrate", 1, "sp²", "плоская треугольная", true, 119, 121, -1, -1},
    {"sulfate", 1, "sp³", "тетраэдрическая", true, 109, 110, -1, -1},
    {"ammonium", 0, "sp³", "тетраэдрическая", true, 109, 110, -1, -1},
    {"dimethyl-sulfoxide", 1, "sp³", "тригонально-пирамидальная", false, 0, 0, -1, -1},
    {"glucose", 2, "sp³", "тетраэдрическая", false, 0, 0, 1, -1},
};

void runExpectations() {
    say("");
    say("");
    say("=== КОНТРОЛЬНЫЕ ПРОВЕРКИ ===");
    say("");

    for (const Expectation& exp : EXPECTATIONS) {
        const DbEntry* entry = entryById(exp.id);
        if (entry == nullptr) {
            problems.fail(std::string("Контроль: молекула «") + exp.id + "» отсутствует в базе");
            continue;
        }

        const Analysis analysis = analyze(moleculeFromEntry(*entry));
        const AtomAnalysis& atom = analysis.atoms[exp.center];
        std::vector<std::string> notes;

        if (atom.hybrid != exp.hybrid) {
            problems.fail(entry->name + ": гибридизация " + atom.hybrid + ", ожидалась " + exp.hybrid);
            notes.push_back("гибридизация " + atom.hybrid + " ≠ " + exp.hybrid);
        }
        if (atom.molecularGeom != exp.molecularGeom) {
            problems.fail(entry->name + ": форма «" + atom.molecularGeom + "», ожидалась «"
                          + exp.molecularGeom + "»");
            notes.push_back("форма «" + atom.molecularGeom + "»");
        }
        if (exp.hasAngle) {
            double minAngle = std::numeric_limits<double>::infinity();
            for (const AngleRecord& a : analysis.angles) {
                if (a.center == exp.center) minAngle = std::min(minAngle, a.actual);
            }
            if (minAngle < exp.angleLo || minAngle > exp.angleHi) {
                problems.fail(entry->name + ": минимальный угол " + fixed(minAngle, 1)
                              + "°, ожидался в диапазоне " + fixed(exp.angleLo, 0) + "–"
                              + fixed(exp.angleHi, 0) + "°");
                notes.push_back("угол " + fixed(minAngle, 1) + "°");
            }
        }
        if (exp.rings >= 0 && static_cast<int>(analysis.rings.size()) != exp.rings) {
            problems.fail(entry->name + ": найдено циклов " + std::to_string(analysis.rings.size())
                          + ", ожидалось " + std::to_string(exp.rings));
            notes.push_back("циклов " + std::to_string(analysis.rings.size()));
        }
        if (exp.planar >= 0) {
            std::set<int> ringIds;
            for (const Ring& r : analysis.rings) {
                for (int a : r.atoms) ringIds.insert(a);
            }
            std::vector<Vec3> points;
            for (const Atom& a : analysis.molecule.atoms) {
                if (ringIds.count(a.id) != 0) points.push_back(a.pos);
            }
            const double spread = maxPlaneDeviation(points);
            const bool isPlanar = spread < 0.12;
            if (isPlanar != (exp.planar == 1)) {
                problems.fail(entry->name + ": кольцо " + (isPlanar ? "плоское" : "неплоское")
                              + " (разброс " + fixed(spread, 2) + " Å), ожидалось "
                              + (exp.planar == 1 ? "плоское" : "неплоское"));
                notes.push_back("разброс " + fixed(spread, 2) + " Å");
            }
        }

        std::string status = "ок";
        if (!notes.empty()) {
            status = "ОШИБКА: ";
            for (std::size_t k = 0; k < notes.size(); k++) {
                if (k > 0) status += ", ";
                status += notes[k];
            }
        }
        say(pad(entry->name, 26) + pad(atom.axe, 8) + pad(atom.hybrid, 7)
            + pad(atom.molecularGeom, 30) + status);
    }
}

// ---------------------------------------------------------------------------
// 4. Сверка с эталоном, снятым с версии на TypeScript
// ---------------------------------------------------------------------------

/** Округление до нужного числа знаков — так же, как это делал выгружающий скрипт. */
double roundTo(double x, int digits) {
    double factor = 1;
    for (int k = 0; k < digits; k++) factor *= 10;
    return std::round(x * factor) / factor;
}

void checkAgainstReference(const std::vector<std::pair<std::string, Analysis>>& results) {
    say("");
    say("");
    say("=== СВЕРКА С ЭТАЛОНОМ (data/reference.json) ===");
    say("");

    Json reference;
    try {
        reference = Json::parseFile(findDataFile("reference.json"));
    } catch (const JsonError& err) {
        problems.fail(std::string("Эталон не прочитан: ") + err.what());
        say("Эталон не прочитан — сверка пропущена.");
        return;
    }

    // разбор по id: порядок молекул в эталоне тот же, но полагаться на это незачем
    const std::vector<Json>& refMolecules = reference["molecules"].items();

    int checkedMolecules = 0, checkedAtoms = 0, checkedBonds = 0, checkedAngles = 0;
    int mismatched = 0;
    // молекулы, где расхождение свелось к перестановке равноценных атомов
    std::set<std::string> degenerate;

    for (const DbEntry& entry : database()) {
        const Json* ref = nullptr;
        for (const Json& r : refMolecules) {
            if (r["id"].str() == entry.id) { ref = &r; break; }
        }
        if (ref == nullptr) {
            problems.fail("Эталон: молекулы «" + entry.id + "» нет в reference.json");
            continue;
        }

        const Analysis* found = nullptr;
        for (const auto& item : results) {
            if (item.second.molecule.meta != nullptr && item.second.molecule.meta->id == entry.id) {
                found = &item.second;
                break;
            }
        }
        if (found == nullptr) continue;  // о неразобранной молекуле уже сообщено
        const Analysis& a = *found;

        const std::size_t before = problems.count();
        const std::string who = "Эталон, " + entry.name + ": ";

        // --- молекула целиком ---
        if (a.formula != (*ref)["formula"].str()) {
            problems.fail(who + "формула " + a.formula + ", в эталоне " + (*ref)["formula"].str());
        }
        if (std::fabs(roundTo(a.mass, 3) - (*ref)["mass"].number()) > 0.0015) {
            problems.fail(who + "масса " + fixed(a.mass, 3) + ", в эталоне "
                          + fixed((*ref)["mass"].number(), 3));
        }
        if (std::fabs(totalCharge(a.molecule) - (*ref)["charge"].number()) > 1e-6) {
            problems.fail(who + "суммарный заряд разошёлся с эталоном");
        }
        if (a.polar != (*ref)["polar"].boolean()) {
            problems.fail(who + std::string("полярность ") + (a.polar ? "да" : "нет")
                          + ", в эталоне " + ((*ref)["polar"].boolean() ? "да" : "нет"));
        }

        // --- циклы (отсортированы по размеру, затем по ароматичности) ---
        std::vector<std::pair<int, bool>> myRings;
        for (const Ring& r : a.rings) myRings.emplace_back(static_cast<int>(r.atoms.size()), r.aromatic);
        std::stable_sort(myRings.begin(), myRings.end(),
                         [](const std::pair<int, bool>& x, const std::pair<int, bool>& y) {
                             if (x.first != y.first) return x.first < y.first;
                             return static_cast<int>(x.second) < static_cast<int>(y.second);
                         });
        const std::vector<Json>& refRings = (*ref)["rings"].items();
        if (myRings.size() != refRings.size()) {
            problems.fail(who + "циклов " + std::to_string(myRings.size()) + ", в эталоне "
                          + std::to_string(refRings.size()));
        } else {
            for (std::size_t k = 0; k < myRings.size(); k++) {
                if (myRings[k].first != static_cast<int>(refRings[k]["size"].number())
                    || myRings[k].second != refRings[k]["aromatic"].boolean()) {
                    problems.fail(who + "цикл " + std::to_string(k) + " не совпал с эталоном");
                }
            }
        }

        // --- атомы ---
        const std::vector<Json>& refAtoms = (*ref)["atoms"].items();
        if (a.atoms.size() != refAtoms.size()) {
            problems.fail(who + "атомов " + std::to_string(a.atoms.size()) + ", в эталоне "
                          + std::to_string(refAtoms.size()));
        } else {
            for (std::size_t k = 0; k < a.atoms.size(); k++) {
                const AtomAnalysis& mine = a.atoms[k];
                const Json& r = refAtoms[k];
                const std::string tag = who + "атом " + std::to_string(k) + ": ";
                if (mine.el != r["el"].str()) problems.fail(tag + "элемент " + mine.el);
                if (std::fabs(roundTo(mine.charge, 4) - r["q"].number()) > 1e-9)
                    problems.fail(tag + "заряд " + fixed(mine.charge, 4));
                if (mine.sigma != static_cast<int>(r["sigma"].number()))
                    problems.fail(tag + "σ-связей " + std::to_string(mine.sigma));
                if (mine.lonePairs != static_cast<int>(r["lp"].number()))
                    problems.fail(tag + "НЭП " + std::to_string(mine.lonePairs));
                if (mine.steric != static_cast<int>(r["sn"].number()))
                    problems.fail(tag + "стерическое число " + std::to_string(mine.steric));
                if (mine.axe != r["axe"].str())
                    problems.fail(tag + "тип " + mine.axe + ", в эталоне " + r["axe"].str());
                if (mine.hybrid != r["hybrid"].str())
                    problems.fail(tag + "гибридизация " + mine.hybrid + ", в эталоне " + r["hybrid"].str());
                if (mine.electronGeom != r["eGeom"].str())
                    problems.fail(tag + "электронная геометрия «" + mine.electronGeom + "»");
                if (mine.molecularGeom != r["mGeom"].str())
                    problems.fail(tag + "геометрия молекулы «" + mine.molecularGeom + "»");
                checkedAtoms++;
            }
        }

        // --- связи ---
        const std::vector<Json>& refBonds = (*ref)["bonds"].items();
        const Molecule& mol = a.molecule;
        bool lengthsInOrder = true;
        std::vector<double> myLengths, refLengths;
        if (mol.bonds.size() != refBonds.size()) {
            problems.fail(who + "связей " + std::to_string(mol.bonds.size()) + ", в эталоне "
                          + std::to_string(refBonds.size()));
        } else {
            for (std::size_t k = 0; k < mol.bonds.size(); k++) {
                const Bond& mine = mol.bonds[k];
                const Json& r = refBonds[k];
                const std::string tag = who + "связь " + std::to_string(k) + ": ";
                if (mine.a != static_cast<int>(r["a"].number()) || mine.b != static_cast<int>(r["b"].number()))
                    problems.fail(tag + "другие концы");
                if (std::fabs(roundTo(mine.order, 4) - r["order"].number()) > 1e-9)
                    problems.fail(tag + "кратность " + fixed(mine.order, 4));
                if (mine.aromatic != r["aromatic"].boolean())
                    problems.fail(tag + "признак ароматичности");
                if (mine.delocalized != r["delocalized"].boolean())
                    problems.fail(tag + "признак делокализации");
                const double length = dist(mol.atoms[mine.a].pos, mol.atoms[mine.b].pos);
                myLengths.push_back(length);
                refLengths.push_back(r["length"].number());
                if (std::fabs(length - r["length"].number()) > 0.01) lengthsInOrder = false;
                checkedBonds++;
            }
        }

        // --- валентные углы (в эталоне отсортированы по центру, затем по концам) ---
        std::vector<const AngleRecord*> myAngles;
        for (const AngleRecord& x : a.angles) myAngles.push_back(&x);
        std::stable_sort(myAngles.begin(), myAngles.end(),
                         [](const AngleRecord* x, const AngleRecord* y) {
                             if (x->center != y->center) return x->center < y->center;
                             if (x->i != y->i) return x->i < y->i;
                             return x->j < y->j;
                         });
        const std::vector<Json>& refAngles = (*ref)["angles"].items();
        bool anglesInOrder = true;
        if (myAngles.size() != refAngles.size()) {
            problems.fail(who + "углов " + std::to_string(myAngles.size()) + ", в эталоне "
                          + std::to_string(refAngles.size()));
        } else {
            for (std::size_t k = 0; k < myAngles.size(); k++) {
                const AngleRecord& mine = *myAngles[k];
                const Json& r = refAngles[k];
                if (mine.label != r["label"].str()
                    || std::fabs(mine.predicted - r["predicted"].number()) > 0.05
                    || std::fabs(mine.actual - r["actual"].number()) > 0.05) {
                    anglesInOrder = false;
                }
                checkedAngles++;
            }
        }

        // Если что-то разошлось поатомно — проверяем вырожденный минимум.
        // У метильной группы три равноценных водорода, и выбор конформации
        // из 36 вариантов при совпадающей энергии решается последним битом:
        // структура получается та же самая, но подписи двух водородов меняются
        // местами. Признак — совпадение НАБОРА значений при том же центре.
        // Лечение описано в разделе 7.1 задания.
        if (!lengthsInOrder) {
            std::vector<double> mineSorted = myLengths, refSorted = refLengths;
            std::sort(mineSorted.begin(), mineSorted.end());
            std::sort(refSorted.begin(), refSorted.end());
            bool sameSet = true;
            for (std::size_t k = 0; k < mineSorted.size(); k++) {
                if (std::fabs(mineSorted[k] - refSorted[k]) > 0.01) sameSet = false;
            }
            if (sameSet) {
                degenerate.insert(entry.name);
            } else {
                for (std::size_t k = 0; k < myLengths.size(); k++) {
                    if (std::fabs(myLengths[k] - refLengths[k]) > 0.01) {
                        problems.fail(who + "связь " + std::to_string(k) + ": длина "
                                      + fixed(myLengths[k], 3) + " Å, в эталоне "
                                      + fixed(refLengths[k], 3) + " Å");
                    }
                }
            }
        }

        if (!anglesInOrder) {
            // сравниваем отсортированные наборы углов при каждом центре
            std::map<int, std::vector<double>> minePredicted, mineActual;
            std::map<int, std::vector<double>> refPredicted, refActual;
            std::map<int, std::vector<std::string>> mineLabels, refLabels;
            for (const AngleRecord* x : myAngles) {
                minePredicted[x->center].push_back(x->predicted);
                mineActual[x->center].push_back(x->actual);
                mineLabels[x->center].push_back(x->label);
            }
            for (const Json& r : refAngles) {
                const int c = static_cast<int>(r["c"].number());
                refPredicted[c].push_back(r["predicted"].number());
                refActual[c].push_back(r["actual"].number());
                refLabels[c].push_back(r["label"].str());
            }

            bool sameSets = true;
            for (auto& item : minePredicted) {
                const int c = item.first;
                const auto compare = [&](std::vector<double> x, std::vector<double> y) {
                    if (x.size() != y.size()) return false;
                    std::sort(x.begin(), x.end());
                    std::sort(y.begin(), y.end());
                    for (std::size_t k = 0; k < x.size(); k++) {
                        if (std::fabs(x[k] - y[k]) > 0.05) return false;
                    }
                    return true;
                };
                std::vector<std::string> ml = mineLabels[c], rl = refLabels[c];
                std::sort(ml.begin(), ml.end());
                std::sort(rl.begin(), rl.end());
                if (ml != rl || !compare(item.second, refPredicted[c])
                    || !compare(mineActual[c], refActual[c])) {
                    sameSets = false;
                    problems.fail(who + "набор углов при центре " + std::to_string(c)
                                  + " не совпал с эталоном");
                }
            }
            if (sameSets) degenerate.insert(entry.name);
        }

        checkedMolecules++;
        if (problems.count() != before) mismatched++;
    }

    say("Сверено молекул: " + std::to_string(checkedMolecules)
        + ", атомов: " + std::to_string(checkedAtoms)
        + ", связей: " + std::to_string(checkedBonds)
        + ", валентных углов: " + std::to_string(checkedAngles) + ".");
    if (!degenerate.empty()) {
        std::string names;
        for (const std::string& n : degenerate) {
            if (!names.empty()) names += ", ";
            names += n;
        }
        say("Вырожденный минимум (набор значений совпал, различаются только подписи");
        say("равноценных атомов — см. раздел 7.1 задания): " + names + ".");
    }
    if (mismatched == 0) {
        say("Расхождений с эталоном нет.");
    } else {
        say("Молекул с расхождениями: " + std::to_string(mismatched) + ".");
    }
}

}  // namespace

int main() {
    testutil::enableUtf8Console();

    try {
        const std::vector<std::pair<std::string, Analysis>> results = runDatabase();
        compareWithExperiment(results);
        runExpectations();
        checkAgainstReference(results);
    } catch (const std::exception& err) {
        std::printf("\nСбой проверки: %s\n", err.what());
        return 1;
    }

    say("");
    say(repeat("=", 78));
    if (problems.count() == 0) {
        say("ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ. Молекул в базе: " + std::to_string(database().size()) + ".");
        return 0;
    }
    say("ОБНАРУЖЕНО ПРОБЛЕМ: " + std::to_string(problems.count()));
    say("");
    std::size_t shown = 0;
    for (const std::string& p : problems.all()) {
        say("  • " + p);
        if (++shown >= 60) {
            say("  … и ещё " + std::to_string(problems.count() - shown));
            break;
        }
    }
    return 1;
}
