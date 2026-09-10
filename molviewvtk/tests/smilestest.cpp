/**
 * Проверка разбора и построения для молекул, которых НЕТ в базе.
 * Запуск: build/smilestest.exe из корня проекта.
 *
 * Именно эта проверка показывает, что программа — не листалка готовых картинок:
 * структура строится по строке SMILES для любой молекулы.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "chem/analyze.h"
#include "chem/periodic.h"
#include "chem/smiles.h"
#include "chem/types.h"
#include "chem/vec.h"
#include "testutil.h"

using namespace chem;
using testutil::fixed;
using testutil::pad;
using testutil::repeat;

namespace {

testutil::Problems problems;

void say(const std::string& line) { std::printf("%s\n", line.c_str()); }

struct Case {
    const char* smiles;
    const char* name;
    int rings;              // -1 — не проверять
    const char* formula;    // nullptr — не проверять
};

const Case CASES[] = {
    {"CC(=O)Oc1ccccc1C(=O)O", "Аспирин", 1, "C9H8O4"},
    {"CC(=O)Nc1ccc(O)cc1", "Парацетамол", 1, "C8H9NO2"},
    {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Кофеин", 2, "C8H10N4O2"},
    {"c1ccc2cc3ccccc3cc2c1", "Антрацен", 3, "C14H10"},
    {"c1ccc2[nH]ccc2c1", "Индол", 2, nullptr},
    {"c1cncnc1", "Пиримидин", 1, "C4H4N2"},
    {"C1C2CC3CC1CC(C2)C3", "Адамантан (мостиковый)", -1, "C10H16"},
    {"C1CC2CCC1C2", "Норборнан (мостиковый)", -1, nullptr},
    {"OCC(O)CO", "Глицерин", 0, "C3H8O3"},
    {"CC(N)C(=O)O", "Аланин", -1, "C3H7NO2"},
    {"CCOC(C)=O", "Этилацетат", -1, "C4H8O2"},
    {"OC(=O)CCC(=O)O", "Янтарная кислота", -1, nullptr},
    {"c1ccccc1[N+](=O)[O-]", "Нитробензол", 1, nullptr},
    {"Cc1ccc(cc1)S(=O)(=O)O", "п-Толуолсульфокислота", 1, nullptr},
    {"F[Br](F)(F)(F)F", "Пентафторид брома (sp³d²)", -1, nullptr},
    {"F[I](F)(F)(F)(F)(F)F", "Гептафторид иода (sp³d³)", -1, nullptr},
    {"[O-][Cl](=O)(=O)=O", "Перхлорат-ион", -1, nullptr},
    {"ClS(Cl)=O", "Хлористый тионил", -1, nullptr},
    {"ClP(Cl)(Cl)=O", "Хлорокись фосфора", -1, nullptr},
    {"[Na+].[Cl-]", "Разорванная запись (две частицы)", -1, nullptr},
    {"C1=CC=CC=C1", "Бензол в форме Кекуле", 1, "C6H6"},
    {"CC(C)(C)C(C)(C)C", "Гексаметилэтан (сильно загромождён)", -1, nullptr},
    {"OCC1OC(O)C(O)C(O)C1O", "Глюкоза", -1, nullptr},
    {"C/C=C/C", "Бутен-2 со стереометками", -1, nullptr},
    {"CC#CC", "Бутин-2", -1, nullptr},
};

struct BadCase {
    const char* smiles;
    const char* why;
};

const BadCase BAD[] = {
    {"C1CC", "незамкнутый цикл"},
    {"CC(", "незакрытая скобка"},
    {"CC)", "лишняя скобка"},
    {"CQC", "неизвестный символ"},
    {"[Zz]", "неизвестный элемент"},
    {"", "пустая строка"},
};

}  // namespace

int main() {
    testutil::enableUtf8Console();

    say("");
    say("=== ПОСТРОЕНИЕ ПО SMILES (молекул нет в базе) ===");
    say("");
    say(pad("Молекула", 32) + pad("Формула", 12) + pad("Атомов", 8) + pad("Циклов", 8)
        + pad("Мин. расст.", 13) + "Макс. ошибка связи");
    say(repeat("-", 96));

    for (const Case& test : CASES) {
        Analysis analysis;
        try {
            analysis = analyze(parseSmiles(test.smiles));
        } catch (const std::exception& err) {
            problems.fail(std::string(test.name) + ": " + err.what());
            say(pad(test.name, 32) + "ОШИБКА: " + err.what());
            continue;
        }

        const Molecule& mol = analysis.molecule;
        double minDist = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < mol.atoms.size(); i++) {
            for (std::size_t j = i + 1; j < mol.atoms.size(); j++) {
                minDist = std::min(minDist, dist(mol.atoms[i].pos, mol.atoms[j].pos));
            }
        }
        double maxBondError = 0;
        for (const Bond& b : mol.bonds) {
            const double expected = bondLength(mol.atoms[b.a].el, mol.atoms[b.b].el, b.order);
            maxBondError = std::max(maxBondError,
                std::fabs(dist(mol.atoms[b.a].pos, mol.atoms[b.b].pos) - expected));
        }

        say(pad(test.name, 32) + pad(analysis.formula, 12)
            + pad(std::to_string(mol.atoms.size()), 8)
            + pad(std::to_string(analysis.rings.size()), 8)
            + pad(fixed(minDist, 2) + " Å", 13) + fixed(maxBondError, 2) + " Å");

        if (minDist < 0.85) {
            problems.fail(std::string(test.name) + ": атомы слиплись (" + fixed(minDist, 2) + " Å)");
        }
        if (maxBondError > 0.25) {
            problems.fail(std::string(test.name) + ": длина связи отклонилась на "
                          + fixed(maxBondError, 2) + " Å");
        }
        for (const Atom& a : mol.atoms) {
            if (!std::isfinite(a.pos.x + a.pos.y + a.pos.z)) {
                problems.fail(std::string(test.name) + ": нечисловые координаты");
                break;
            }
        }
        if (test.rings >= 0 && static_cast<int>(analysis.rings.size()) != test.rings) {
            problems.fail(std::string(test.name) + ": циклов "
                          + std::to_string(analysis.rings.size())
                          + ", ожидалось " + std::to_string(test.rings));
        }
        if (test.formula != nullptr && analysis.formula != test.formula) {
            problems.fail(std::string(test.name) + ": формула " + analysis.formula
                          + ", ожидалась " + test.formula);
        }
    }

    say("");
    say("=== ОБРАБОТКА ОШИБОК ВВОДА ===");
    say("");
    for (const BadCase& bad : BAD) {
        try {
            parseSmiles(bad.smiles);
            problems.fail(std::string("«") + bad.smiles + "» (" + bad.why + ") разобралось без ошибки");
            say(pad(std::string("«") + bad.smiles + "»", 16)
                + "ОШИБКА: разбор прошёл, хотя не должен был");
        } catch (const SmilesError& err) {
            say(pad(std::string("«") + bad.smiles + "»", 16) + pad(bad.why, 22) + "→ " + err.what());
        } catch (const std::exception& err) {
            problems.fail(std::string("«") + bad.smiles + "»: исключение не того типа");
            say(pad(std::string("«") + bad.smiles + "»", 16) + pad(bad.why, 22)
                + "→ (не SmilesError) " + err.what());
        }
    }

    say("");
    say(repeat("=", 96));
    if (problems.count() == 0) {
        say("ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ. Разобрано молекул вне базы: "
            + std::to_string(sizeof(CASES) / sizeof(CASES[0])) + ".");
        return 0;
    }
    say("ОБНАРУЖЕНО ПРОБЛЕМ: " + std::to_string(problems.count()));
    say("");
    for (const std::string& p : problems.all()) say("  • " + p);
    return 1;
}
