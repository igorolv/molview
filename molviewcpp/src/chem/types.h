#pragma once

// Структуры данных вычислительного ядра: граф молекулы и результат разбора.
// Никакой графики и никакого Windows — этот заголовок одинаково пригоден
// и для консольных тестов, и для оконной программы.

#include <string>
#include <utility>
#include <vector>

#include "chem/vec.h"

namespace chem {

/** Атом в графе молекулы. */
struct Atom {
    int id = 0;
    std::string el;
    /** Формальный заряд. При делокализации становится дробным. */
    double charge = 0.0;
    /** Входит ли атом в ароматическую систему. */
    bool aromatic = false;
    /** Число неявных атомов водорода (до их развёртывания в явные). */
    int hCount = 0;
    /** Координаты, Å. Заполняются построителем геометрии. */
    Vec3 pos;
    bool placed = false;
    /** true — атом появился при развёртывании неявных водородов. */
    bool fromImplicitH = false;
};

/** Связь. order: 1 / 2 / 3, для ароматических — 1,5, для делокализованных дробный. */
struct Bond {
    int a = 0;
    int b = 0;
    double order = 1.0;
    bool aromatic = false;
    /** Признак замыкающей цикл связи (для отладки построителя). */
    bool ringClosure = false;
    /** Кратность усреднена из-за делокализации (нитрат-ион, карбоксилат). */
    bool delocalized = false;
};

/** Запись справочника molecules.json. */
struct MoleculeMeta {
    std::string id;
    std::string name;
    std::string smiles;
    std::string cat;
    std::vector<std::string> syn;
    /**
     * Экспериментальные валентные углы, ключ вида «H-O-H».
     * Для одной подписи значений может быть несколько (в SF₆ это и 90°, и 180°).
     * Порядок хранения — как в файле: от него зависит привязка к углам.
     */
    std::vector<std::pair<std::string, std::vector<double>>> exp;
    /** Экспериментальный дипольный момент, Д. Отрицательное — «не указан». */
    double dipole = -1.0;
    bool hasDipole = false;
    std::string note;
};

struct Molecule {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    /** Заголовок из базы, если молекула найдена в справочнике. */
    const MoleculeMeta* meta = nullptr;
};

/** Результат применения теории Гиллеспи к одному атому. */
struct AtomAnalysis {
    int id = 0;
    std::string el;
    double charge = 0.0;
    /** Число σ-связей = число соседей. */
    int sigma = 0;
    /** Число неподелённых электронных пар. */
    int lonePairs = 0;
    /** Стерическое число = σ + НЭП. */
    int steric = 0;
    /** Обозначение типа: AX₄E₂. */
    std::string axe;
    std::string hybrid;
    /** Геометрия электронных пар. */
    std::string electronGeom;
    /** Геометрия молекулы (по положениям ядер). */
    std::string molecularGeom;
    /** Идеальный угол для электронной геометрии. */
    double idealAngle = 0.0;
    bool hasIdealAngle = false;
    /** Угол с учётом поправки на неподелённые пары. */
    double predictedAngle = 0.0;
    bool hasPredictedAngle = false;
    /** Является ли атом центральным (два и более соседа). */
    bool isCentral = false;
    /** Направления неподелённых пар (единичные векторы). */
    std::vector<Vec3> lonePairDirs;
    /** Направления гибридных орбиталей (все, включая связывающие). */
    std::vector<Vec3> orbitalDirs;
    /** Пояснение, если модель даёт нетипичный результат. Пусто — предупреждения нет. */
    std::string warning;
};

/** Один валентный угол. */
struct AngleRecord {
    int i = 0;
    int center = 0;
    int j = 0;
    std::string label;
    /** Угол, измеренный по построенным координатам. */
    double actual = 0.0;
    /** Предсказание ОЭПВО. */
    double predicted = 0.0;
    bool hasPredicted = false;
    /** Экспериментальное значение из справочника. */
    double experimental = 0.0;
    bool hasExperimental = false;
    /** Оба атома входят в один цикл — угол задан геометрией кольца. */
    bool inRing = false;
};

struct Ring {
    std::vector<int> atoms;
    bool aromatic = false;
    bool planar = false;
};

struct Analysis {
    Molecule molecule;
    std::string formula;
    std::string formulaHtml;
    double mass = 0.0;
    std::vector<AtomAnalysis> atoms;
    std::vector<AngleRecord> angles;
    std::vector<Ring> rings;
    /** Расчётный вектор дипольного момента (качественная оценка). */
    Vec3 dipoleVec;
    double dipoleValue = 0.0;
    bool polar = false;
    /** Сводка по гибридизациям в порядке появления, например sp³×2, sp²×1. */
    std::vector<std::pair<std::string, int>> hybridSummary;
};

}  // namespace chem
