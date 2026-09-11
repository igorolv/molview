/**
 * Проверка волновых функций и гибридных орбиталей (src/chem/orbital.cpp).
 * Запуск: build/orbitaltest.exe из корня проекта.
 *
 * Главное здесь — углы и ортогональность. Направления гибридов программа
 * берёт у теории Гиллеспи и ниоткуда больше; углы между ними считаются по
 * ПОСЧИТАННЫМ функциям, и то, что получается 109,47° / 120° / 180°, а разные
 * гибриды при этом ортогональны, — не подгонка, а следствие математики.
 * Поставь неверные коэффициенты при s и p, и ортогональность рассыплется.
 */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chem/orbital.h"
#include "chem/vec.h"
#include "testutil.h"

using namespace chem;
using testutil::fixed;
using testutil::pad;
using testutil::repeat;

namespace {

testutil::Problems problems;

void say(const std::string& line) { std::printf("%s\n", line.c_str()); }

bool close(double got, double want, double tol) { return std::abs(got - want) <= tol; }

void check(bool ok, const std::string& what, const std::string& detail) {
    say(pad(what, 42) + pad(detail, 32) + (ok ? "да" : "НЕТ"));
    if (!ok) problems.fail(what + ": " + detail);
}

/**
 * Ящик для проверок: ±8 Å с шагом 0,25 Å.
 *
 * Меньше нельзя. Водородоподобные функции при Z = 1 очень размазаны: внутри
 * сферы радиусом 4 Å лежит лишь 77 % плотности 2s, и на таком ящике нормировка
 * дала бы 0,9 вместо единицы — проверка перестала бы проверять что-либо, кроме
 * размера ящика. Сетке показа этого не требуется: там уровень изоповерхности
 * берётся в долях от наибольшего |ψ|, а не в абсолютных единицах.
 */
constexpr int BOX = 65;
constexpr double HALF = 8.0;

/** Сетка проверок для выбранной радиальной части. */
OrbitalOptions box(Radial radial, double zeff = 1.0) {
    OrbitalOptions o;
    o.size = BOX;
    // Ящик сжимается вместе с орбиталью, иначе для сжатой функции сетка
    // окажется слишком грубой.
    o.half = HALF / zeff;
    o.zeff = zeff;
    o.radial = radial;
    return o;
}

std::string nameOf(Radial radial) {
    return radial == Radial::Slater ? "Слейтер" : "водородная";
}

/** Тетраэдрический угол: arccos(−1/3). */
const double TETRAHEDRAL = std::acos(-1.0 / 3.0) * 180.0 / PI;

/** Направления гибридов идеальной электронной геометрии. */
std::vector<Vec3> tetrahedral() {
    const double t = 1.0 / std::sqrt(3.0);
    return {v3(t, t, t), v3(t, -t, -t), v3(-t, t, -t), v3(-t, -t, t)};
}

std::vector<Vec3> trigonal() {
    const double h = std::sqrt(3.0) / 2;
    return {v3(1, 0, 0), v3(-0.5, h, 0), v3(-0.5, -h, 0)};
}

std::vector<Vec3> linear() { return {v3(0, 0, 1), v3(0, 0, -1)}; }

/** Наибольшее отклонение угла между осями от идеального. */
double worstAngle(const std::vector<OrbitalGrid>& set, double ideal) {
    double worst = 0;
    for (std::size_t i = 0; i < set.size(); i++) {
        for (std::size_t j = i + 1; j < set.size(); j++) {
            worst = std::max(worst, std::abs(angleDeg(set[i].axis(), set[j].axis()) - ideal));
        }
    }
    return worst;
}

/** Наибольшее скалярное произведение РАЗНЫХ орбиталей — должно быть нулём. */
double worstCross(const std::vector<OrbitalGrid>& set) {
    double worst = 0;
    for (std::size_t i = 0; i < set.size(); i++) {
        for (std::size_t j = i + 1; j < set.size(); j++) {
            worst = std::max(worst, std::abs(overlap(set[i], set[j])));
        }
    }
    return worst;
}

/** Наибольшее отклонение нормировки от единицы. */
double worstNorm(const std::vector<OrbitalGrid>& set) {
    double worst = 0;
    for (const OrbitalGrid& g : set) worst = std::max(worst, std::abs(overlap(g, g) - 1.0));
    return worst;
}

/** Наибольшее расхождение оси орбитали с направлением, которое ей задали. */
double worstAxis(const std::vector<OrbitalGrid>& set, const std::vector<Vec3>& dirs) {
    double worst = 0;
    for (std::size_t i = 0; i < set.size(); i++) {
        worst = std::max(worst, angleDeg(set[i].axis(), dirs[i]));
    }
    return worst;
}

void reportHybrid(const std::string& title, Orbital kind, const std::vector<Vec3>& dirs,
                  double ideal, Radial radial) {
    const std::vector<OrbitalGrid> set = orbitalGrids(kind, dirs, v3(), box(radial));
    const std::size_t pairs = set.size() * (set.size() - 1) / 2;

    // Ось считается через первый момент, а у куба он ТОЧНО кратен направлению
    // гибрида (см. OrbitalGrid::axis), поэтому допуск нужен только на
    // разрядность double, а не на грубость сетки.
    check(worstAngle(set, ideal) < 1e-6, title + ": углы между осями",
          "все " + std::to_string(pairs) + " = " + fixed(ideal, 2) + "°");
    check(worstAxis(set, dirs) < 1e-4, title + ": ось = направление ОЭПВО",
          "расхождение " + fixed(worstAxis(set, dirs), 6) + "°");
    check(worstNorm(set) < 2e-3, title + ": нормировка каждого",
          "1 ± " + fixed(worstNorm(set), 5));
    check(worstCross(set) < 2e-3, title + ": разные ортогональны",
          "|⟨i|j⟩| ≤ " + fixed(worstCross(set), 5));
}

}  // namespace

int main() {
    testutil::enableUtf8Console();

    // Обе радиальные части проверяются одинаково: в задании записана
    // водородоподобная, а рисуется слейтеровская (см. chem::Radial), и
    // перечисленные свойства обязаны выполняться у обеих.
    for (Radial radial : {Radial::Hydrogen, Radial::Slater}) {
        const std::string tag = " [" + nameOf(radial) + "]";

        say("");
        say("=== 2s И 2p, радиальная часть: " + nameOf(radial) + " ===");
        say("");
        say(pad("что проверяется", 42) + pad("значение", 32) + "сошлось");
        say(repeat("-", 80));

        const OrbitalGrid s = orbitalGrid(Orbital::S, v3(0, 0, 1), v3(), box(radial));
        const OrbitalGrid pz = orbitalGrid(Orbital::P, v3(0, 0, 1), v3(), box(radial));
        const OrbitalGrid px = orbitalGrid(Orbital::P, v3(1, 0, 0), v3(), box(radial));

        check(close(overlap(s, s), 1, 2e-3), "2s: ∫ψ²dV = 1" + tag, fixed(overlap(s, s), 5));
        check(close(overlap(pz, pz), 1, 2e-3), "2p: ∫ψ²dV = 1" + tag, fixed(overlap(pz, pz), 5));
        check(std::abs(overlap(s, pz)) < 1e-9, "2s и 2p ортогональны" + tag,
              "⟨s|p⟩ = " + fixed(overlap(s, pz), 9));
        check(std::abs(overlap(px, pz)) < 1e-9, "2p_x и 2p_z ортогональны" + tag,
              "⟨p_x|p_z⟩ = " + fixed(overlap(px, pz), 9));

        // Узловая плоскость: при z = 0 у p_z РОВНО нули, а не «почти нули».
        // Средний слой сетки существует только потому, что число узлов по
        // ребру нечётное; на чётном проверять было бы нечего.
        const int middle = pz.size / 2;
        double node = 0;
        for (int j = 0; j < pz.size; j++) {
            for (int i = 0; i < pz.size; i++) node = std::max(node, std::abs(pz.at(i, j, middle)));
        }
        check(node == 0.0, "2p_z: узловая плоскость z = 0" + tag, "max|ψ| = " + fixed(node, 9));
        // А вне плоскости значения есть — иначе проверка выше проходила бы и
        // на тождественном нуле.
        check(pz.peak() > 0.05, "2p_z: вне плоскости не ноль" + tag,
              "peak = " + fixed(pz.peak(), 4));

        say("");
        say("=== ГИБРИДЫ: УГЛЫ И ОРТОНОРМИРОВАННОСТЬ, " + nameOf(radial) + " ===");
        say("");
        say("Направления взяты те, что даёт теория Гиллеспи для линейной,");
        say("тригональной и тетраэдрической электронной геометрии.");
        say("");
        say(pad("что проверяется", 42) + pad("значение", 32) + "сошлось");
        say(repeat("-", 80));

        reportHybrid("sp³" + tag, Orbital::Sp3, tetrahedral(), TETRAHEDRAL, radial);
        say("");
        reportHybrid("sp²" + tag, Orbital::Sp2, trigonal(), 120.0, radial);
        say("");
        reportHybrid("sp" + tag, Orbital::Sp, linear(), 180.0, radial);
    }

    say("");
    say("=== НЕЗАВИСИМОСТЬ ОТ ОСЕЙ СЕТКИ И ОТ ЗАРЯДА ЯДРА ===");
    say("");
    say(pad("что проверяется", 42) + pad("значение", 32) + "сошлось");
    say(repeat("-", 80));

    // Ни углы, ни ортогональность не должны зависеть от того, как гибриды
    // повёрнуты относительно куба. Повернём тетраэдр так, чтобы ни одно
    // направление не легло вдоль оси сетки.
    std::vector<Vec3> tilted;
    const Mat3 turn = rotationAxis(norm(v3(0.31, 1.0, -0.73)), 0.7);
    for (const Vec3& d : tetrahedral()) tilted.push_back(applyMat(turn, d));
    const std::vector<OrbitalGrid> turned =
        orbitalGrids(Orbital::Sp3, tilted, v3(), box(Radial::Hydrogen));
    check(worstAngle(turned, TETRAHEDRAL) < 1e-6, "sp³ повёрнутый: углы те же",
          fixed(TETRAHEDRAL, 2) + "° ± " + fixed(worstAngle(turned, TETRAHEDRAL), 6));
    check(worstAxis(turned, tilted) < 1e-4, "sp³ повёрнутый: оси на месте",
          "расхождение " + fixed(worstAxis(turned, tilted), 6) + "°");
    check(worstCross(turned) < 2e-3, "sp³ повёрнутый: ортогональны",
          "|⟨i|j⟩| ≤ " + fixed(worstCross(turned), 5));

    // Эффективный заряд ядра — параметр показа, он лишь сжимает орбиталь по
    // длине; на свойствах сказываться не должен.
    const double zeff = 3.2;
    const std::vector<OrbitalGrid> small =
        orbitalGrids(Orbital::Sp3, tetrahedral(), v3(), box(Radial::Slater, zeff));
    check(worstAngle(small, TETRAHEDRAL) < 1e-6, "sp³ при Z = 3,2: углы те же",
          fixed(TETRAHEDRAL, 2) + "° ± " + fixed(worstAngle(small, TETRAHEDRAL), 6));
    check(worstNorm(small) < 2e-3, "sp³ при Z = 3,2: нормировка",
          "1 ± " + fixed(worstNorm(small), 5));
    check(worstCross(small) < 2e-3, "sp³ при Z = 3,2: ортогональны",
          "|⟨i|j⟩| ≤ " + fixed(worstCross(small), 5));

    say("");
    say("=== ЧЕМ ОТЛИЧАЮТСЯ ДВЕ РАДИАЛЬНЫЕ ЧАСТИ ===");
    say("");
    say("Свойства у них одинаковы, а форма — нет, и от формы зависит,");
    say("можно ли на орбиталь смотреть. Сравниваются доли объёма сетки,");
    say("которые занимают положительная и отрицательная части гибрида sp³.");
    say("");
    say(pad("радиальная часть 2s", 42) + pad("ψ < 0 против ψ > 0", 32) + "сошлось");
    say(repeat("-", 80));

    // У водородоподобной 2s есть радиальный узел, за ним она меняет знак, и
    // отрицательная часть гибрида расплывается широким конусом. У
    // слейтеровской узла нет, и обратный лепесток выходит меньше переднего —
    // такой, каким его рисуют в учебнике.
    for (Radial radial : {Radial::Hydrogen, Radial::Slater}) {
        const OrbitalGrid g = orbitalGrid(Orbital::Sp3, v3(0, 0, 1), v3(), box(radial));
        const double iso = 0.22 * g.peak();
        std::size_t plus = 0, minus = 0;
        for (double v : g.values) {
            if (v >= iso) plus++;
            if (v <= -iso) minus++;
        }
        const double ratio = plus > 0 ? double(minus) / double(plus) : 0.0;
        const bool compact = radial == Radial::Slater ? ratio < 1.0 : ratio > 1.0;
        check(compact, "sp³, " + nameOf(radial),
              "отношение объёмов " + fixed(ratio, 2));
    }

    say("");
    say("=== ГИБРИД ПО СТЕРИЧЕСКОМУ ЧИСЛУ ===");
    say("");
    say(pad("что проверяется", 42) + pad("значение", 32) + "сошлось");
    say(repeat("-", 80));
    check(orbitalForSteric(2) == Orbital::Sp, "стерическое число 2 → sp", "линейная");
    check(orbitalForSteric(3) == Orbital::Sp2, "стерическое число 3 → sp²", "тригональная");
    check(orbitalForSteric(4) == Orbital::Sp3, "стерическое число 4 → sp³", "тетраэдрическая");
    check(orbitalForSteric(6) == Orbital::Sp3, "стерическое число 6 → форма sp³",
          "d-орбиталей в задании нет");

    say("");
    say(repeat("=", 80));
    if (problems.count() == 0) {
        say("ВСЕ ПРОВЕРКИ ПРОЙДЕНЫ. Тетраэдрический угол не подставлен, а получен:");
        say("гибриды с любым другим углом между осями не были бы ортогональны.");
        return 0;
    }
    say("ОБНАРУЖЕНО ПРОБЛЕМ: " + std::to_string(problems.count()));
    say("");
    for (const std::string& p : problems.all()) say("  • " + p);
    return 1;
}
