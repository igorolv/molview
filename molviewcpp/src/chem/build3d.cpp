#include "chem/build3d.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <string>

#include "chem/periodic.h"
#include "chem/rings.h"
#include "chem/smiles.h"

namespace chem {

namespace {

std::string pairKey(int a, int b) {
    return a < b ? std::to_string(a) + "-" + std::to_string(b)
                 : std::to_string(b) + "-" + std::to_string(a);
}

/** Корневой атом: самый «нагруженный» тяжёлый атом — от него удобнее расти. */
int pickRoot(const Molecule& mol, const std::vector<AtomEnvironment>& env) {
    int best = 0;
    int bestScore = -1;
    for (const Atom& atom : mol.atoms) {
        const int score = (atom.el == "H" ? 0 : 100) + static_cast<int>(env[atom.id].neighbors.size());
        if (score > bestScore) { bestScore = score; best = atom.id; }
    }
    return best;
}

double spanOf(const Molecule& mol) {
    double maxLen = 0;
    for (const Atom& a : mol.atoms) {
        if (a.placed) maxLen = std::max(maxLen, len(a.pos));
    }
    return maxLen;
}

/**
 * Подбор поворота вокруг единственной известной связи — выбор конформации.
 * Перебираем 36 положений и берём то, при котором новые атомы дальше всего
 * от уже построенных. Так сам собой получается заторможенный этан.
 */
Mat3 chooseTwist(const Molecule& mol, const std::vector<AtomEnvironment>& env,
                 int id, const Mat3& base, const Vec3& axis) {
    const AtomEnvironment& e = env[id];
    const Atom& me = mol.atoms[id];

    std::vector<std::pair<int, int>> pending;  // (nbId, slot)
    for (std::size_t slot = 0; slot < e.neighbors.size(); slot++) {
        const int nbId = e.neighbors[slot];
        if (!mol.atoms[nbId].placed) pending.emplace_back(nbId, static_cast<int>(slot));
    }
    if (pending.empty()) return base;

    std::vector<const Atom*> others;
    for (const Atom& a : mol.atoms) {
        if (a.placed && a.id != id) others.push_back(&a);
    }
    if (others.size() <= 1) return base;

    Mat3 bestRot = base;
    double bestScore = std::numeric_limits<double>::infinity();

    for (int step = 0; step < 36; step++) {
        const double angle = (step * PI * 2) / 36;
        const Mat3 rot = matMul(rotationAxis(axis, angle), base);

        double score = 0;
        for (const auto& item : pending) {
            const int nbId = item.first;
            const int slot = item.second;
            const Vec3 dir = norm(applyMat(rot, e.localDirs[slot]));
            const double L = bondLength(me.el, mol.atoms[nbId].el, e.bondOrders[slot]);
            const Vec3 p = add(me.pos, mul(dir, L));
            for (const Atom* other : others) {
                const double d = std::max(0.4, dist(p, other->pos));
                score += 1 / (d * d * d * d);
            }
        }
        if (score < bestScore) { bestScore = score; bestRot = rot; }
    }
    return bestRot;
}

Mat3 resolveFrame(const Molecule& mol, const std::vector<AtomEnvironment>& env, int id) {
    const AtomEnvironment& e = env[id];
    const Atom& me = mol.atoms[id];

    struct Placed { int slot; Vec3 local; Vec3 world; };
    std::vector<Placed> placed;

    for (std::size_t slot = 0; slot < e.neighbors.size(); slot++) {
        const int nbId = e.neighbors[slot];
        if (!mol.atoms[nbId].placed) continue;
        const Vec3 world = sub(mol.atoms[nbId].pos, me.pos);
        if (len(world) < 1e-6) continue;
        placed.push_back({static_cast<int>(slot), e.localDirs[slot], norm(world)});
    }

    if (placed.empty()) return IDENTITY;

    if (placed.size() == 1) {
        const Mat3 base = rotationFromTo(placed[0].local, placed[0].world);
        return chooseTwist(mol, env, id, base, placed[0].world);
    }

    // выбираем пару направлений с наибольшим углом между ними — так устойчивее
    std::size_t bi = 0, bj = 1;
    double bestSep = -2;
    for (std::size_t i = 0; i < placed.size(); i++) {
        for (std::size_t j = i + 1; j < placed.size(); j++) {
            const double sep = -dot(placed[i].world, placed[j].world);
            if (sep > bestSep) { bestSep = sep; bi = i; bj = j; }
        }
    }

    // Совмещаем по биссектрисе и нормали: это симметрично и не «перекашивает»
    // окружение, когда реальный угол в цикле отличается от идеального.
    Vec3 localSum = v3(), worldSum = v3();
    for (const Placed& p : placed) {
        localSum = add(localSum, p.local);
        worldSum = add(worldSum, p.world);
    }

    const Vec3 localNormal = cross(placed[bi].local, placed[bj].local);
    const Vec3 worldNormal = cross(placed[bi].world, placed[bj].world);

    if (len(localSum) > 0.15 && len(worldSum) > 0.15
        && len(localNormal) > 0.05 && len(worldNormal) > 0.05) {
        return alignPair(localSum, localNormal, worldSum, worldNormal);
    }
    return alignPair(placed[bi].local, placed[bj].local, placed[bi].world, placed[bj].world);
}

// ---------------------------------------------------------------------------
// Укладка циклов по шаблонам
// ---------------------------------------------------------------------------

/** Есть ли атом, принадлежащий трём и более циклам группы. */
bool hasBridgehead(const std::vector<Ring>& group) {
    std::map<int, int> count;
    for (const Ring& ring : group) {
        for (int a : ring.atoms) count[a]++;
    }
    for (const auto& item : count) {
        if (item.second >= 3) return true;
    }
    return false;
}

/** Средняя длина связи в кольце. */
double ringBondLength(const Molecule& mol, const Ring& ring) {
    double sum = 0;
    for (std::size_t i = 0; i < ring.atoms.size(); i++) {
        const int a = ring.atoms[i];
        const int b = ring.atoms[(i + 1) % ring.atoms.size()];
        const Bond* bond = findBond(mol, a, b);
        sum += bondLength(mol.atoms[a].el, mol.atoms[b].el, bond != nullptr ? bond->order : 1.0);
    }
    return sum / static_cast<double>(ring.atoms.size());
}

double polygonRadius(int sides, double side) {
    return side / (2 * std::sin(PI / sides));
}

bool tryAttachPlanarRing(Molecule& mol, const Ring& ring, std::set<int>& templatePlaced) {
    const int size = static_cast<int>(ring.atoms.size());
    std::vector<int> placedIdx;
    for (int k = 0; k < size; k++) {
        if (mol.atoms[ring.atoms[k]].placed) placedIdx.push_back(k);
    }
    if (placedIdx.size() < 2) return false;

    // ищем два соседних в кольце уже размещённых атома — общую связь
    int anchorA = -1, anchorB = -1;
    for (int k : placedIdx) {
        const int next = (k + 1) % size;
        if (mol.atoms[ring.atoms[next]].placed) { anchorA = k; anchorB = next; break; }
    }
    if (anchorA < 0) return false;

    const Vec3 pA = mol.atoms[ring.atoms[anchorA]].pos;
    const Vec3 pB = mol.atoms[ring.atoms[anchorB]].pos;
    const double L = ringBondLength(mol, ring);
    const double R = polygonRadius(size, L);
    const double apothem = std::sqrt(std::max(0.0, R * R - (L / 2) * (L / 2)));

    const Vec3 mid = mul(add(pA, pB), 0.5);
    const Vec3 along = norm(sub(pB, pA));
    // нормаль в плоскости z = 0
    Vec3 outward = norm(v3(-along.y, along.x, 0));

    // направление «наружу» от уже построенной части
    std::vector<Vec3> built;
    for (const Atom& a : mol.atoms) {
        if (a.placed) built.push_back(a.pos);
    }
    const Vec3 bulk = centroid(built);
    if (dot(outward, sub(mid, bulk)) < 0) outward = mul(outward, -1);

    const Vec3 center = add(mid, mul(outward, apothem));

    const auto angleOf = [&center](const Vec3& p) { return std::atan2(p.y - center.y, p.x - center.x); };
    const double aAngle = angleOf(pA);
    const double bAngle = angleOf(pB);
    const double stepSize = (2 * PI) / size;
    double delta = bAngle - aAngle;
    while (delta > PI) delta -= 2 * PI;
    while (delta < -PI) delta += 2 * PI;
    const double direction = delta > 0 ? 1 : -1;

    for (int step = 1; step < size; step++) {
        const int ringIndex = (anchorB + step) % size;
        const int id = ring.atoms[ringIndex];
        if (mol.atoms[id].placed) continue;
        const double a = bAngle + direction * stepSize * step;
        mol.atoms[id].pos = v3(center.x + R * std::cos(a), center.y + R * std::sin(a), 0);
        mol.atoms[id].placed = true;
        templatePlaced.insert(id);
    }
    return true;
}

/** Плоская (ароматическая или сопряжённая) система колец — всё в плоскости z = 0. */
void placePlanarSystem(Molecule& mol, const std::vector<Ring>& group, std::set<int>& templatePlaced) {
    const Ring& first = group[0];
    const double L = ringBondLength(mol, first);
    const double R = polygonRadius(static_cast<int>(first.atoms.size()), L);

    for (std::size_t k = 0; k < first.atoms.size(); k++) {
        const int id = first.atoms[k];
        const double a = (2 * PI * k) / first.atoms.size() + PI / 2;
        mol.atoms[id].pos = v3(R * std::cos(a), R * std::sin(a), 0);
        mol.atoms[id].placed = true;
        templatePlaced.insert(id);
    }

    std::vector<Ring> remaining(group.begin() + 1, group.end());
    bool progress = true;
    while (!remaining.empty() && progress) {
        progress = false;
        for (std::size_t idx = 0; idx < remaining.size(); idx++) {
            if (tryAttachPlanarRing(mol, remaining[idx], templatePlaced)) {
                remaining.erase(remaining.begin() + static_cast<long>(idx));
                progress = true;
                break;
            }
        }
    }
    // непристроенные кольца достроит общий обход
}

/** Насыщенный цикл: «кресло» для шестичленного, «конверт» для пятичленного. */
void placeSaturatedRing(Molecule& mol, const Ring& ring, std::set<int>& templatePlaced) {
    const int size = static_cast<int>(ring.atoms.size());
    const double L = ringBondLength(mol, ring);

    if (size == 6) {
        // Классическое кресло: радиус 0,948·L, отклонение по оси ±0,162·L.
        // При L = 1,54 Å это даёт длину связи 1,54 Å и валентный угол 111,5°.
        const double R = 0.948 * L;
        const double h = 0.162 * L;
        for (int k = 0; k < size; k++) {
            const int id = ring.atoms[k];
            const double a = (2 * PI * k) / 6;
            mol.atoms[id].pos = v3(R * std::cos(a), R * std::sin(a), k % 2 == 0 ? h : -h);
            mol.atoms[id].placed = true;
            templatePlaced.insert(id);
        }
    } else if (size == 5) {
        // Конверт: четыре атома в плоскости, пятый приподнят.
        const double R = polygonRadius(5, L);
        for (int k = 0; k < size; k++) {
            const int id = ring.atoms[k];
            const double a = (2 * PI * k) / 5;
            mol.atoms[id].pos = v3(R * std::cos(a), R * std::sin(a), k == 4 ? 0.28 * L : 0);
            mol.atoms[id].placed = true;
            templatePlaced.insert(id);
        }
    } else {
        const double R = polygonRadius(size, L);
        for (int k = 0; k < size; k++) {
            const int id = ring.atoms[k];
            const double a = (2 * PI * k) / size;
            // лёгкая складка для четырёхчленного цикла
            const double z = size == 4 ? (k % 2 == 0 ? 0.05 * L : -0.05 * L) : 0;
            mol.atoms[id].pos = v3(R * std::cos(a), R * std::sin(a), z);
            mol.atoms[id].placed = true;
            templatePlaced.insert(id);
        }
    }
}

void placeRings(Molecule& mol, const std::vector<AtomEnvironment>& env,
                const std::vector<Ring>& rings, std::set<int>& templatePlaced) {
    if (rings.empty()) return;

    for (const std::vector<Ring>& group : fuseRingSystems(rings)) {
        bool allPlanar = true;
        for (const Ring& r : group) {
            for (int a : r.atoms) {
                if (!(env[a].steric <= 3 || mol.atoms[a].aromatic)) { allPlanar = false; break; }
            }
            if (!allPlanar) break;
        }

        // Мостиковые системы (адамантан, норборнан): есть атом, входящий сразу
        // в три цикла — ни плоским многоугольником, ни «креслом» такое не уложить.
        // Их полностью строит обход графа, а расхождения правит релаксация.
        if (hasBridgehead(group)) continue;

        if (allPlanar) {
            placePlanarSystem(mol, group, templatePlaced);
        } else {
            // насыщенная система: аккуратно кладём только первое кольцо,
            // остальное достроит обход графа
            placeSaturatedRing(mol, group[0], templatePlaced);
        }
    }
}

// ---------------------------------------------------------------------------
// Релаксация: разводим атомы, которые оказались слишком близко
// ---------------------------------------------------------------------------

std::set<std::string> buildProximity(const Molecule& mol, int maxBonds) {
    const std::size_t n = mol.atoms.size();
    const std::vector<std::vector<int>> nb = neighborLists(mol);

    std::set<std::string> near;
    for (std::size_t start = 0; start < n; start++) {
        std::map<int, int> depth;
        depth.emplace(static_cast<int>(start), 0);
        std::deque<int> queue{static_cast<int>(start)};
        while (!queue.empty()) {
            const int cur = queue.front();
            queue.pop_front();
            const int d = depth[cur];
            if (d >= maxBonds) continue;
            for (int next : nb[cur]) {
                if (depth.count(next) != 0) continue;
                depth.emplace(next, d + 1);
                near.insert(pairKey(static_cast<int>(start), next));
                queue.push_back(next);
            }
        }
    }
    return near;
}

void relax(Molecule& mol, const std::vector<AtomEnvironment>& env, const std::set<int>& frozen) {
    const std::size_t n = mol.atoms.size();
    if (n < 4) return;

    std::map<std::string, double> targetBond;
    for (const Bond& b : mol.bonds) {
        targetBond[pairKey(b.a, b.b)] =
            bondLength(mol.atoms[b.a].el, mol.atoms[b.b].el, b.order);
    }

    // целевые расстояния 1–3 (через один атом) — они и держат валентные углы
    struct Target13 { int i; int j; double d; };
    std::vector<Target13> target13;
    for (std::size_t id = 0; id < n; id++) {
        const AtomEnvironment& e = env[id];
        for (std::size_t a = 0; a < e.neighbors.size(); a++) {
            for (std::size_t b = a + 1; b < e.neighbors.size(); b++) {
                const int i = e.neighbors[a], j = e.neighbors[b];
                // угол внутри уложенного по шаблону кольца задан геометрией самого кольца
                if (frozen.count(static_cast<int>(id)) != 0 && frozen.count(i) != 0
                    && frozen.count(j) != 0) continue;
                const double angle = std::acos(std::max(-1.0, std::min(1.0,
                    dot(e.localDirs[a], e.localDirs[b]))));
                const double la = bondLength(mol.atoms[id].el, mol.atoms[i].el, e.bondOrders[a]);
                const double lb = bondLength(mol.atoms[id].el, mol.atoms[j].el, e.bondOrders[b]);
                target13.push_back({i, j, std::sqrt(la * la + lb * lb - 2 * la * lb * std::cos(angle))});
            }
        }
    }

    const std::set<std::string> bondedWithin3 = buildProximity(mol, 3);
    std::vector<bool> movable(n);
    for (std::size_t i = 0; i < n; i++) movable[i] = frozen.count(static_cast<int>(i)) == 0;

    for (int iter = 0; iter < 600; iter++) {
        std::vector<Vec3> shift(n, v3());

        const auto applyPair = [&](int i, int j, double targetDist, double k) {
            const Vec3 d = sub(mol.atoms[j].pos, mol.atoms[i].pos);
            const double l = len(d);
            if (l < 1e-6) return;
            const Vec3 correction = mul(d, ((l - targetDist) / l) * k * 0.5);
            shift[i] = add(shift[i], correction);
            shift[j] = sub(shift[j], correction);
        };

        for (const Bond& b : mol.bonds) applyPair(b.a, b.b, targetBond[pairKey(b.a, b.b)], 0.85);
        for (const Target13& t : target13) applyPair(t.i, t.j, t.d, 0.35);

        // мягкое отталкивание несвязанных атомов
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = i + 1; j < n; j++) {
                if (bondedWithin3.count(pairKey(static_cast<int>(i), static_cast<int>(j))) != 0) continue;
                const double d = dist(mol.atoms[i].pos, mol.atoms[j].pos);
                const double minDist = 2.4;
                if (d < minDist && d > 1e-6) {
                    applyPair(static_cast<int>(i), static_cast<int>(j), minDist, 0.25);
                }
            }
        }

        for (std::size_t i = 0; i < n; i++) {
            if (!movable[i]) continue;
            mol.atoms[i].pos = add(mol.atoms[i].pos, mul(shift[i], 0.6));
        }
    }
}

// ---------------------------------------------------------------------------
// Разворот молекулы по главным осям
// ---------------------------------------------------------------------------

Mat3 orientByPrincipalAxes(Molecule& mol) {
    std::vector<int> heavy;
    for (const Atom& a : mol.atoms) {
        if (a.el != "H") heavy.push_back(a.id);
    }
    std::vector<int> reference;
    if (heavy.size() >= 2) {
        reference = heavy;
    } else {
        for (const Atom& a : mol.atoms) reference.push_back(a.id);
    }

    std::vector<Vec3> refPositions;
    for (int id : reference) refPositions.push_back(mol.atoms[id].pos);
    const Vec3 center = centroid(refPositions);

    for (Atom& a : mol.atoms) a.pos = sub(a.pos, center);

    if (mol.atoms.size() < 3) return IDENTITY;

    // матрица ковариации положений
    double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
    for (int id : reference) {
        const Vec3 p = mol.atoms[id].pos;
        xx += p.x * p.x; xy += p.x * p.y; xz += p.x * p.z;
        yy += p.y * p.y; yz += p.y * p.z; zz += p.z * p.z;
    }
    const Eigen3 eigen = jacobiEigen(Mat3{xx, xy, xz, xy, yy, yz, xz, yz, zz});

    std::vector<int> order{0, 1, 2};
    std::stable_sort(order.begin(), order.end(),
                     [&eigen](int a, int b) { return eigen.values[a] > eigen.values[b]; });

    const Vec3 ex = norm(eigen.vectors[order[0]]);
    Vec3 ey = norm(eigen.vectors[order[1]]);
    const Vec3 ez = cross(ex, ey);
    ey = cross(ez, ex);

    for (Atom& a : mol.atoms) {
        const Vec3 p = a.pos;
        a.pos = v3(dot(p, ex), dot(p, ey), dot(p, ez));
    }

    // небольшой наклон, чтобы объём читался лучше, чем строго анфас
    const Mat3 tilt = rotationAxis(v3(1, 0, 0), -0.32);
    const Mat3 spin = rotationAxis(v3(0, 1, 0), 0.42);
    for (Atom& a : mol.atoms) a.pos = applyMat(spin, applyMat(tilt, a.pos));

    const Mat3 axes = {ex.x, ex.y, ex.z, ey.x, ey.y, ey.z, ez.x, ez.y, ez.z};
    return matMul(spin, matMul(tilt, axes));
}

}  // namespace

BuildResult buildGeometry(Molecule& mol, const std::vector<AtomEnvironment>& env,
                          const std::vector<Ring>& rings) {
    const std::size_t n = mol.atoms.size();
    std::vector<Mat3> frames(n, IDENTITY);
    std::vector<bool> frameKnown(n, false);
    for (Atom& atom : mol.atoms) { atom.placed = false; atom.pos = v3(); }

    // Атомы, уложенные по геометрическому шаблону кольца. Только их имеет смысл
    // держать неподвижными при релаксации: их геометрия заведомо правильная.
    // Мостиковые и прочие сложные системы шаблоном не покрываются, их достраивает
    // обход графа — и такие атомы релаксация обязана поправлять.
    std::set<int> templatePlaced;
    placeRings(mol, env, rings, templatePlaced);

    // --- наращивание остальной молекулы ---------------------------------------
    std::deque<int> queue;
    for (const Atom& a : mol.atoms) {
        if (a.placed) queue.push_back(a.id);
    }

    if (queue.empty()) {
        const int root = pickRoot(mol, env);
        mol.atoms[root].pos = v3();
        mol.atoms[root].placed = true;
        queue.push_back(root);
    }

    std::size_t guard = 0;
    while (guard++ < n * 8) {
        if (queue.empty()) {
            // отдельная несвязанная часть молекулы — отодвигаем её в сторону
            Atom* next = nullptr;
            for (Atom& a : mol.atoms) {
                if (!a.placed) { next = &a; break; }
            }
            if (next == nullptr) break;
            const double shift = spanOf(mol) + 3;
            next->pos = v3(shift, 0, 0);
            next->placed = true;
            queue.push_back(next->id);
        }

        const int id = queue.front();
        queue.pop_front();
        const AtomEnvironment& e = env[id];

        std::vector<int> pending;
        for (int k : e.neighbors) {
            if (!mol.atoms[k].placed) pending.push_back(k);
        }
        if (pending.empty() && frameKnown[id]) continue;

        const Mat3 rot = resolveFrame(mol, env, id);
        frames[id] = rot;
        frameKnown[id] = true;

        for (int nbId : pending) {
            const int slot = static_cast<int>(
                std::find(e.neighbors.begin(), e.neighbors.end(), nbId) - e.neighbors.begin());
            const Vec3 dir = norm(applyMat(rot, e.localDirs[slot]));
            const double L = bondLength(mol.atoms[id].el, mol.atoms[nbId].el, e.bondOrders[slot]);
            mol.atoms[nbId].pos = add(mol.atoms[id].pos, mul(dir, L));
            mol.atoms[nbId].placed = true;
            queue.push_back(nbId);
        }
    }

    // атомы, до которых обход не добрался (защита от зависания)
    for (Atom& atom : mol.atoms) {
        if (!atom.placed) { atom.pos = v3(spanOf(mol) + 2, 0, 0); atom.placed = true; }
    }

    relax(mol, env, templatePlaced);

    // --- направления неподелённых пар в мировых координатах --------------------
    BuildResult result;
    result.frames = frames;
    for (std::size_t id = 0; id < n; id++) {
        if (!frameKnown[id]) {
            result.frames[id] = resolveFrame(mol, env, static_cast<int>(id));
            frameKnown[id] = true;
        }
        const AtomEnvironment& e = env[id];
        std::vector<Vec3> dirs;
        for (std::size_t k = e.neighbors.size(); k < e.localDirs.size(); k++) {
            dirs.push_back(norm(applyMat(result.frames[id], e.localDirs[k])));
        }
        result.lonePairDirs.push_back(dirs);
    }

    // Разворот по главным осям двигает координаты атомов, поэтому тот же поворот
    // надо применить и к системам координат атомов. В оригинале на TypeScript
    // этого не делалось, и направления неподелённых пар оставались в старой
    // системе: облака НЭП у воды рисовались не в тетраэдрических позициях,
    // а дипольный момент получал случайную величину. Расхождение согласовано
    // с заказчиком; на сверку с эталоном оно не влияет — там ни направлений
    // пар, ни величины μ нет, а признак полярности совпадает.
    const Mat3 orientation = orientByPrincipalAxes(mol);
    for (std::size_t id = 0; id < n; id++) {
        result.frames[id] = matMul(orientation, result.frames[id]);
        for (Vec3& d : result.lonePairDirs[id]) d = norm(applyMat(orientation, d));
    }

    return result;
}

}  // namespace chem
