#include "chem/rings.h"

#include <algorithm>
#include <deque>
#include <map>
#include <set>
#include <string>

#include "chem/smiles.h"

namespace chem {

namespace {

/** Ключ неупорядоченной пары атомов. */
std::string pairKey(int a, int b) {
    return a < b ? std::to_string(a) + "-" + std::to_string(b)
                 : std::to_string(b) + "-" + std::to_string(a);
}

std::vector<std::string> ringBondKeys(const std::vector<int>& ring) {
    std::vector<std::string> keys;
    for (std::size_t i = 0; i < ring.size(); i++) {
        keys.push_back(pairKey(ring[i], ring[(i + 1) % ring.size()]));
    }
    return keys;
}

/** Кратчайший путь из from в to обходом в ширину, минуя запрещённую связь. */
std::vector<int> shortestPath(const std::vector<std::vector<int>>& nb,
                              const std::set<int>& isHeavy,
                              int from, int to, const Bond& forbidden) {
    std::map<int, int> prev;
    std::deque<int> queue{from};
    std::set<int> visited{from};

    while (!queue.empty()) {
        const int cur = queue.front();
        queue.pop_front();
        if (cur == to) break;
        for (int next : nb[cur]) {
            if (isHeavy.count(next) == 0) continue;
            const bool isForbidden =
                (cur == forbidden.a && next == forbidden.b) || (cur == forbidden.b && next == forbidden.a);
            if (isForbidden) continue;
            if (visited.count(next) != 0) continue;
            visited.insert(next);
            prev.emplace(next, cur);
            queue.push_back(next);
        }
    }

    if (visited.count(to) == 0) return {};

    std::vector<int> path;
    int cur = to;
    for (;;) {
        path.push_back(cur);
        if (cur == from) break;
        const auto it = prev.find(cur);
        if (it == prev.end()) break;
        cur = it->second;
    }
    if (path.size() < 3) return {};
    std::reverse(path.begin(), path.end());
    return path;
}

int countComponents(const std::vector<std::vector<int>>& nb, const std::vector<int>& heavy,
                    const std::set<int>& isHeavy) {
    std::set<int> visited;
    int count = 0;
    for (int start : heavy) {
        if (visited.count(start) != 0) continue;
        count++;
        std::vector<int> stack{start};
        visited.insert(start);
        while (!stack.empty()) {
            const int cur = stack.back();
            stack.pop_back();
            for (int next : nb[cur]) {
                if (isHeavy.count(next) == 0 || visited.count(next) != 0) continue;
                visited.insert(next);
                stack.push_back(next);
            }
        }
    }
    return count;
}

Ring describeRing(const Molecule& mol, const std::vector<int>& atoms) {
    bool aromatic = true;
    for (std::size_t i = 0; i < atoms.size(); i++) {
        const Bond* bond = findBond(mol, atoms[i], atoms[(i + 1) % atoms.size()]);
        if (bond == nullptr || !bond->aromatic) { aromatic = false; break; }
    }
    Ring ring;
    ring.atoms = atoms;
    ring.aromatic = aromatic;
    ring.planar = aromatic;
    return ring;
}

}  // namespace

std::vector<Ring> findRings(const Molecule& mol) {
    const std::vector<std::vector<int>> nb = neighborLists(mol);

    std::vector<int> heavy;
    std::set<int> isHeavy;
    for (const Atom& a : mol.atoms) {
        if (a.el != "H") { heavy.push_back(a.id); isHeavy.insert(a.id); }
    }

    std::vector<std::vector<int>> candidates;
    for (const Bond& bond : mol.bonds) {
        if (isHeavy.count(bond.a) == 0 || isHeavy.count(bond.b) == 0) continue;
        const std::vector<int> path = shortestPath(nb, isHeavy, bond.a, bond.b, bond);
        if (path.size() >= 3) candidates.push_back(path);
    }

    // уникальные циклы по набору атомов, от меньших к большим;
    // сортировка устойчивая — при равной длине порядок задан порядком связей
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const std::vector<int>& x, const std::vector<int>& y) {
                         return x.size() < y.size();
                     });

    std::set<std::string> seen;
    std::vector<std::vector<int>> unique;
    for (const std::vector<int>& c : candidates) {
        std::vector<int> sorted = c;
        std::sort(sorted.begin(), sorted.end());
        std::string key;
        for (std::size_t k = 0; k < sorted.size(); k++) {
            if (k > 0) key += ',';
            key += std::to_string(sorted[k]);
        }
        if (seen.count(key) != 0) continue;
        seen.insert(key);
        unique.push_back(c);
    }

    // цикломатическое число: сколько независимых циклов есть в графе
    int heavyBondCount = 0;
    for (const Bond& b : mol.bonds) {
        if (isHeavy.count(b.a) != 0 && isHeavy.count(b.b) != 0) heavyBondCount++;
    }
    const int components = countComponents(nb, heavy, isHeavy);
    const int cycleCount = heavyBondCount - static_cast<int>(heavy.size()) + components;

    std::vector<std::vector<int>> chosen;
    std::set<std::string> coveredBonds;
    for (const std::vector<int>& c : unique) {
        if (static_cast<int>(chosen.size()) >= cycleCount) break;
        const std::vector<std::string> bondKeys = ringBondKeys(c);
        bool hasNew = false;
        for (const std::string& k : bondKeys) {
            if (coveredBonds.count(k) == 0) { hasNew = true; break; }
        }
        if (hasNew || chosen.empty()) {
            chosen.push_back(c);
            for (const std::string& k : bondKeys) coveredBonds.insert(k);
        }
    }

    std::vector<Ring> result;
    for (const std::vector<int>& atoms : chosen) result.push_back(describeRing(mol, atoms));
    return result;
}

std::vector<std::vector<Ring>> fuseRingSystems(const std::vector<Ring>& rings) {
    std::vector<std::vector<Ring>> groups;
    std::vector<bool> assigned(rings.size(), false);

    for (std::size_t i = 0; i < rings.size(); i++) {
        if (assigned[i]) continue;
        std::vector<Ring> group{rings[i]};
        assigned[i] = true;
        bool grew = true;
        while (grew) {
            grew = false;
            std::set<int> inGroup;
            for (const Ring& r : group) {
                for (int a : r.atoms) inGroup.insert(a);
            }
            for (std::size_t j = 0; j < rings.size(); j++) {
                if (assigned[j]) continue;
                bool shares = false;
                for (int a : rings[j].atoms) {
                    if (inGroup.count(a) != 0) { shares = true; break; }
                }
                if (shares) {
                    group.push_back(rings[j]);
                    assigned[j] = true;
                    grew = true;
                }
            }
        }
        groups.push_back(group);
    }
    return groups;
}

}  // namespace chem
