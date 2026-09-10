#include "chem/resonance.h"

#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

#include "chem/smiles.h"

namespace chem {

void applyResonance(Molecule& mol) {
    const std::vector<std::vector<int>> nb = neighborLists(mol);

    for (std::size_t centerId = 0; centerId < mol.atoms.size(); centerId++) {
        const int center = static_cast<int>(centerId);
        if (nb[center].size() < 2) continue;

        // концевые соседи, сгруппированные по элементу (порядок групп — как у соседей)
        std::vector<std::pair<std::string, std::vector<int>>> groups;
        for (int nbId : nb[center]) {
            const Atom& neighbor = mol.atoms[nbId];
            if (neighbor.el == "H") continue;
            if (nb[nbId].size() != 1) continue;  // не концевой — не равноценен

            auto it = groups.begin();
            while (it != groups.end() && it->first != neighbor.el) ++it;
            if (it == groups.end()) {
                groups.emplace_back(neighbor.el, std::vector<int>{nbId});
            } else {
                it->second.push_back(nbId);
            }
        }

        for (const auto& group : groups) {
            const std::vector<int>& members = group.second;
            if (members.size() < 2) continue;

            std::vector<Bond*> bonds;
            for (int id : members) bonds.push_back(findBond(mol, center, id));

            bool anyAromatic = false;
            for (const Bond* b : bonds) if (b->aromatic) anyAromatic = true;
            if (anyAromatic) continue;

            bool sameOrder = true;
            bool sameCharge = true;
            for (std::size_t k = 0; k < members.size(); k++) {
                if (std::fabs(bonds[k]->order - bonds[0]->order) >= 1e-9) sameOrder = false;
                if (mol.atoms[members[k]].charge != mol.atoms[members[0]].charge) sameCharge = false;
            }
            if (sameOrder && sameCharge) continue;  // равноценны и без усреднения

            double sumOrder = 0;
            double sumCharge = 0;
            for (std::size_t k = 0; k < members.size(); k++) {
                sumOrder += bonds[k]->order;
                sumCharge += mol.atoms[members[k]].charge;
            }
            const double avgOrder = sumOrder / static_cast<double>(members.size());
            const double avgCharge = sumCharge / static_cast<double>(members.size());

            for (Bond* b : bonds) { b->order = avgOrder; b->delocalized = true; }
            for (int id : members) mol.atoms[id].charge = avgCharge;
        }
    }
}

std::string formatBondOrder(double order) {
    static const std::pair<double, const char*> table[] = {
        {1.0, "1"}, {4.0 / 3.0, "1⅓"}, {1.25, "1¼"}, {1.5, "1½"}, {5.0 / 3.0, "1⅔"},
        {2.0, "2"}, {2.5, "2½"}, {3.0, "3"},
    };
    for (const auto& row : table) {
        if (std::fabs(order - row.first) < 0.02) return row.second;
    }
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.2f", order);
    std::string s(buffer);
    if (!s.empty() && s.back() == '0') s.pop_back();
    return s;
}

}  // namespace chem
