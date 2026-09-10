#pragma once

// Разбор строки SMILES в граф молекулы.
//
// Поддерживается: органическое подмножество (B, C, N, O, P, S, F, Cl, Br, I),
// атомы в квадратных скобках с зарядом и явным числом водородов,
// кратные связи - = # :, ветвления, замыкания циклов (в том числе %nn),
// ароматические атомы в нижнем регистре, разрыв «.».
//
// Не поддерживается (сознательно, для школьного проекта не нужно):
// стереохимия @/@@ и /\ — эти символы просто пропускаются.

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "chem/types.h"

namespace chem {

class SmilesError : public std::runtime_error {
public:
    SmilesError(const std::string& message, std::size_t position)
        : std::runtime_error(message), position(position) {}
    std::size_t position;
};

Molecule parseSmiles(const std::string& input);

/** Списки соседей для каждого атома. */
std::vector<std::vector<int>> neighborLists(const Molecule& mol);

/** Сумма кратностей связей у атома. */
double bondOrderSum(const Molecule& mol, int id);

const Bond* findBond(const Molecule& mol, int a, int b);
Bond* findBond(Molecule& mol, int a, int b);

}  // namespace chem
