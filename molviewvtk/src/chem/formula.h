#pragma once

// Брутто-формулы и разбор формулы, введённой пользователем.

#include <string>
#include <utility>
#include <vector>

#include "chem/types.h"

namespace chem {

/** Состав: элемент → количество, в порядке первого появления в молекуле. */
std::vector<std::pair<std::string, int>> composition(const Molecule& mol);

double totalCharge(const Molecule& mol);
double molecularMass(const Molecule& mol);

/** Формула обычными символами: «C2H6O». Используется для поиска. */
std::string plainFormula(const Molecule& mol);

/** Формула с настоящими подстрочными цифрами и зарядом: «C₂H₆O», «NO₃⁻». */
std::string prettyFormula(const Molecule& mol);

/**
 * Разбор формулы, введённой пользователем: «C2H6O», «h2o», «CH3COOH».
 * Возвращает нормализованный ключ состава либо пустую строку, если это не формула.
 */
std::string parseFormulaQuery(const std::string& input);

/** Канонический ключ состава для сравнения формул: элементы по алфавиту. */
std::string compositionKey(const std::vector<std::pair<std::string, int>>& counts);

std::string moleculeCompositionKey(const Molecule& mol);

/** Нормализация подписи валентного угла: внешние атомы сортируются по алфавиту. */
std::string angleLabel(const std::string& outerA, const std::string& center,
                       const std::string& outerB);

}  // namespace chem
