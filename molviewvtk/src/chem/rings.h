#pragma once

// Поиск наименьшего набора наименьших циклов (SSSR).

#include <vector>

#include "chem/types.h"

namespace chem {

/**
 * Приём простой и надёжный: для каждой связи мысленно её удаляем и ищем
 * кратчайший путь между её концами. Путь плюс сама связь и есть цикл,
 * причём наименьший из проходящих через эту связь. Дальше отбираем
 * независимые циклы — столько, сколько даёт цикломатическое число графа.
 */
std::vector<Ring> findRings(const Molecule& mol);

/** Группировка циклов в конденсированные системы (имеющие общие атомы). */
std::vector<std::vector<Ring>> fuseRingSystems(const std::vector<Ring>& rings);

}  // namespace chem
