#pragma once

// База молекул: загрузка, поиск, изомеры, разбор пользовательского запроса.

#include <string>
#include <vector>

#include "chem/types.h"

namespace chem {

/** Запись базы вместе с вычисленными при загрузке полями. */
struct DbEntry : MoleculeMeta {
    std::string compositionKey;
    std::string formula;
    std::string formulaPretty;
    /** Название, синонимы, формула и SMILES одной строкой в нижнем регистре. */
    std::string searchText;
};

const std::vector<DbEntry>& database();

const DbEntry* entryById(const std::string& id);

struct Category {
    std::string name;
    std::vector<const DbEntry*> entries;
};

/** Категории в порядке появления в базе. */
std::vector<Category> categories();

struct Resolution {
    enum class Kind { Entries, Smiles, Error };
    enum class Reason { Name, Formula };

    Kind kind = Kind::Error;
    /** kind == Entries */
    std::vector<const DbEntry*> entries;
    Reason reason = Reason::Name;
    /** kind == Smiles */
    Molecule molecule;
    std::string smiles;
    /** kind == Error */
    std::string message;
    std::string hint;
};

/**
 * Разбор пользовательского запроса.
 *
 * Порядок попыток: название → брутто-формула → вхождение в название → SMILES.
 * Формула специально идёт раньше SMILES: на запрос «C2H6O» программа должна
 * показать ОБА изомера, а не молчаливо разобрать строку как SMILES.
 */
Resolution resolveQuery(const std::string& input);

/** Подсказки автодополнения. */
std::vector<const DbEntry*> suggestions(const std::string& input, std::size_t limit = 8);

/** Молекулы с тем же составом — изомеры. */
std::vector<const DbEntry*> isomersOf(const DbEntry& entry);

Molecule moleculeFromEntry(const DbEntry& entry);

/** Приведение строки к нижнему регистру (латиница и кириллица). */
std::string toLower(const std::string& s);

}  // namespace chem
