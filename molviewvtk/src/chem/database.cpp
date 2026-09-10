#include "chem/database.h"

#include <cstdio>

#include <utility>

#include "chem/formula.h"
#include "chem/json.h"
#include "chem/smiles.h"

namespace chem {

namespace {

/** Разбор UTF-8 в кодовые точки: без этого кириллицу не привести к нижнему регистру. */
std::vector<unsigned> decodeUtf8(const std::string& s) {
    std::vector<unsigned> out;
    std::size_t i = 0;
    while (i < s.size()) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        unsigned code;
        int extra;
        if (c < 0x80) { code = c; extra = 0; }
        else if ((c & 0xE0) == 0xC0) { code = c & 0x1Fu; extra = 1; }
        else if ((c & 0xF0) == 0xE0) { code = c & 0x0Fu; extra = 2; }
        else { code = c & 0x07u; extra = 3; }
        i++;
        for (int k = 0; k < extra && i < s.size(); k++, i++) {
            code = (code << 6) | (static_cast<unsigned char>(s[i]) & 0x3Fu);
        }
        out.push_back(code);
    }
    return out;
}

void encodeUtf8(std::string& out, unsigned code) {
    if (code < 0x80) {
        out += static_cast<char>(code);
    } else if (code < 0x800) {
        out += static_cast<char>(0xC0 | (code >> 6));
        out += static_cast<char>(0x80 | (code & 0x3F));
    } else {
        out += static_cast<char>(0xE0 | (code >> 12));
        out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
        out += static_cast<char>(0x80 | (code & 0x3F));
    }
}

bool isSpace(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

std::string trim(const std::string& s) {
    std::size_t from = 0, to = s.size();
    while (from < to && isSpace(s[from])) from++;
    while (to > from && isSpace(s[to - 1])) to--;
    return s.substr(from, to - from);
}

/** Ключ поиска: обрезка, нижний регистр и «ё» вместо «е», как в оригинале. */
std::string normalizeQuery(const std::string& s) {
    const std::string lowered = toLower(trim(s));
    std::string out;
    for (unsigned code : decodeUtf8(lowered)) {
        encodeUtf8(out, code == 0x451 ? 0x435 : code);  // ё → е
    }
    return out;
}

/** Название, синонимы, формула и SMILES одной строкой — по ней идёт поиск подстроки. */
std::string buildSearchText(const MoleculeMeta& meta, const std::string& formula) {
    std::string joined = meta.name;
    for (const std::string& s : meta.syn) joined += " " + s;
    joined += " " + formula + " " + meta.smiles;
    return normalizeQuery(joined);
}

std::vector<DbEntry> buildIndex() {
    std::vector<DbEntry> entries;
    const Json root = Json::parseFile(findDataFile("molecules.json"));

    for (const Json& item : root["molecules"].items()) {
        MoleculeMeta meta;
        meta.id = item["id"].str();
        meta.name = item["name"].str();
        meta.smiles = item["smiles"].str();
        meta.cat = item["cat"].str();
        for (const Json& s : item["syn"].items()) meta.syn.push_back(s.str());
        for (const auto& field : item["exp"].fields()) {
            std::vector<double> values;
            if (field.second.isArray()) {
                for (const Json& v : field.second.items()) values.push_back(v.number());
            } else {
                values.push_back(field.second.number());
            }
            meta.exp.emplace_back(field.first, values);
        }
        if (item["dipole"].isNumber()) {
            meta.dipole = item["dipole"].number();
            meta.hasDipole = true;
        }
        meta.note = item["note"].str();

        Molecule mol;
        try {
            mol = parseSmiles(meta.smiles);
        } catch (const SmilesError& err) {
            // молекула с ошибкой в базе не должна ронять всю программу
            std::fprintf(stderr, "Молекула «%s» пропущена: %s\n", meta.name.c_str(), err.what());
            continue;
        }

        DbEntry entry;
        static_cast<MoleculeMeta&>(entry) = meta;
        entry.compositionKey = moleculeCompositionKey(mol);
        entry.formula = plainFormula(mol);
        entry.formulaPretty = prettyFormula(mol);
        entry.searchText = buildSearchText(meta, entry.formula);
        entries.push_back(std::move(entry));
    }
    return entries;
}

/** Похоже ли на брутто-формулу: только буквы и цифры, первый символ — буква. */
bool looksLikeFormula(const std::string& query) {
    if (query.empty()) return false;
    const char first = query[0];
    if (!((first >= 'A' && first <= 'Z') || (first >= 'a' && first <= 'z'))) return false;
    for (char c : query) {
        const bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
        if (!ok) return false;
    }
    return true;
}

}  // namespace

std::string toLower(const std::string& s) {
    std::string out;
    for (unsigned code : decodeUtf8(s)) {
        if (code >= 'A' && code <= 'Z') code += 32;                  // латиница
        else if (code >= 0x410 && code <= 0x42F) code += 32;         // А–Я
        else if (code == 0x401) code = 0x451;                        // Ё
        encodeUtf8(out, code);
    }
    return out;
}

const std::vector<DbEntry>& database() {
    static const std::vector<DbEntry> index = buildIndex();
    return index;
}

const DbEntry* entryById(const std::string& id) {
    for (const DbEntry& e : database()) {
        if (e.id == id) return &e;
    }
    return nullptr;
}

std::vector<Category> categories() {
    std::vector<Category> result;
    for (const DbEntry& e : database()) {
        auto it = result.begin();
        while (it != result.end() && it->name != e.cat) ++it;
        if (it == result.end()) {
            result.push_back(Category{e.cat, {&e}});
        } else {
            it->entries.push_back(&e);
        }
    }
    return result;
}

Resolution resolveQuery(const std::string& input) {
    Resolution result;
    const std::string query = trim(input);
    if (query.empty()) {
        result.kind = Resolution::Kind::Error;
        result.message = "Введите формулу, название или SMILES";
        return result;
    }

    const std::vector<DbEntry>& db = database();
    const std::string needle = normalizeQuery(query);

    // 1. точное совпадение названия или синонима
    std::vector<const DbEntry*> exact;
    for (const DbEntry& e : db) {
        bool match = normalizeQuery(e.name) == needle;
        if (!match) {
            for (const std::string& s : e.syn) {
                if (normalizeQuery(s) == needle) { match = true; break; }
            }
        }
        if (match) exact.push_back(&e);
    }
    if (!exact.empty()) {
        result.kind = Resolution::Kind::Entries;
        result.entries = exact;
        result.reason = Resolution::Reason::Name;
        return result;
    }

    // 2. брутто-формула — сюда попадают все изомеры
    const std::string key = parseFormulaQuery(query);
    if (!key.empty()) {
        std::vector<const DbEntry*> byFormula;
        for (const DbEntry& e : db) {
            if (e.compositionKey == key) byFormula.push_back(&e);
        }
        if (!byFormula.empty()) {
            result.kind = Resolution::Kind::Entries;
            result.entries = byFormula;
            result.reason = Resolution::Reason::Formula;
            return result;
        }
    }

    // 3. вхождение в название
    std::vector<const DbEntry*> partial;
    for (const DbEntry& e : db) {
        if (e.searchText.find(needle) != std::string::npos) partial.push_back(&e);
    }
    if (!partial.empty()) {
        result.kind = Resolution::Kind::Entries;
        result.entries = partial;
        result.reason = Resolution::Reason::Name;
        return result;
    }

    // 4. SMILES
    try {
        result.molecule = parseSmiles(query);
        result.kind = Resolution::Kind::Smiles;
        result.smiles = query;
        return result;
    } catch (const SmilesError& err) {
        result.kind = Resolution::Kind::Error;
        if (!key.empty()) {
            result.message = "Формулы " + query + " нет в справочнике";
            result.hint = "Строение можно задать напрямую в виде SMILES — например, CCO для этанола "
                          "или c1ccccc1 для бензола. Подсказка по синтаксису — кнопка «?» справа вверху.";
            return result;
        }
        // запрос без скобок и знаков связи скорее всего задумывался как формула
        if (looksLikeFormula(query)) {
            result.message = "Не похоже ни на формулу, ни на название: " + toLower(err.what());
            result.hint = "Если это брутто-формула, проверьте символы элементов: они пишутся "
                          "с заглавной буквы (Cl, Br, Xe). Можно также ввести название по-русски "
                          "или строение в виде SMILES.";
        } else {
            result.message = std::string("Не удалось разобрать запрос: ") + err.what();
            result.hint = "Введите брутто-формулу (H2O, C2H6O), название («вода», «бензол») "
                          "или строение в виде SMILES (CCO, c1ccccc1).";
        }
        return result;
    }
}

std::vector<const DbEntry*> suggestions(const std::string& input, std::size_t limit) {
    const std::string needle = normalizeQuery(input);
    if (needle.empty()) return {};

    const std::vector<DbEntry>& db = database();
    std::vector<const DbEntry*> starts;
    std::vector<const DbEntry*> contains;
    for (const DbEntry& e : db) {
        if (normalizeQuery(e.name).rfind(needle, 0) == 0) starts.push_back(&e);
        else if (e.searchText.find(needle) != std::string::npos) contains.push_back(&e);
    }

    std::vector<const DbEntry*> result = starts;
    result.insert(result.end(), contains.begin(), contains.end());
    if (result.size() > limit) result.resize(limit);
    return result;
}

std::vector<const DbEntry*> isomersOf(const DbEntry& entry) {
    std::vector<const DbEntry*> result;
    for (const DbEntry& e : database()) {
        if (e.compositionKey == entry.compositionKey && e.id != entry.id) result.push_back(&e);
    }
    return result;
}

Molecule moleculeFromEntry(const DbEntry& entry) {
    Molecule mol = parseSmiles(entry.smiles);
    mol.meta = &entry;
    return mol;
}

}  // namespace chem
