#include "chem/formula.h"

#include <algorithm>
#include <cmath>

#include "chem/periodic.h"

namespace chem {

namespace {

const char* const SUB[] = {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};

std::string toSub(int n) {
    std::string out;
    for (char c : std::to_string(n)) out += SUB[c - '0'];
    return out;
}

/** Надстрочные знаки для заряда. */
std::string toSup(const std::string& s) {
    std::string out;
    for (char c : s) {
        switch (c) {
            case '0': out += "⁰"; break;
            case '1': out += "¹"; break;
            case '2': out += "²"; break;
            case '3': out += "³"; break;
            case '4': out += "⁴"; break;
            case '5': out += "⁵"; break;
            case '6': out += "⁶"; break;
            case '7': out += "⁷"; break;
            case '8': out += "⁸"; break;
            case '9': out += "⁹"; break;
            case '+': out += "⁺"; break;
            case '-': out += "⁻"; break;
            default: out += c; break;
        }
    }
    return out;
}

/**
 * Электроотрицательностный ряд элементов (IUPAC, Red Book, табл. VI).
 * В формуле неорганического вещества первым пишется элемент, стоящий в ряду
 * левее (более электроположительный), поэтому получается SO₂, PCl₃, XeF₂,
 * а кислород и фтор оказываются в конце. Водород стоит между азотом и серой:
 * отсюда NH₃, PH₃, SiH₄, но H₂O, H₂S, HCl.
 */
const char* const EN_ORDER[] = {
    "Rn", "Xe", "Kr", "Ar", "Ne", "He",
    "Fr", "Cs", "Rb", "K",  "Na", "Li",
    "Ra", "Ba", "Sr", "Ca", "Mg", "Be",
    "Tl", "In", "Ga", "Al", "B",
    "Pb", "Sn", "Ge", "Si", "C",
    "Bi", "Sb", "As", "P",  "N",
    "H",
    "Po", "Te", "Se", "S",
    "At", "I",  "Br", "Cl",
    "O",  "F",
};

/** Место элемента в ряду; для элементов вне списка — −1, то есть в самое начало. */
int enRank(const std::string& symbol) {
    const int n = static_cast<int>(sizeof(EN_ORDER) / sizeof(EN_ORDER[0]));
    for (int i = 0; i < n; ++i) {
        if (symbol == EN_ORDER[i]) return i;
    }
    return -1;
}

/**
 * Порядок элементов в брутто-формуле.
 *
 * Для органики — система Хилла (C, H, остальное по алфавиту): C₂H₆O, C₆H₆.
 * Для неорганики — электроотрицательностный ряд: SO₂, PCl₃, NH₃, а не
 * алфавитные O₂S, Cl₃P, H₃N. У кислородсодержащих кислот водород по традиции
 * выносится вперёд: H₂SO₄, HNO₃, H₃PO₄. Немногие исключения (OH⁻ вместо HO⁻,
 * HCN вместо хилловского CHN) перечислены явно.
 *
 * Сравнение строк обычное побайтовое: символы элементов — заглавная буква плюс
 * необязательная строчная, и для такого набора оно совпадает с localeCompare.
 */
std::vector<std::string> formulaOrder(const std::vector<std::pair<std::string, int>>& counts) {
    const std::string key = compositionKey(counts);
    if (key == "H1O1") return {"O", "H"};
    if (key == "C1H1N1") return {"H", "C", "N"};

    std::vector<std::string> symbols;
    for (const auto& item : counts) symbols.push_back(item.first);

    const auto present = [&counts](const std::string& s) {
        for (const auto& item : counts) {
            if (item.first == s) return true;
        }
        return false;
    };

    if (present("C")) {
        std::vector<std::string> rest;
        for (const std::string& s : symbols) {
            if (s != "C" && s != "H") rest.push_back(s);
        }
        std::sort(rest.begin(), rest.end());

        std::vector<std::string> order{"C"};
        if (present("H")) order.push_back("H");
        order.insert(order.end(), rest.begin(), rest.end());
        return order;
    }

    std::sort(symbols.begin(), symbols.end(),
              [](const std::string& a, const std::string& b) {
                  const int ra = enRank(a);
                  const int rb = enRank(b);
                  return ra != rb ? ra < rb : a < b;
              });

    // кислота: есть и водород, и кислород, и ещё хотя бы один элемент
    if (symbols.size() > 2 && present("H") && present("O")) {
        std::vector<std::string> order{"H"};
        for (const std::string& s : symbols) {
            if (s != "H") order.push_back(s);
        }
        return order;
    }
    return symbols;
}

int countOf(const std::vector<std::pair<std::string, int>>& counts, const std::string& symbol) {
    for (const auto& item : counts) {
        if (item.first == symbol) return item.second;
    }
    return 0;
}

bool isDigit(char c) { return c >= '0' && c <= '9'; }
bool isLower(char c) { return c >= 'a' && c <= 'z'; }
bool isUpper(char c) { return c >= 'A' && c <= 'Z'; }
char upper(char c) { return isLower(c) ? static_cast<char>(c - 'a' + 'A') : c; }
char lower(char c) { return isUpper(c) ? static_cast<char>(c - 'A' + 'a') : c; }

}  // namespace

std::vector<std::pair<std::string, int>> composition(const Molecule& mol) {
    std::vector<std::pair<std::string, int>> counts;
    for (const Atom& atom : mol.atoms) {
        auto it = counts.begin();
        while (it != counts.end() && it->first != atom.el) ++it;
        if (it == counts.end()) counts.emplace_back(atom.el, 1);
        else it->second++;
    }
    return counts;
}

double totalCharge(const Molecule& mol) {
    double s = 0;
    for (const Atom& a : mol.atoms) s += a.charge;
    return s;
}

double molecularMass(const Molecule& mol) {
    double s = 0;
    for (const Atom& a : mol.atoms) s += element(a.el).mass;
    return s;
}

std::string plainFormula(const Molecule& mol) {
    const auto counts = composition(mol);
    std::string out;
    for (const std::string& s : formulaOrder(counts)) {
        const int n = countOf(counts, s);
        out += s;
        if (n > 1) out += std::to_string(n);
    }
    return out;
}

std::string prettyFormula(const Molecule& mol) {
    const auto counts = composition(mol);
    std::string s;
    for (const std::string& sym : formulaOrder(counts)) {
        const int n = countOf(counts, sym);
        s += sym;
        if (n > 1) s += toSub(n);
    }
    const double q = totalCharge(mol);
    if (q != 0) {
        // после усреднения зарядов сумма бывает вида 0,9999999999 — целое берём округлением
        const double magnitude = std::fabs(q);
        std::string prefix;
        if (magnitude > 1) prefix = std::to_string(static_cast<int>(jsRound(magnitude)));
        s += toSup(prefix + (q > 0 ? "+" : "-"));
    }
    return s;
}

std::string parseFormulaQuery(const std::string& input) {
    // выбрасываем пробелы, скобки, точку и знак умножения кристаллогидрата
    std::string text;
    for (std::size_t k = 0; k < input.size(); k++) {
        const unsigned char c = static_cast<unsigned char>(input[k]);
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '(' || c == ')' || c == '.') continue;
        if (c == 0xC2 && k + 1 < input.size() && static_cast<unsigned char>(input[k + 1]) == 0xB7) {
            k++;  // «·» в UTF-8 занимает два байта
            continue;
        }
        text += input[k];
    }
    if (text.empty()) return "";

    std::vector<std::pair<std::string, int>> counts;
    const auto bump = [&counts](const std::string& symbol, int n) {
        auto it = counts.begin();
        while (it != counts.end() && it->first != symbol) ++it;
        if (it == counts.end()) counts.emplace_back(symbol, n);
        else it->second += n;
    };

    std::size_t i = 0;
    while (i < text.size()) {
        // элемент: заглавная + необязательная строчная
        std::string symbol;
        if (i + 1 < text.size()
            && (isUpper(text[i]) || isLower(text[i])) && isLower(text[i + 1])) {
            const std::string candidate = std::string(1, upper(text[i])) + lower(text[i + 1]);
            if (isKnownElement(candidate)) { symbol = candidate; i += 2; }
        }
        if (symbol.empty()) {
            const std::string one(1, upper(text[i]));
            if (!isKnownElement(one)) return "";
            symbol = one;
            i += 1;
        }
        std::string digits;
        while (i < text.size() && isDigit(text[i])) digits += text[i++];
        const int n = digits.empty() ? 1 : std::stoi(digits);
        bump(symbol, n);
    }
    if (counts.empty()) return "";
    return compositionKey(counts);
}

std::string compositionKey(const std::vector<std::pair<std::string, int>>& counts) {
    std::vector<std::pair<std::string, int>> filtered;
    for (const auto& item : counts) {
        if (item.second > 0) filtered.push_back(item);
    }
    std::sort(filtered.begin(), filtered.end(),
              [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
                  return a.first < b.first;
              });
    std::string out;
    for (const auto& item : filtered) out += item.first + std::to_string(item.second);
    return out;
}

std::string moleculeCompositionKey(const Molecule& mol) {
    return compositionKey(composition(mol));
}

std::string angleLabel(const std::string& outerA, const std::string& center,
                       const std::string& outerB) {
    const std::string& x = outerA <= outerB ? outerA : outerB;
    const std::string& y = outerA <= outerB ? outerB : outerA;
    return x + "-" + center + "-" + y;
}

}  // namespace chem
