#include "chem/smiles.h"

#include <set>
#include <utility>

#include "chem/periodic.h"

namespace chem {

namespace {

const std::set<std::string>& organicSubset() {
    static const std::set<std::string> set = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"};
    return set;
}

bool isAromaticSymbol(char c) {
    return c == 'b' || c == 'c' || c == 'n' || c == 'o' || c == 'p' || c == 's';
}

bool isDigit(char c) { return c >= '0' && c <= '9'; }
bool isLower(char c) { return c >= 'a' && c <= 'z'; }
bool isUpper(char c) { return c >= 'A' && c <= 'Z'; }
bool isSpace(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v'; }

std::string trim(const std::string& s) {
    std::size_t from = 0;
    std::size_t to = s.size();
    while (from < to && isSpace(s[from])) from++;
    while (to > from && isSpace(s[to - 1])) to--;
    return s.substr(from, to - from);
}

char upper(char c) { return isLower(c) ? static_cast<char>(c - 'a' + 'A') : c; }

/** Атом «как написан», до развёртывания неявных водородов. */
struct RawAtom {
    std::string el;
    int charge = 0;
    bool aromatic = false;
    int hCount = 0;
    bool explicitH = false;
};

/** Открытое замыкание цикла: ждёт второго появления своего номера. */
struct OpenRing {
    int atom = 0;
    bool hasOrder = false;
    double order = 0.0;
    bool aromatic = false;
};

RawAtom parseBracketAtom(const std::string& body, std::size_t position) {
    std::size_t p = 0;
    // изотоп — пропускаем
    while (p < body.size() && isDigit(body[p])) p++;

    std::string el;
    bool aromatic = false;
    if (p < body.size() && isLower(body[p]) && isAromaticSymbol(body[p])) {
        aromatic = true;
        el = std::string(1, upper(body[p]));
        p++;
    } else if (p < body.size() && isUpper(body[p])) {
        el = std::string(1, body[p]);
        p++;
        if (p < body.size() && isLower(body[p]) && isKnownElement(el + body[p])) {
            el += body[p];
            p++;
        }
    } else {
        throw SmilesError("Не удалось определить элемент в «[" + body + "]»", position);
    }

    // хиральность — пропускаем
    while (p < body.size() && body[p] == '@') p++;

    int hCount = 0;
    if (p < body.size() && body[p] == 'H') {
        p++;
        std::string digits;
        while (p < body.size() && isDigit(body[p])) digits += body[p++];
        hCount = digits.empty() ? 1 : std::stoi(digits);
    }

    int charge = 0;
    while (p < body.size() && (body[p] == '+' || body[p] == '-')) {
        const int sign = body[p] == '+' ? 1 : -1;
        const char signChar = body[p];
        p++;
        std::string digits;
        while (p < body.size() && isDigit(body[p])) digits += body[p++];
        if (!digits.empty()) {
            charge += sign * std::stoi(digits);
        } else {
            charge += sign;
            while (p < body.size() && body[p] == signChar) { charge += sign; p++; }
        }
    }

    // класс атома «:12» — пропускаем
    if (p < body.size() && body[p] == ':') {
        p++;
        while (p < body.size() && isDigit(body[p])) p++;
    }

    if (p != body.size()) {
        throw SmilesError("Лишние символы в «[" + body + "]»", position);
    }
    if (!isKnownElement(el)) {
        throw SmilesError("Элемент «" + el + "» отсутствует в справочнике программы", position);
    }

    RawAtom atom;
    atom.el = el;
    atom.charge = charge;
    atom.aromatic = aromatic;
    atom.hCount = hCount;
    atom.explicitH = true;
    return atom;
}

/** Развёртывание неявных водородов в явные атомы. */
Molecule buildMolecule(const std::vector<RawAtom>& raw, const std::vector<Bond>& bonds) {
    std::vector<double> bondSum(raw.size(), 0.0);
    // для ароматических атомов: число связей плюс «лишняя» кратность
    // неароматических связей — см. пояснение в implicitHydrogens
    std::vector<double> aromaticSum(raw.size(), 0.0);
    for (const Bond& b : bonds) {
        bondSum[b.a] += b.order;
        bondSum[b.b] += b.order;
        const double extra = b.aromatic ? 0.0 : b.order - 1.0;
        aromaticSum[b.a] += 1.0 + extra;
        aromaticSum[b.b] += 1.0 + extra;
    }

    Molecule mol;
    mol.atoms.reserve(raw.size());
    for (std::size_t id = 0; id < raw.size(); id++) {
        const RawAtom& r = raw[id];
        Atom atom;
        atom.id = static_cast<int>(id);
        atom.el = r.el;
        atom.charge = r.charge;
        atom.aromatic = r.aromatic;
        atom.hCount = r.explicitH
            ? r.hCount
            : (r.aromatic
                   ? implicitHydrogens(r.el, aromaticSum[id] + 1.0, r.charge, true)
                   : implicitHydrogens(r.el, bondSum[id], r.charge));
        mol.atoms.push_back(atom);
    }

    mol.bonds = bonds;
    const std::size_t heavyCount = mol.atoms.size();
    for (std::size_t id = 0; id < heavyCount; id++) {
        const int n = mol.atoms[id].hCount;
        for (int k = 0; k < n; k++) {
            Atom h;
            h.id = static_cast<int>(mol.atoms.size());
            h.el = "H";
            h.fromImplicitH = true;
            mol.atoms.push_back(h);

            Bond b;
            b.a = static_cast<int>(id);
            b.b = h.id;
            b.order = 1.0;
            mol.bonds.push_back(b);
        }
    }
    return mol;
}

}  // namespace

Molecule parseSmiles(const std::string& input) {
    const std::string src = trim(input);
    if (src.empty()) throw SmilesError("Пустая строка", 0);

    std::vector<RawAtom> raw;
    std::vector<Bond> bonds;
    // порядок открытия важен: в сообщении об ошибке называется первый
    // незамкнутый цикл, а не цикл с наименьшим номером
    std::vector<std::pair<int, OpenRing>> rings;
    std::vector<int> branchStack;

    std::size_t i = 0;
    int prev = -1;  // -1 означает «предыдущего атома нет» (null в оригинале)
    bool hasPendingOrder = false;
    double pendingOrder = 0.0;
    bool pendingAromatic = false;

    const auto fail = [&i](const std::string& message) {
        throw SmilesError(message, i);
    };

    // связь с уже существующим атомом (замыкание цикла)
    const auto linkTo = [&](int target, bool hasOrder, double order, bool aromatic, bool ringClosure) {
        if (prev < 0) return;
        const bool bothAromatic = raw[prev].aromatic && raw[target].aromatic;
        Bond b;
        b.a = prev;
        b.b = target;
        if (!hasOrder) {
            b.aromatic = bothAromatic;
            b.order = bothAromatic ? 1.5 : 1.0;
        } else {
            b.aromatic = aromatic;
            b.order = order;
        }
        b.ringClosure = ringClosure;
        bonds.push_back(b);
    };

    const auto pushAtom = [&](const RawAtom& a) {
        const int id = static_cast<int>(raw.size());
        raw.push_back(a);
        if (prev >= 0) {
            const bool bothAromatic = raw[prev].aromatic && a.aromatic;
            Bond b;
            b.a = prev;
            b.b = id;
            b.order = hasPendingOrder ? pendingOrder : (bothAromatic ? 1.5 : 1.0);
            b.aromatic = hasPendingOrder ? pendingAromatic : bothAromatic;
            bonds.push_back(b);
        }
        hasPendingOrder = false;
        pendingAromatic = false;
        prev = id;
    };

    while (i < src.size()) {
        const char ch = src[i];

        // --- ветвления -------------------------------------------------------
        if (ch == '(') {
            if (prev < 0) fail("Ветвление «(» не может открывать формулу");
            branchStack.push_back(prev);
            i++;
            continue;
        }
        if (ch == ')') {
            if (branchStack.empty()) fail("Лишняя закрывающая скобка «)»");
            prev = branchStack.back();
            branchStack.pop_back();
            i++;
            continue;
        }

        // --- символы связей --------------------------------------------------
        if (ch == '-') { hasPendingOrder = true; pendingOrder = 1.0; pendingAromatic = false; i++; continue; }
        if (ch == '=') { hasPendingOrder = true; pendingOrder = 2.0; pendingAromatic = false; i++; continue; }
        if (ch == '#') { hasPendingOrder = true; pendingOrder = 3.0; pendingAromatic = false; i++; continue; }
        if (ch == ':') { hasPendingOrder = true; pendingOrder = 1.5; pendingAromatic = true; i++; continue; }
        if (ch == '/' || ch == '\\') { i++; continue; }  // стереохимия двойной связи — игнорируем
        if (ch == '~') { hasPendingOrder = true; pendingOrder = 1.0; i++; continue; }

        // --- разрыв цепи -----------------------------------------------------
        if (ch == '.') { prev = -1; hasPendingOrder = false; i++; continue; }

        // --- замыкания циклов ------------------------------------------------
        if (ch == '%' || isDigit(ch)) {
            int label = 0;
            if (ch == '%') {
                const std::string digits = src.substr(i + 1, 2);
                if (digits.size() != 2 || !isDigit(digits[0]) || !isDigit(digits[1])) {
                    fail("После «%» должны идти две цифры номера цикла");
                }
                label = std::stoi(digits);
                i += 3;
            } else {
                label = ch - '0';
                i += 1;
            }
            if (prev < 0) fail("Номер цикла не может стоять до первого атома");

            auto open = rings.begin();
            while (open != rings.end() && open->first != label) ++open;
            if (open == rings.end()) {
                OpenRing r;
                r.atom = prev;
                r.hasOrder = hasPendingOrder;
                r.order = pendingOrder;
                r.aromatic = pendingAromatic;
                rings.emplace_back(label, r);
            } else {
                if (open->second.atom == prev) {
                    fail("Цикл " + std::to_string(label) + " замыкается сам на себя");
                }
                const bool hasOrder = hasPendingOrder || open->second.hasOrder;
                const double order = hasPendingOrder ? pendingOrder : open->second.order;
                const bool aromatic = hasPendingOrder ? pendingAromatic : open->second.aromatic;
                linkTo(open->second.atom, hasOrder, order, aromatic, true);
                rings.erase(open);
            }
            hasPendingOrder = false;
            pendingAromatic = false;
            continue;
        }

        // --- атом в квадратных скобках ---------------------------------------
        if (ch == '[') {
            const std::size_t close = src.find(']', i);
            if (close == std::string::npos) fail("Не закрыта квадратная скобка «[»");
            const std::string body = src.substr(i + 1, close - i - 1);
            const RawAtom atom = parseBracketAtom(body, i);
            i = close + 1;
            pushAtom(atom);
            continue;
        }

        // --- атом органического подмножества ---------------------------------
        const std::string two = src.substr(i, 2);
        if (two.size() == 2 && organicSubset().count(two) != 0) {
            i += 2;
            RawAtom a; a.el = two;
            pushAtom(a);
            continue;
        }
        if (organicSubset().count(std::string(1, ch)) != 0) {
            i += 1;
            RawAtom a; a.el = std::string(1, ch);
            pushAtom(a);
            continue;
        }
        if (isAromaticSymbol(ch)) {
            i += 1;
            RawAtom a; a.el = std::string(1, upper(ch)); a.aromatic = true;
            pushAtom(a);
            continue;
        }
        if (ch == '*') {
            i += 1;
            RawAtom a; a.el = "C";
            pushAtom(a);
            continue;
        }

        // Символ может быть многобайтовым: если забыли переключить раскладку,
        // сюда попадает кириллица. Вырезаем его целиком, иначе в сообщении
        // окажется обрубок UTF-8 и пользователь увидит «кракозябру».
        {
            std::size_t width = 1;
            const unsigned char lead = static_cast<unsigned char>(ch);
            if ((lead & 0xE0) == 0xC0) width = 2;
            else if ((lead & 0xF0) == 0xE0) width = 3;
            else if ((lead & 0xF8) == 0xF0) width = 4;
            width = std::min(width, src.size() - i);
            fail("Непонятный символ «" + src.substr(i, width) + "» в позиции "
                 + std::to_string(i + 1));
        }
    }

    if (!branchStack.empty()) fail("Не закрыта скобка «(»");
    if (!rings.empty()) {
        const int label = rings.front().first;
        throw SmilesError("Цикл с номером " + std::to_string(label) + " открыт, но не замкнут",
                          src.size());
    }
    if (raw.empty()) throw SmilesError("Не найдено ни одного атома", 0);

    return buildMolecule(raw, bonds);
}

std::vector<std::vector<int>> neighborLists(const Molecule& mol) {
    std::vector<std::vector<int>> nb(mol.atoms.size());
    for (const Bond& b : mol.bonds) {
        nb[b.a].push_back(b.b);
        nb[b.b].push_back(b.a);
    }
    return nb;
}

double bondOrderSum(const Molecule& mol, int id) {
    double s = 0;
    for (const Bond& b : mol.bonds) {
        if (b.a == id || b.b == id) s += b.order;
    }
    return s;
}

const Bond* findBond(const Molecule& mol, int a, int b) {
    for (const Bond& x : mol.bonds) {
        if ((x.a == a && x.b == b) || (x.a == b && x.b == a)) return &x;
    }
    return nullptr;
}

Bond* findBond(Molecule& mol, int a, int b) {
    for (Bond& x : mol.bonds) {
        if ((x.a == a && x.b == b) || (x.a == b && x.b == a)) return &x;
    }
    return nullptr;
}

}  // namespace chem
