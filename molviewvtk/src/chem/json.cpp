#include "chem/json.h"

#include <cstdlib>
#include <fstream>
#include <sstream>

namespace chem {

namespace {
const Json NULL_VALUE;
}

namespace {
std::string programDirectory;

bool readable(const std::string& path) {
    std::ifstream probe(path, std::ios::binary);
    return static_cast<bool>(probe);
}
}  // namespace

void setProgramDirectory(const std::string& path) { programDirectory = path; }

std::string findDataFile(const std::string& name) {
    // сначала рядом с программой и на уровень выше (сборка кладёт .exe в build)
    if (!programDirectory.empty()) {
        std::string prefix = programDirectory + "/";
        for (int up = 0; up < 3; up++) {
            const std::string candidate = prefix + "data/" + name;
            if (readable(candidate)) return candidate;
            prefix += "../";
        }
    }
    // затем относительно текущего каталога — так удобнее запускать тесты
    std::string prefix;
    for (int up = 0; up < 4; up++) {
        const std::string candidate = prefix + "data/" + name;
        if (readable(candidate)) return candidate;
        prefix += "../";
    }
    // ничего не нашли — вернём обычный путь, чтобы сообщение об ошибке было понятным
    return "data/" + name;
}

std::size_t Json::size() const {
    if (type_ == Type::Array) return array_.size();
    if (type_ == Type::Object) return object_.size();
    return 0;
}

bool Json::has(const std::string& key) const {
    if (type_ != Type::Object) return false;
    for (const auto& field : object_) {
        if (field.first == key) return true;
    }
    return false;
}

const Json& Json::operator[](const std::string& key) const {
    if (type_ != Type::Object) return NULL_VALUE;
    for (const auto& field : object_) {
        if (field.first == key) return field.second;
    }
    return NULL_VALUE;
}

const Json& Json::operator[](std::size_t index) const {
    if (type_ != Type::Array || index >= array_.size()) return NULL_VALUE;
    return array_[index];
}

// ---------------------------------------------------------------------------
// Разбор
// ---------------------------------------------------------------------------

class JsonParser {
public:
    explicit JsonParser(const std::string& text) : src(text) {}

    Json run() {
        skipSpace();
        Json value = parseValue();
        skipSpace();
        if (pos != src.size()) fail("лишние символы после конца документа");
        return value;
    }

private:
    const std::string& src;
    std::size_t pos = 0;

    [[noreturn]] void fail(const std::string& message) const {
        throw JsonError("JSON, позиция " + std::to_string(pos + 1) + ": " + message, pos);
    }

    char peek() const { return pos < src.size() ? src[pos] : '\0'; }
    bool eof() const { return pos >= src.size(); }

    void skipSpace() {
        while (!eof()) {
            const char c = src[pos];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') pos++;
            else break;
        }
    }

    void expect(char c) {
        if (peek() != c) fail(std::string("ожидался символ «") + c + "»");
        pos++;
    }

    bool literal(const char* word) {
        const std::size_t n = std::char_traits<char>::length(word);
        if (src.compare(pos, n, word) != 0) return false;
        pos += n;
        return true;
    }

    Json parseValue() {
        if (eof()) fail("документ оборвался");
        switch (peek()) {
            case '{': return parseObject();
            case '[': return parseArray();
            case '"': {
                Json v;
                v.type_ = Json::Type::String;
                v.string_ = parseString();
                return v;
            }
            case 't': {
                if (!literal("true")) fail("непонятное значение");
                Json v; v.type_ = Json::Type::Bool; v.bool_ = true; return v;
            }
            case 'f': {
                if (!literal("false")) fail("непонятное значение");
                Json v; v.type_ = Json::Type::Bool; v.bool_ = false; return v;
            }
            case 'n': {
                if (!literal("null")) fail("непонятное значение");
                return Json();
            }
            default: return parseNumber();
        }
    }

    Json parseObject() {
        expect('{');
        Json v;
        v.type_ = Json::Type::Object;
        skipSpace();
        if (peek() == '}') { pos++; return v; }
        for (;;) {
            skipSpace();
            if (peek() != '"') fail("имя поля должно быть строкой в кавычках");
            std::string key = parseString();
            skipSpace();
            expect(':');
            skipSpace();
            v.object_.emplace_back(std::move(key), parseValue());
            skipSpace();
            if (peek() == ',') { pos++; continue; }
            if (peek() == '}') { pos++; break; }
            fail("ожидалась запятая или «}»");
        }
        return v;
    }

    Json parseArray() {
        expect('[');
        Json v;
        v.type_ = Json::Type::Array;
        skipSpace();
        if (peek() == ']') { pos++; return v; }
        for (;;) {
            skipSpace();
            v.array_.push_back(parseValue());
            skipSpace();
            if (peek() == ',') { pos++; continue; }
            if (peek() == ']') { pos++; break; }
            fail("ожидалась запятая или «]»");
        }
        return v;
    }

    Json parseNumber() {
        const std::size_t start = pos;
        if (peek() == '-' || peek() == '+') pos++;
        while (!eof() && src[pos] >= '0' && src[pos] <= '9') pos++;
        if (peek() == '.') {
            pos++;
            while (!eof() && src[pos] >= '0' && src[pos] <= '9') pos++;
        }
        if (peek() == 'e' || peek() == 'E') {
            pos++;
            if (peek() == '-' || peek() == '+') pos++;
            while (!eof() && src[pos] >= '0' && src[pos] <= '9') pos++;
        }
        if (pos == start) fail("ожидалось число");

        // strtod округляет так же, как разбор чисел в JavaScript
        const std::string token = src.substr(start, pos - start);
        Json v;
        v.type_ = Json::Type::Number;
        v.number_ = std::strtod(token.c_str(), nullptr);
        return v;
    }

    /** Строка вместе с раскрытием escape-последовательностей. Результат — UTF-8. */
    std::string parseString() {
        expect('"');
        std::string out;
        for (;;) {
            if (eof()) fail("строка не закрыта кавычкой");
            const char c = src[pos];
            if (c == '"') { pos++; break; }
            if (c != '\\') { out.push_back(src[pos++]); continue; }

            pos++;
            if (eof()) fail("строка обрывается на знаке экранирования");
            const char esc = src[pos++];
            switch (esc) {
                case '"':  out.push_back('"');  break;
                case '\\': out.push_back('\\'); break;
                case '/':  out.push_back('/');  break;
                case 'b':  out.push_back('\b'); break;
                case 'f':  out.push_back('\f'); break;
                case 'n':  out.push_back('\n'); break;
                case 'r':  out.push_back('\r'); break;
                case 't':  out.push_back('\t'); break;
                case 'u': {
                    unsigned code = readHex4();
                    // суррогатная пара: старший заместитель дополняется младшим
                    if (code >= 0xD800 && code <= 0xDBFF && src.compare(pos, 2, "\\u") == 0) {
                        const std::size_t save = pos;
                        pos += 2;
                        const unsigned low = readHex4();
                        if (low >= 0xDC00 && low <= 0xDFFF) {
                            code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
                        } else {
                            pos = save;
                        }
                    }
                    appendUtf8(out, code);
                    break;
                }
                default: fail("неизвестная escape-последовательность");
            }
        }
        return out;
    }

    unsigned readHex4() {
        unsigned value = 0;
        for (int k = 0; k < 4; k++) {
            if (eof()) fail("оборванный код символа");
            const char c = src[pos++];
            unsigned digit = 0;
            if (c >= '0' && c <= '9') digit = static_cast<unsigned>(c - '0');
            else if (c >= 'a' && c <= 'f') digit = static_cast<unsigned>(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') digit = static_cast<unsigned>(c - 'A' + 10);
            else fail("ожидалась шестнадцатеричная цифра");
            value = value * 16 + digit;
        }
        return value;
    }

    static void appendUtf8(std::string& out, unsigned code) {
        if (code < 0x80) {
            out.push_back(static_cast<char>(code));
        } else if (code < 0x800) {
            out.push_back(static_cast<char>(0xC0 | (code >> 6)));
            out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
        } else if (code < 0x10000) {
            out.push_back(static_cast<char>(0xE0 | (code >> 12)));
            out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
        } else {
            out.push_back(static_cast<char>(0xF0 | (code >> 18)));
            out.push_back(static_cast<char>(0x80 | ((code >> 12) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (code & 0x3F)));
        }
    }
};

Json Json::parse(const std::string& text) {
    // BOM в начале файла разбору мешает, поэтому снимаем его заранее
    if (text.size() >= 3 && static_cast<unsigned char>(text[0]) == 0xEF
        && static_cast<unsigned char>(text[1]) == 0xBB
        && static_cast<unsigned char>(text[2]) == 0xBF) {
        const std::string trimmed = text.substr(3);
        return JsonParser(trimmed).run();
    }
    return JsonParser(text).run();
}

Json Json::parseFile(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw JsonError("не удалось открыть файл «" + path + "»", 0);
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return parse(buffer.str());
}

}  // namespace chem
