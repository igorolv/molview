#pragma once

// Мелочи, общие для обеих проверок: выравнивание колонок и вывод в UTF-8.
// Дополнять таблицы пробелами по числу БАЙТОВ нельзя — русский текст в UTF-8
// занимает два байта на букву, и колонки разъезжаются.

#include <cstdio>
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace testutil {

/** Число символов (а не байтов) в строке UTF-8. */
inline std::size_t charCount(const std::string& s) {
    std::size_t n = 0;
    for (char c : s) {
        if ((static_cast<unsigned char>(c) & 0xC0) != 0x80) n++;
    }
    return n;
}

/** Первые n символов строки. */
inline std::string charSlice(const std::string& s, std::size_t n) {
    std::size_t seen = 0;
    for (std::size_t i = 0; i < s.size(); i++) {
        if ((static_cast<unsigned char>(s[i]) & 0xC0) != 0x80) {
            if (seen == n) return s.substr(0, i);
            seen++;
        }
    }
    return s;
}

/** Дополнение справа до n символов (или обрезка). */
inline std::string pad(const std::string& s, std::size_t n) {
    const std::size_t have = charCount(s);
    if (have >= n) return charSlice(s, n);
    return s + std::string(n - have, ' ');
}

/** Дополнение слева до n символов. */
inline std::string padLeft(const std::string& s, std::size_t n) {
    const std::size_t have = charCount(s);
    if (have >= n) return s;
    return std::string(n - have, ' ') + s;
}

inline std::string repeat(const std::string& s, std::size_t n) {
    std::string out;
    for (std::size_t i = 0; i < n; i++) out += s;
    return out;
}

/** Число с фиксированным числом знаков после запятой. */
inline std::string fixed(double value, int digits) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", digits, value);
    return buffer;
}

/** Консоль Windows по умолчанию не в UTF-8, и русский текст превращается в кашу. */
inline void enableUtf8Console() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
}

/** Копилка обнаруженных расхождений. */
class Problems {
public:
    void fail(const std::string& message) { list.push_back(message); }
    std::size_t count() const { return list.size(); }
    const std::vector<std::string>& all() const { return list; }

private:
    std::vector<std::string> list;
};

}  // namespace testutil
