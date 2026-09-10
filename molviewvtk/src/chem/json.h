#pragma once

// Собственный разбор JSON методом рекурсивного спуска.
//
// Внешние библиотеки в проекте запрещены, а файлы данных простые: объекты,
// массивы, строки, числа, true/false/null. Порядок полей объекта сохраняется —
// на него опирается и порядок молекул в базе, и порядок справочных углов.

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chem {

class JsonError : public std::runtime_error {
public:
    JsonError(const std::string& message, std::size_t position)
        : std::runtime_error(message), position(position) {}
    std::size_t position;
};

/**
 * Поиск файла данных. Программу запускают то из корня проекта, то из build,
 * то двойным щелчком по .exe, поэтому каталог data ищется вверх по дереву —
 * сначала рядом с самой программой, потом относительно текущего каталога.
 */
std::string findDataFile(const std::string& name);

/**
 * Каталог, рядом с которым лежит программа. Оконная версия сообщает его при
 * запуске: рабочий каталог у ярлыка может быть каким угодно, а data лежит
 * рядом с .exe. Ядру при этом не нужно знать ни про Windows, ни про пути.
 */
void setProgramDirectory(const std::string& path);

class Json {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Json() = default;

    static Json parse(const std::string& text);
    /** Читает файл целиком и разбирает. Бросает JsonError, если файла нет. */
    static Json parseFile(const std::string& path);

    Type type() const { return type_; }
    bool isNull() const { return type_ == Type::Null; }
    bool isNumber() const { return type_ == Type::Number; }
    bool isString() const { return type_ == Type::String; }
    bool isArray() const { return type_ == Type::Array; }
    bool isObject() const { return type_ == Type::Object; }

    double number(double fallback = 0.0) const { return type_ == Type::Number ? number_ : fallback; }
    bool boolean(bool fallback = false) const { return type_ == Type::Bool ? bool_ : fallback; }
    const std::string& str() const { return string_; }

    const std::vector<Json>& items() const { return array_; }
    const std::vector<std::pair<std::string, Json>>& fields() const { return object_; }

    /** Число элементов массива или полей объекта. */
    std::size_t size() const;

    bool has(const std::string& key) const;
    /** Отсутствующее поле возвращается как Null — это удобнее проверок на каждом шагу. */
    const Json& operator[](const std::string& key) const;
    const Json& operator[](std::size_t index) const;

private:
    Type type_ = Type::Null;
    bool bool_ = false;
    double number_ = 0.0;
    std::string string_;
    std::vector<Json> array_;
    std::vector<std::pair<std::string, Json>> object_;

    friend class JsonParser;
};

}  // namespace chem
