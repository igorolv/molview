#pragma once

// Оформление сцены. Все цвета собраны здесь, чтобы тему можно было менять целиком.
// Заголовок намеренно не знает про GDI+: это просто числа.

#include <cmath>
#include <cstdint>
#include <string>

namespace render {

struct Rgb {
    double r = 0;
    double g = 0;
    double b = 0;
};

/** Разбор записи вида «#7dd3fc» или «#abc». */
inline Rgb hexToRgb(const std::string& hex) {
    std::string h = hex;
    if (!h.empty() && h[0] == '#') h.erase(0, 1);
    if (h.size() == 3) {
        std::string full;
        for (char c : h) { full += c; full += c; }
        h = full;
    }
    const auto part = [&h](std::size_t from) {
        return static_cast<double>(std::stoi(h.substr(from, 2), nullptr, 16));
    };
    if (h.size() < 6) return Rgb{154, 160, 170};
    return Rgb{part(0), part(2), part(4)};
}

inline Rgb mix(const Rgb& a, const Rgb& b, double t) {
    return Rgb{a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t};
}

inline Rgb lighten(const Rgb& c, double t) { return mix(c, Rgb{255, 255, 255}, t); }
inline Rgb darken(const Rgb& c, double t) { return mix(c, Rgb{8, 10, 20}, t); }

/** Яркость — по ней выбирается, писать подпись тёмным или светлым. */
inline double luminance(const Rgb& c) {
    return (0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b) / 255.0;
}

namespace theme {
const Rgb BACKGROUND = {8, 11, 20};        // #080b14
const Rgb BACKGROUND_GLOW = {20, 28, 51};  // #141c33
const Rgb BACKGROUND_MID = {11, 16, 32};   // #0b1020
const Rgb FOG = {10, 14, 26};              // #0a0e1a
const Rgb ACCENT = {94, 234, 212};         // #5eead4
const Rgb ACCENT_WARM = {251, 191, 36};    // #fbbf24
const Rgb ACCENT_VIOLET = {167, 139, 250}; // #a78bfa
const Rgb LONE_PAIR = {125, 211, 252};     // #7dd3fc
const Rgb ORBITAL = {192, 132, 252};       // #c084fc
const Rgb DIPOLE = {251, 191, 36};         // #fbbf24
const Rgb ANGLE_ARC = {94, 234, 212};      // #5eead4
const Rgb ANGLE_ARC_RING = {251, 113, 133};// #fb7185
const Rgb TEXT = {230, 236, 255};          // #e6ecff
const Rgb TEXT_DIM = {139, 151, 184};      // #8b97b8
const Rgb TEXT_FAINT = {104, 114, 145};

// --- оформление окна ---
const Rgb SURFACE = {13, 17, 30};          // фон панелей
const Rgb SURFACE_RAISED = {19, 25, 42};   // карточки
const Rgb BORDER = {35, 45, 72};
const Rgb POSITIVE = {252, 165, 165};      // положительный заряд
const Rgb NEGATIVE = {147, 197, 253};      // отрицательный заряд
const Rgb GOOD = {110, 231, 183};
const Rgb WARN = {251, 191, 36};
const Rgb BAD = {248, 113, 113};
}  // namespace theme

}  // namespace render
