#pragma once

// Тонкая прослойка над GDI+: перевод строк в UTF-16, шрифты, скруглённые
// прямоугольники, эмуляция свечения. Общая и для сцены, и для интерфейса.
//
// ВАЖНО про текст. Фигуры рисует GDI+, а текст — обычный GDI, и рисует его
// БЕЗ СГЛАЖИВАНИЯ. Это выяснено пробой, а не выбрано из общих соображений:
// на тёмном фоне мелкий текст от сглаживания мылится — штрих размазывается
// по двум пикселям вполсилы, буквы слипаются в пятно. Без сглаживания каждый
// штрих попадает в пиксель целиком и полной яркостью. ClearType дал бы ту же
// резкость, но на тёмном фоне у него цветная бахрома, и раздел 8 задания его
// прямо запрещает.

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#include <windows.h>
#include <objidl.h>
#include <gdiplus.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "render/palette.h"

namespace render {

/** Прямоугольник в пикселях окна. */
struct RectD {
    double x = 0, y = 0, w = 0, h = 0;
    bool contains(double px, double py) const {
        return px >= x && px < x + w && py >= y && py < y + h;
    }
    double right() const { return x + w; }
    double bottom() const { return y + h; }
};

/** UTF-8 → UTF-16: Windows API работает с широкими строками. */
std::wstring toWide(const std::string& utf8);

/** UTF-16 → UTF-8: обратно, для строки поиска из поля ввода. */
std::string toUtf8(const std::wstring& wide);

inline Gdiplus::Color toColor(const Rgb& c, double alpha = 1.0) {
    const auto clamp255 = [](double v) {
        if (v < 0) return 0;
        if (v > 255) return 255;
        return static_cast<int>(v + 0.5);
    };
    const auto clampA = [](double v) {
        if (v < 0) return 0;
        if (v > 1) return 255;
        return static_cast<int>(v * 255 + 0.5);
    };
    return Gdiplus::Color(static_cast<BYTE>(clampA(alpha)),
                          static_cast<BYTE>(clamp255(c.r)),
                          static_cast<BYTE>(clamp255(c.g)),
                          static_cast<BYTE>(clamp255(c.b)));
}

/** Шрифт GDI. Обёртка нужна только затем, чтобы владеть HFONT. */
class Font {
public:
    explicit Font(HFONT handle) : handle_(handle) {}
    ~Font();
    Font(const Font&) = delete;
    Font& operator=(const Font&) = delete;
    HFONT handle() const { return handle_; }

private:
    HFONT handle_ = nullptr;
};

/**
 * Хранилище шрифтов. Создавать шрифт на каждый кадр слишком дорого,
 * поэтому они кэшируются по размеру и начертанию. Размер округляется
 * до целого числа пикселей: на дробном хинтинг не может выровнять
 * штрихи по пиксельной сетке, и текст снова становится мутным.
 */
class Fonts {
public:
    /** Обычный шрифт интерфейса. */
    const Font* ui(double px, bool bold = false);
    /** Моноширинный — для чисел в таблицах и подписей углов. */
    const Font* mono(double px, bool bold = false);

private:
    const Font* get(const wchar_t* family, double px, bool bold);
    std::map<std::wstring, std::unique_ptr<Font>> cache;
};

/** Ширина и высота строки в пикселях. */
Gdiplus::SizeF measure(Gdiplus::Graphics& g, const std::wstring& text, const Font& font);

enum class Align { Left, Center, Right };

/**
 * Текст с базовой точкой (x, y) — y задаёт ВЕРХ строки.
 *
 * У GDI нет прозрачности, поэтому полупрозрачный текст заранее смешивается
 * с цветом подложки. Подложка всегда известна: под подписью угла лежит тёмная
 * плашка, под символом элемента — сам шар. Для непрозрачного текста параметр
 * не играет роли.
 */
void drawText(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
              const Gdiplus::Color& color, double x, double y, Align align = Align::Left,
              const Rgb& backdrop = theme::SURFACE);

/** Текст, у которого (x, y) — центр по обеим осям. */
void drawTextCentered(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
                      const Gdiplus::Color& color, double x, double y,
                      const Rgb& backdrop = theme::SURFACE);

/** Высота абзаца, свёрстанного в колонку заданной ширины. */
double measureWrapped(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
                      double width);

/** Абзац с переносом по словам. Возвращает занятую высоту. */
double drawWrapped(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
                   const Gdiplus::Color& color, double x, double y, double width);

/** Путь скруглённого прямоугольника. */
void addRoundRect(Gdiplus::GraphicsPath& path, double x, double y, double w, double h, double r);

void fillRoundRect(Gdiplus::Graphics& g, double x, double y, double w, double h, double r,
                   const Gdiplus::Color& color);
void strokeRoundRect(Gdiplus::Graphics& g, double x, double y, double w, double h, double r,
                     const Gdiplus::Color& color, double width = 1.0);

/**
 * Эмуляция shadowBlur из Canvas: несколько обводок с растущей толщиной
 * и падающей прозрачностью. В GDI+ размытия тени нет.
 */
void glowPath(Gdiplus::Graphics& g, const Gdiplus::GraphicsPath& path, const Rgb& color,
              double baseWidth, double alpha, int layers = 3);

/** Остановка градиента. Позиция 0 — центр (или начало), 1 — край (или конец). */
struct Stop {
    double position;
    Gdiplus::Color color;
};

/**
 * Радиальный градиент — замена createRadialGradient.
 * Роль смещённого внутреннего круга играет SetCenterPoint у PathGradientBrush:
 * именно так получается блик, сдвинутый к верхнему левому краю шара.
 */
void fillRadial(Gdiplus::Graphics& g, double cx, double cy, double rx, double ry,
                double focusX, double focusY, const std::vector<Stop>& stops);

/** Продольный градиент по заданному пути — замена createLinearGradient. */
void fillPathLinear(Gdiplus::Graphics& g, const Gdiplus::GraphicsPath& path,
                    double x0, double y0, double x1, double y1, const std::vector<Stop>& stops);

}  // namespace render
