#include "render/draw.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace render {

std::wstring toWide(const std::string& utf8) {
    if (utf8.empty()) return std::wstring();
    const int size = MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(),
                                         static_cast<int>(utf8.size()), nullptr, 0);
    std::wstring out(static_cast<std::size_t>(size), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(), static_cast<int>(utf8.size()),
                        &out[0], size);
    return out;
}

std::string toUtf8(const std::wstring& wide) {
    if (wide.empty()) return std::string();
    const int size = WideCharToMultiByte(CP_UTF8, 0, wide.c_str(), static_cast<int>(wide.size()),
                                         nullptr, 0, nullptr, nullptr);
    std::string out(static_cast<std::size_t>(size), '\0');
    WideCharToMultiByte(CP_UTF8, 0, wide.c_str(), static_cast<int>(wide.size()),
                        &out[0], size, nullptr, nullptr);
    return out;
}

// ---------------------------------------------------------------------------
// Шрифты
// ---------------------------------------------------------------------------

Font::~Font() {
    if (handle_ != nullptr) DeleteObject(handle_);
}

const Font* Fonts::get(const wchar_t* family, double px, bool bold) {
    // Размер округляется до целого: на дробном хинтинг не может выровнять
    // штрихи по пиксельной сетке, и буквы снова расплываются.
    const int height = std::max(7, static_cast<int>(px + 0.5));
    const std::wstring key = std::wstring(family) + L"|" + std::to_wstring(height)
                           + (bold ? L"|b" : L"|n");
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second.get();

    const bool mono = family[0] == L'C';  // Consolas
    HFONT handle = CreateFontW(
        -height, 0, 0, 0,
        bold ? FW_SEMIBOLD : FW_NORMAL,
        FALSE, FALSE, FALSE,
        DEFAULT_CHARSET, OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS,
        // Сглаживания нет намеренно. Мелкий текст интерфейса на тёмном фоне
        // от серого сглаживания мылится: штрих размазывается по двум пикселям
        // вполсилы, и буквы слипаются. Без сглаживания каждый штрих ложится
        // в пиксель целиком и полной яркостью. ClearType дал бы ту же резкость,
        // но на тёмном фоне у него цветная бахрома (замер: разброс R/G/B до 175
        // при 48 у чистого цвета), и раздел 8 задания его прямо запрещает.
        NONANTIALIASED_QUALITY,
        static_cast<DWORD>(mono ? (FIXED_PITCH | FF_MODERN) : (DEFAULT_PITCH | FF_DONTCARE)),
        family);

    auto font = std::make_unique<Font>(handle);
    const Font* raw = font.get();
    cache.emplace(key, std::move(font));
    return raw;
}

const Font* Fonts::ui(double px, bool bold) { return get(L"Segoe UI", px, bold); }
const Font* Fonts::mono(double px, bool bold) { return get(L"Consolas", px, bold); }

// ---------------------------------------------------------------------------
// Текст средствами GDI
// ---------------------------------------------------------------------------

namespace {

/**
 * Отдельный контекст для измерений. Обмеры идут десятками на кадр, а
 * Graphics::GetHDC каждый раз синхронизирует GDI+; для измерения настоящий
 * контекст окна не нужен — достаточно совместимого с экраном.
 */
HDC measureDc() {
    static HDC dc = nullptr;
    if (dc == nullptr) {
        HDC screen = GetDC(nullptr);
        dc = CreateCompatibleDC(screen);
        ReleaseDC(nullptr, screen);
        SetBkMode(dc, TRANSPARENT);
    }
    return dc;
}

/** Захват контекста устройства: пока он занят, обращаться к GDI+ нельзя. */
class GdiText {
public:
    GdiText(Gdiplus::Graphics& g, const Font& font) : graphics(&g) {
        dc = graphics->GetHDC();
        oldFont = static_cast<HFONT>(SelectObject(dc, font.handle()));
        oldBkMode = SetBkMode(dc, TRANSPARENT);
    }
    ~GdiText() {
        SetBkMode(dc, oldBkMode);
        SelectObject(dc, oldFont);
        graphics->ReleaseHDC(dc);
    }
    GdiText(const GdiText&) = delete;
    GdiText& operator=(const GdiText&) = delete;
    HDC handle() const { return dc; }

private:
    Gdiplus::Graphics* graphics;
    HDC dc = nullptr;
    HFONT oldFont = nullptr;
    int oldBkMode = 0;
};

/** Прозрачности у GDI нет: смешиваем цвет текста с подложкой заранее. */
COLORREF flatten(const Gdiplus::Color& color, const Rgb& backdrop) {
    const Rgb front{static_cast<double>(color.GetR()),
                    static_cast<double>(color.GetG()),
                    static_cast<double>(color.GetB())};
    const double alpha = color.GetA() / 255.0;
    const Rgb result = alpha >= 0.999 ? front : mix(backdrop, front, alpha);
    const auto byte = [](double v) {
        return static_cast<int>(std::max(0.0, std::min(255.0, v)) + 0.5);
    };
    return RGB(byte(result.r), byte(result.g), byte(result.b));
}

int rounded(double v) { return static_cast<int>(std::floor(v + 0.5)); }

}  // namespace

Gdiplus::SizeF measure(Gdiplus::Graphics& g, const std::wstring& text, const Font& font) {
    (void)g;
    HDC dc = measureDc();
    HFONT oldFont = static_cast<HFONT>(SelectObject(dc, font.handle()));
    SIZE size = {0, 0};
    GetTextExtentPoint32W(dc, text.c_str(), static_cast<int>(text.size()), &size);
    SelectObject(dc, oldFont);
    return Gdiplus::SizeF(static_cast<Gdiplus::REAL>(size.cx),
                          static_cast<Gdiplus::REAL>(size.cy));
}

void drawText(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
              const Gdiplus::Color& color, double x, double y, Align align,
              const Rgb& backdrop) {
    if (text.empty() || color.GetA() == 0) return;
    double left = x;
    if (align != Align::Left) {
        const Gdiplus::SizeF size = measure(g, text, font);
        left = align == Align::Center ? x - size.Width / 2 : x - size.Width;
    }

    GdiText target(g, font);
    SetTextColor(target.handle(), flatten(color, backdrop));
    // координаты округляются до целых: иначе глиф снова окажется между пикселями
    TextOutW(target.handle(), rounded(left), rounded(y), text.c_str(), static_cast<int>(text.size()));
}

void drawTextCentered(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
                      const Gdiplus::Color& color, double x, double y, const Rgb& backdrop) {
    if (text.empty() || color.GetA() == 0) return;
    const Gdiplus::SizeF size = measure(g, text, font);
    drawText(g, text, font, color, x - size.Width / 2, y - size.Height / 2, Align::Left, backdrop);
}

namespace {

/** Ширина строки в пикселях на контексте для измерений (шрифт уже выбран). */
double widthOf(HDC dc, const std::wstring& text) {
    SIZE size = {0, 0};
    GetTextExtentPoint32W(dc, text.c_str(), static_cast<int>(text.size()), &size);
    return size.cx;
}

/**
 * Разбиение абзаца на строки.
 *
 * Свой перенос, а не DrawTextW с DT_WORDBREAK: тот на этом тексте рвёт слова
 * посередине («поля/рна»), а здесь всё под контролем — и, главное, замер
 * и отрисовка гарантированно разбивают текст одинаково.
 */
std::vector<std::wstring> wrapLines(const Font& font, const std::wstring& text, double width) {
    std::vector<std::wstring> lines;
    if (text.empty()) return lines;

    HDC dc = measureDc();
    HFONT oldFont = static_cast<HFONT>(SelectObject(dc, font.handle()));

    std::wstring line;
    std::size_t i = 0;
    while (i < text.size()) {
        // слово вместе с идущими за ним пробелами
        std::size_t end = text.find(L' ', i);
        if (end == std::wstring::npos) end = text.size();
        std::wstring word = text.substr(i, end - i);
        std::size_t next = end;
        while (next < text.size() && text[next] == L' ') next++;

        const std::wstring candidate = line.empty() ? word : line + L" " + word;
        if (!line.empty() && widthOf(dc, candidate) > width) {
            lines.push_back(line);
            line = word;
        } else {
            line = candidate;
        }

        // слово длиннее всей колонки приходится резать по буквам
        while (widthOf(dc, line) > width && line.size() > 1) {
            std::size_t fit = line.size() - 1;
            while (fit > 1 && widthOf(dc, line.substr(0, fit)) > width) fit--;
            lines.push_back(line.substr(0, fit));
            line = line.substr(fit);
        }
        i = next;
    }
    if (!line.empty()) lines.push_back(line);

    SelectObject(dc, oldFont);
    return lines;
}

/** Высота строки текста для данного шрифта. */
double lineHeight(const Font& font) {
    HDC dc = measureDc();
    HFONT oldFont = static_cast<HFONT>(SelectObject(dc, font.handle()));
    TEXTMETRICW metrics = {};
    GetTextMetricsW(dc, &metrics);
    SelectObject(dc, oldFont);
    return metrics.tmHeight;
}

}  // namespace

double measureWrapped(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
                      double width) {
    (void)g;
    if (text.empty()) return 0;
    return static_cast<double>(wrapLines(font, text, width).size()) * lineHeight(font);
}

double drawWrapped(Gdiplus::Graphics& g, const std::wstring& text, const Font& font,
                   const Gdiplus::Color& color, double x, double y, double width) {
    if (text.empty() || color.GetA() == 0) return 0;
    const std::vector<std::wstring> lines = wrapLines(font, text, width);
    const double step = lineHeight(font);

    GdiText target(g, font);
    SetTextColor(target.handle(), flatten(color, theme::SURFACE));
    double ly = y;
    for (const std::wstring& line : lines) {
        TextOutW(target.handle(), rounded(x), rounded(ly), line.c_str(),
                 static_cast<int>(line.size()));
        ly += step;
    }
    return static_cast<double>(lines.size()) * step;
}

// ---------------------------------------------------------------------------
// Фигуры (по-прежнему GDI+)
// ---------------------------------------------------------------------------

void addRoundRect(Gdiplus::GraphicsPath& path, double x, double y, double w, double h, double r) {
    const double radius = std::max(0.0, std::min(r, std::min(w, h) / 2));
    const auto f = [](double v) { return static_cast<Gdiplus::REAL>(v); };
    if (radius <= 0.01) {
        path.AddRectangle(Gdiplus::RectF(f(x), f(y), f(w), f(h)));
        return;
    }
    const double d = radius * 2;
    path.AddArc(f(x), f(y), f(d), f(d), 180, 90);
    path.AddArc(f(x + w - d), f(y), f(d), f(d), 270, 90);
    path.AddArc(f(x + w - d), f(y + h - d), f(d), f(d), 0, 90);
    path.AddArc(f(x), f(y + h - d), f(d), f(d), 90, 90);
    path.CloseFigure();
}

void fillRoundRect(Gdiplus::Graphics& g, double x, double y, double w, double h, double r,
                   const Gdiplus::Color& color) {
    Gdiplus::GraphicsPath path;
    addRoundRect(path, x, y, w, h, r);
    Gdiplus::SolidBrush brush(color);
    g.FillPath(&brush, &path);
}

void strokeRoundRect(Gdiplus::Graphics& g, double x, double y, double w, double h, double r,
                     const Gdiplus::Color& color, double width) {
    Gdiplus::GraphicsPath path;
    addRoundRect(path, x, y, w, h, r);
    Gdiplus::Pen pen(color, static_cast<Gdiplus::REAL>(width));
    g.DrawPath(&pen, &path);
}

void fillRadial(Gdiplus::Graphics& g, double cx, double cy, double rx, double ry,
                double focusX, double focusY, const std::vector<Stop>& stops) {
    if (rx <= 0.05 || ry <= 0.05 || stops.empty()) return;
    const auto f = [](double v) { return static_cast<Gdiplus::REAL>(v); };

    Gdiplus::GraphicsPath path;
    path.AddEllipse(f(cx - rx), f(cy - ry), f(rx * 2), f(ry * 2));
    Gdiplus::PathGradientBrush brush(&path);
    brush.SetCenterPoint(Gdiplus::PointF(f(focusX), f(focusY)));

    // У PathGradientBrush позиция 0 — это КРАЙ фигуры, а 1 — центр,
    // то есть порядок остановок обратный привычному по Canvas.
    std::vector<Gdiplus::Color> colors;
    std::vector<Gdiplus::REAL> positions;
    for (auto it = stops.rbegin(); it != stops.rend(); ++it) {
        colors.push_back(it->color);
        positions.push_back(f(1.0 - it->position));
    }
    if (positions.front() > 0.0001f) {
        colors.insert(colors.begin(), stops.back().color);
        positions.insert(positions.begin(), 0.0f);
    }
    if (positions.back() < 0.9999f) {
        colors.push_back(stops.front().color);
        positions.push_back(1.0f);
    }

    brush.SetInterpolationColors(colors.data(), positions.data(), static_cast<INT>(colors.size()));
    g.FillPath(&brush, &path);
}

void fillPathLinear(Gdiplus::Graphics& g, const Gdiplus::GraphicsPath& path,
                    double x0, double y0, double x1, double y1, const std::vector<Stop>& stops) {
    if (stops.size() < 2) return;
    const auto f = [](double v) { return static_cast<Gdiplus::REAL>(v); };
    // вырожденный отрезок GDI+ не переваривает — чуть раздвигаем
    if (std::abs(x1 - x0) < 0.01 && std::abs(y1 - y0) < 0.01) x1 += 0.05;

    Gdiplus::LinearGradientBrush brush(
        Gdiplus::PointF(f(x0), f(y0)), Gdiplus::PointF(f(x1), f(y1)),
        stops.front().color, stops.back().color);

    std::vector<Gdiplus::Color> colors;
    std::vector<Gdiplus::REAL> positions;
    for (const Stop& s : stops) {
        colors.push_back(s.color);
        positions.push_back(f(s.position));
    }
    positions.front() = 0.0f;
    positions.back() = 1.0f;
    brush.SetInterpolationColors(colors.data(), positions.data(), static_cast<INT>(colors.size()));
    g.FillPath(&brush, &path);
}

void glowPath(Gdiplus::Graphics& g, const Gdiplus::GraphicsPath& path, const Rgb& color,
              double baseWidth, double alpha, int layers) {
    // от широкой и почти прозрачной обводки к узкой и яркой
    for (int k = layers; k >= 1; k--) {
        const double width = baseWidth + k * 2.0;
        const double a = alpha * 0.16 / k;
        Gdiplus::Pen pen(toColor(color, a), static_cast<Gdiplus::REAL>(width));
        pen.SetLineJoin(Gdiplus::LineJoinRound);
        g.DrawPath(&pen, &path);
    }
}

}  // namespace render
