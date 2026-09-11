#include "ui/panel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>

#include "chem/formula.h"
#include "chem/periodic.h"
#include "chem/resonance.h"
#include "chem/vsepr.h"

namespace ui {

using render::Align;
using render::drawText;
using render::drawTextCentered;
using render::drawWrapped;
using render::fillRoundRect;
using render::Font;
using render::Fonts;
using render::measure;
using render::measureWrapped;
using render::Rgb;
using render::strokeRoundRect;
using render::toColor;
using render::toWide;
namespace theme = render::theme;

namespace {

// --- размеры раскладки ------------------------------------------------------
constexpr double TOP_BAR_HEIGHT = 60;
constexpr double CATALOG_WIDTH = 252;
constexpr double PANEL_WIDTH = 350;
constexpr double TOOLBAR_HEIGHT = 42;
constexpr double PADDING = 14;

/** Русское склонение по числу: 1 цикл, 2 цикла, 5 циклов. */
std::string plural(int n, const std::string& one, const std::string& few, const std::string& many) {
    const int mod100 = n % 100;
    const int mod10 = n % 10;
    if (mod100 >= 11 && mod100 <= 14) return many;
    if (mod10 == 1) return one;
    if (mod10 >= 2 && mod10 <= 4) return few;
    return many;
}

/** Заглавная первая буква с учётом кириллицы в UTF-8. */
std::string capitalize(const std::string& s) {
    if (s.empty()) return s;
    const unsigned char c0 = static_cast<unsigned char>(s[0]);
    if (c0 >= 'a' && c0 <= 'z') {
        std::string out = s;
        out[0] = static_cast<char>(c0 - 'a' + 'A');
        return out;
    }
    if (s.size() >= 2 && (c0 == 0xD0 || c0 == 0xD1)) {
        const unsigned code = ((c0 & 0x1Fu) << 6) | (static_cast<unsigned char>(s[1]) & 0x3Fu);
        unsigned upper = code;
        if (code >= 0x430 && code <= 0x44F) upper = code - 32;   // а–я
        else if (code == 0x451) upper = 0x401;                   // ё
        if (upper != code) {
            std::string out = s;
            out[0] = static_cast<char>(0xC0 | (upper >> 6));
            out[1] = static_cast<char>(0x80 | (upper & 0x3F));
            return out;
        }
    }
    return s;
}

std::string number(double value, int digits) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", digits, value);
    return buffer;
}

/** Окончание для «σ-связь / связи / связей». */
std::string bondEnding(int n) { return plural(n, "ь", "и", "ей"); }
std::string pairEnding(int n) { return plural(n, "а", "ы", ""); }

/** Подпись заряда рядом с названием атома. */
std::string chargeSuffix(double charge) {
    if (std::fabs(charge) < 0.01) return "";
    if (std::fabs(charge - std::floor(charge + 0.5)) > 1e-9) {
        return " (q = " + number(charge, 2) + ")";
    }
    const std::string sign = charge > 0 ? "+" : "−";
    const double magnitude = std::fabs(charge);
    return " " + (magnitude > 1 ? std::to_string(static_cast<int>(magnitude + 0.5)) : "") + sign;
}

/** Прямоугольная плашка-«чип» с текстом. Возвращает ширину. */
double drawChip(Gdiplus::Graphics& g, Fonts& fonts, const std::string& text,
                double x, double y, const Rgb& color) {
    const Font* font = fonts.ui(11);
    const std::wstring wide = toWide(text);
    const double w = measure(g, wide, *font).Width + 16;
    fillRoundRect(g, x, y, w, 20, 10, toColor(color, 0.14));
    strokeRoundRect(g, x, y, w, 20, 10, toColor(color, 0.35));
    drawTextCentered(g, wide, *font, toColor(color), x + w / 2, y + 10);
    return w;
}

/** Кнопка панели инструментов. */
void drawButton(Gdiplus::Graphics& g, Fonts& fonts, const std::string& text, const RectD& rect,
                bool active, bool hovered) {
    const Rgb fill = active ? theme::ACCENT : theme::TEXT_DIM;
    const double bgAlpha = active ? 0.16 : (hovered ? 0.10 : 0.0);
    if (bgAlpha > 0) fillRoundRect(g, rect.x, rect.y, rect.w, rect.h, 7, toColor(fill, bgAlpha));
    if (active) strokeRoundRect(g, rect.x, rect.y, rect.w, rect.h, 7, toColor(theme::ACCENT, 0.5));
    drawTextCentered(g, toWide(text), *fonts.ui(12, active),
                     toColor(active ? theme::ACCENT : (hovered ? theme::TEXT : theme::TEXT_DIM)),
                     rect.x + rect.w / 2, rect.y + rect.h / 2);
}

/** Ширина кнопки по её подписи. */
double buttonWidth(Gdiplus::Graphics& g, Fonts& fonts, const std::string& text) {
    return measure(g, toWide(text), *fonts.ui(12, true)).Width + 20;
}

/**
 * Одинаковые по окружению атомы показываем один раз — иначе бензол
 * даёт шесть одинаковых карточек.
 */
std::vector<const chem::AtomAnalysis*> dedupeByShape(
    const std::vector<const chem::AtomAnalysis*>& atoms, const chem::Analysis& a) {
    std::vector<std::string> seen;
    std::vector<const chem::AtomAnalysis*> result;
    for (const chem::AtomAnalysis* atom : atoms) {
        std::vector<std::string> neighbors;
        for (const chem::Bond& b : a.molecule.bonds) {
            if (b.a != atom->id && b.b != atom->id) continue;
            const int other = b.a == atom->id ? b.b : b.a;
            neighbors.push_back(a.molecule.atoms[other].el + number(b.order, 2));
        }
        std::sort(neighbors.begin(), neighbors.end());
        std::string key = atom->el + "|" + number(atom->charge, 4) + "|" + atom->axe + "|";
        for (const std::string& n : neighbors) key += n + ",";
        if (std::find(seen.begin(), seen.end(), key) != seen.end()) continue;
        seen.push_back(key);
        result.push_back(atom);
    }
    return result;
}

}  // namespace

// ---------------------------------------------------------------------------
// Раскладка
// ---------------------------------------------------------------------------

void Panels::layout(double width, double height) {
    windowRect = RectD{0, 0, width, height};
    topBarRect = RectD{0, 0, width, TOP_BAR_HEIGHT};

    // на узком окне боковые панели сжимаются, но не исчезают
    const double catalogW = std::min(CATALOG_WIDTH, std::max(150.0, width * 0.22));
    const double panelW = std::min(PANEL_WIDTH, std::max(200.0, width * 0.30));

    catalogRect = RectD{0, TOP_BAR_HEIGHT, catalogW, height - TOP_BAR_HEIGHT};
    panelRect = RectD{width - panelW, TOP_BAR_HEIGHT, panelW, height - TOP_BAR_HEIGHT};
    sceneRect = RectD{catalogW, TOP_BAR_HEIGHT,
                      std::max(120.0, width - catalogW - panelW), height - TOP_BAR_HEIGHT};

    // поле поиска начинается там, где кончается подзаголовок программы
    const double searchLeft = 330;
    searchRect = RectD{searchLeft, 14, std::max(140.0, width - searchLeft - 190), 32};
}

void Panels::addHotspot(const RectD& rect, Action action, int value, const std::string& id) {
    hotspots.push_back(Hotspot{rect, action, value, id});
}

const Hotspot* Panels::hit(double x, double y) const {
    // последними нарисованы верхние слои — их и проверяем первыми
    for (auto it = hotspots.rbegin(); it != hotspots.rend(); ++it) {
        if (it->rect.contains(x, y)) return &*it;
    }
    return nullptr;
}

void Panels::scrollCatalog(double delta) {
    const double maxScroll = std::max(0.0, catalogContent - catalogRect.h + PADDING);
    catalogScroll = std::max(0.0, std::min(maxScroll, catalogScroll - delta));
}

void Panels::scrollPanel(double delta) {
    const double maxScroll = std::max(0.0, panelContent - panelRect.h + PADDING);
    panelScroll = std::max(0.0, std::min(maxScroll, panelScroll - delta));
}

void Panels::revealActive() { revealRequested = true; }

// ---------------------------------------------------------------------------
// Отрисовка целиком
// ---------------------------------------------------------------------------

void Panels::draw(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                  const render::ViewOptions& options, double mouseX, double mouseY) {
    hotspots.clear();

    drawTopBar(g, fonts, state, mouseX, mouseY);
    drawCatalog(g, fonts, state, mouseX, mouseY);
    drawAnalysisPanel(g, fonts, state, mouseX, mouseY);
    drawHeading(g, fonts, state);
    drawLegend(g, fonts, state, options);
    drawToolbar(g, fonts, options, mouseX, mouseY);
    if (!state.suggestions.empty()) drawSuggestions(g, fonts, state, mouseX, mouseY);
    if (state.helpVisible) drawHelp(g, fonts, mouseX, mouseY);
}

// ---------------------------------------------------------------------------
// Верхняя строка: заголовок, поиск, кнопки
// ---------------------------------------------------------------------------

void Panels::drawTopBar(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                        double mouseX, double mouseY) {
    (void)state;
    Gdiplus::SolidBrush bg(toColor(theme::SURFACE));
    g.FillRectangle(&bg, 0.0f, 0.0f, static_cast<Gdiplus::REAL>(topBarRect.w),
                    static_cast<Gdiplus::REAL>(topBarRect.h));
    Gdiplus::Pen line(toColor(theme::BORDER), 1.0f);
    g.DrawLine(&line, 0.0f, static_cast<Gdiplus::REAL>(topBarRect.h),
               static_cast<Gdiplus::REAL>(topBarRect.w), static_cast<Gdiplus::REAL>(topBarRect.h));

    drawText(g, L"Строение молекул", *fonts.ui(17, true), toColor(theme::TEXT), PADDING, 10);
    drawText(g, L"валентные углы и гибридизация по теории Гиллеспи", *fonts.ui(11),
             toColor(theme::TEXT_FAINT), PADDING, 33);

    // рамка вокруг нативного поля ввода — само поле рисует Windows
    strokeRoundRect(g, searchRect.x - 6, searchRect.y - 4, searchRect.w + 12, searchRect.h + 8, 8,
                    toColor(theme::BORDER));

    const double buttonY = 14;
    double x = windowRect.w - PADDING;

    const auto topButton = [&](const std::string& text, Action action) {
        const double w = buttonWidth(g, fonts, text) + 8;
        const RectD rect{x - w, buttonY, w, 32};
        drawButton(g, fonts, text, rect, false, rect.contains(mouseX, mouseY));
        strokeRoundRect(g, rect.x, rect.y, rect.w, rect.h, 7, toColor(theme::BORDER));
        addHotspot(rect, action);
        x -= w + 8;
    };

    topButton("?", Action::ShowHelp);
    topButton("Наугад", Action::Random);
}

// ---------------------------------------------------------------------------
// Каталог
// ---------------------------------------------------------------------------

void Panels::drawCatalog(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                         double mouseX, double mouseY) {
    Gdiplus::SolidBrush bg(toColor(theme::SURFACE));
    g.FillRectangle(&bg, static_cast<Gdiplus::REAL>(catalogRect.x), static_cast<Gdiplus::REAL>(catalogRect.y),
                    static_cast<Gdiplus::REAL>(catalogRect.w), static_cast<Gdiplus::REAL>(catalogRect.h));
    Gdiplus::Pen line(toColor(theme::BORDER), 1.0f);
    g.DrawLine(&line, static_cast<Gdiplus::REAL>(catalogRect.right()), static_cast<Gdiplus::REAL>(catalogRect.y),
               static_cast<Gdiplus::REAL>(catalogRect.right()), static_cast<Gdiplus::REAL>(catalogRect.bottom()));

    Gdiplus::Region saved;
    g.GetClip(&saved);
    g.SetClip(Gdiplus::RectF(static_cast<Gdiplus::REAL>(catalogRect.x), static_cast<Gdiplus::REAL>(catalogRect.y),
                             static_cast<Gdiplus::REAL>(catalogRect.w), static_cast<Gdiplus::REAL>(catalogRect.h)),
              Gdiplus::CombineModeIntersect);

    const std::string activeId = state.entry != nullptr ? state.entry->id : std::string();
    double y = catalogRect.y + 10 - catalogScroll;
    const double x = catalogRect.x + 10;
    const double w = catalogRect.w - 20;
    activeCatalogY = -1;

    // Положение считается для всех записей без исключения — иначе прокрутка
    // к активной записи не сработает, когда та ушла далеко вниз. Рисуются
    // при этом только видимые: восемьдесят девять строк стоят копейки.
    for (const chem::Category& group : chem::categories()) {
        if (y > catalogRect.y - 30 && y < catalogRect.bottom() + 10) {
            drawText(g, toWide(group.name), *fonts.ui(10.5, true), toColor(theme::TEXT_FAINT), x, y + 4);
        }
        y += 22;

        for (const chem::DbEntry* entry : group.entries) {
            const RectD rect{x, y, w, 24};
            const bool active = !activeId.empty() && entry->id == activeId;
            if (active) activeCatalogY = y + catalogScroll - catalogRect.y;

            if (y > catalogRect.y - 30 && y < catalogRect.bottom() + 10) {
                const bool hovered = rect.contains(mouseX, mouseY);
                if (active || hovered) {
                    fillRoundRect(g, rect.x, rect.y, rect.w, rect.h, 6,
                                  toColor(active ? theme::ACCENT : theme::TEXT, active ? 0.14 : 0.06));
                }
                drawText(g, toWide(entry->name), *fonts.ui(12, active),
                         toColor(active ? theme::ACCENT : theme::TEXT), rect.x + 8, rect.y + 5);
                drawText(g, toWide(entry->formulaPretty), *fonts.ui(11),
                         toColor(theme::TEXT_FAINT), rect.right() - 8, rect.y + 6, Align::Right);
                addHotspot(rect, Action::LoadMolecule, 0, entry->id);
            }
            y += 26;
        }
        y += 8;
    }
    catalogContent = y + catalogScroll - catalogRect.y;

    g.SetClip(&saved, Gdiplus::CombineModeReplace);

    // прокрутить к активной записи, если её попросили показать
    if (revealRequested && activeCatalogY >= 0) {
        revealRequested = false;
        const double visibleTop = catalogScroll;
        const double visibleBottom = catalogScroll + catalogRect.h - 40;
        if (activeCatalogY < visibleTop || activeCatalogY > visibleBottom) {
            catalogScroll = std::max(0.0, activeCatalogY - catalogRect.h / 2);
        }
    }
}

// ---------------------------------------------------------------------------
// Заголовок сцены и легенда
// ---------------------------------------------------------------------------

void Panels::drawHeading(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state) {
    if (state.analysis == nullptr) return;
    const chem::Analysis& a = *state.analysis;

    double y = sceneRect.y + 14;
    const double x = sceneRect.x + 18;

    const std::string name = state.entry != nullptr ? state.entry->name : "Построено по SMILES";
    drawText(g, toWide(name), *fonts.ui(22, true), toColor(theme::TEXT), x, y);
    y += 30;
    drawText(g, toWide(a.formulaHtml), *fonts.ui(14), toColor(theme::ACCENT), x, y);
    y += 21;
    if (state.entry == nullptr && !state.customSmiles.empty()) {
        drawText(g, toWide(state.customSmiles), *fonts.mono(11.5), toColor(theme::TEXT_FAINT), x, y);
        y += 17;
    }

    // главный центральный атом — по нему подписывается заголовок
    const chem::AtomAnalysis* main = nullptr;
    for (const chem::AtomAnalysis& atom : a.atoms) {
        if (!atom.isCentral || atom.el == "H") continue;
        if (main == nullptr || atom.sigma + atom.lonePairs > main->sigma + main->lonePairs) main = &atom;
    }
    const std::string geometry = main != nullptr
        ? main->molecularGeom + ", " + main->hybrid + "-гибридизация центрального атома"
        : "двухатомная частица";
    drawText(g, toWide(geometry), *fonts.ui(11.5), toColor(theme::TEXT_DIM), x, y);
    y += 20;
    drawText(g, L"перетаскивание — поворот · колесо — масштаб · щелчок по атому — разбор",
             *fonts.ui(10.5), toColor(theme::TEXT_FAINT), x, y);
}

void Panels::drawLegend(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                        const render::ViewOptions& options) {
    std::vector<std::pair<Rgb, std::string>> rows;
    rows.emplace_back(theme::ANGLE_ARC, "валентный угол");
    if (state.analysis != nullptr && !state.analysis->rings.empty()) {
        rows.emplace_back(theme::ANGLE_ARC_RING, "угол задан циклом");
    }
    if (options.showLonePairs) rows.emplace_back(theme::LONE_PAIR, "неподелённая пара");
    // Передние лепестки покрашены по назначению — связь или пара, — а
    // обратные все одним цветом; без подписи три цвета читаются как три
    // разных предмета.
    if (options.showOrbitals) {
        rows.emplace_back(theme::ORBITAL, "орбиталь связи");
        rows.emplace_back(theme::LONE_PAIR, "орбиталь пары");
        rows.emplace_back(theme::ORBITAL_MINUS, "обратный лепесток, ψ < 0");
    }
    if (options.showDipole) rows.emplace_back(theme::DIPOLE, "дипольный момент");

    const Font* font = fonts.ui(11);
    double y = sceneRect.y + 14;
    for (const auto& row : rows) {
        const double w = measure(g, toWide(row.second), *font).Width + 34;
        const double x = sceneRect.right() - 16 - w;
        fillRoundRect(g, x, y, w, 22, 11, toColor(theme::SURFACE, 0.75));
        Gdiplus::SolidBrush dot(toColor(row.first));
        g.FillEllipse(&dot, static_cast<Gdiplus::REAL>(x + 9), static_cast<Gdiplus::REAL>(y + 8),
                      7.0f, 7.0f);
        drawText(g, toWide(row.second), *font, toColor(theme::TEXT_DIM), x + 22, y + 5);
        y += 26;
    }
}

// ---------------------------------------------------------------------------
// Панель инструментов
// ---------------------------------------------------------------------------

void Panels::drawToolbar(Gdiplus::Graphics& g, Fonts& fonts, const render::ViewOptions& options,
                         double mouseX, double mouseY) {
    struct Item {
        std::string text;
        Action action;
        int value;
        bool active;
        bool separatorAfter;
    };

    // все элементы панели — кнопки: подписи групп нет ни у одной, иначе
    // единственная подпись читается как ещё один переключатель
    const std::vector<Item> items = {
        {"Шары", Action::SetStyle, static_cast<int>(render::Style::BallStick),
         options.style == render::Style::BallStick, false},
        {"Объём", Action::SetStyle, static_cast<int>(render::Style::SpaceFill),
         options.style == render::Style::SpaceFill, false},
        {"Каркас", Action::SetStyle, static_cast<int>(render::Style::Wire),
         options.style == render::Style::Wire, true},

        {"НЭП", Action::ToggleOption, static_cast<int>(Option::LonePairs), options.showLonePairs, false},
        {"Орбитали", Action::ToggleOption, static_cast<int>(Option::Orbitals), options.showOrbitals, false},
        {"Как считается", Action::ToggleOption, static_cast<int>(Option::OrbitalView),
         options.orbitalView == render::OrbitalView::Computed, false},
        {"Диполь", Action::ToggleOption, static_cast<int>(Option::Dipole), options.showDipole, false},
        {"Подписи", Action::ToggleOption, static_cast<int>(Option::Labels), options.showLabels, true},

        {"Без углов", Action::SetAngleMode, static_cast<int>(render::AngleMode::None),
         options.showAngles == render::AngleMode::None, false},
        {"У атома", Action::SetAngleMode, static_cast<int>(render::AngleMode::Selected),
         options.showAngles == render::AngleMode::Selected, false},
        {"Все углы", Action::SetAngleMode, static_cast<int>(render::AngleMode::All),
         options.showAngles == render::AngleMode::All, true},

        {"Вращение", Action::ToggleOption, static_cast<int>(Option::AutoRotate), options.autoRotate, false},
        {"Сброс вида", Action::ResetView, 0, false, false},
        {"Снимок", Action::SaveImage, 0, false, false},
    };

    // ширины считаются заранее: панель переносится по строкам, чтобы на узком
    // окне не залезать на панель разбора
    std::vector<double> widths;
    for (const Item& item : items) widths.push_back(buttonWidth(g, fonts, item.text));

    const double available = std::max(160.0, sceneRect.w - 28);
    const double rowHeight = 34;

    // Переносятся ЦЕЛЫЕ группы, а не отдельные кнопки: иначе переключатели
    // одного и того же режима расползаются по разным строкам.
    std::vector<std::vector<std::size_t>> groups;
    std::vector<double> groupWidths;
    {
        std::vector<std::size_t> group;
        double groupWidth = 0;
        for (std::size_t k = 0; k < items.size(); k++) {
            group.push_back(k);
            groupWidth += widths[k] + 4 + (items[k].separatorAfter ? 10 : 0);
            if (items[k].separatorAfter || k + 1 == items.size()) {
                groups.push_back(group);
                groupWidths.push_back(groupWidth);
                group.clear();
                groupWidth = 0;
            }
        }
    }

    std::vector<std::vector<std::size_t>> rows;
    std::vector<double> rowWidths;
    std::vector<std::size_t> current;
    double currentWidth = 0;
    for (std::size_t k = 0; k < groups.size(); k++) {
        if (!current.empty() && currentWidth + groupWidths[k] > available) {
            rows.push_back(current);
            rowWidths.push_back(currentWidth);
            current.clear();
            currentWidth = 0;
        }
        current.insert(current.end(), groups[k].begin(), groups[k].end());
        currentWidth += groupWidths[k];
    }
    if (!current.empty()) { rows.push_back(current); rowWidths.push_back(currentWidth); }

    double barW = 0;
    for (double w : rowWidths) barW = std::max(barW, w);
    barW += 20;
    const double barH = rowHeight * static_cast<double>(rows.size()) + 10;
    const double barX = sceneRect.x + (sceneRect.w - barW) / 2;
    const double barY = sceneRect.bottom() - barH - 14;
    toolbarRect = RectD{barX, barY, barW, barH};

    fillRoundRect(g, barX, barY, barW, barH, 12, toColor(theme::SURFACE, 0.92));
    strokeRoundRect(g, barX, barY, barW, barH, 12, toColor(theme::BORDER));

    double rowY = barY + 5;
    for (std::size_t r = 0; r < rows.size(); r++) {
        double x = barX + (barW - rowWidths[r]) / 2;
        for (std::size_t k : rows[r]) {
            const RectD rect{x, rowY + 4, widths[k], rowHeight - 8};
            drawButton(g, fonts, items[k].text, rect, items[k].active, rect.contains(mouseX, mouseY));
            addHotspot(rect, items[k].action, items[k].value);
            x += widths[k] + 4;
            // группы разделены зазором, а не линейкой: при переносе строк
            // линейка так и норовит остаться висеть в конце строки
            if (items[k].separatorAfter) x += 10;
        }
        rowY += rowHeight;
    }
}

// ---------------------------------------------------------------------------
// Панель разбора
// ---------------------------------------------------------------------------

void Panels::drawAnalysisPanel(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                               double mouseX, double mouseY) {
    Gdiplus::SolidBrush bg(toColor(theme::SURFACE));
    g.FillRectangle(&bg, static_cast<Gdiplus::REAL>(panelRect.x), static_cast<Gdiplus::REAL>(panelRect.y),
                    static_cast<Gdiplus::REAL>(panelRect.w), static_cast<Gdiplus::REAL>(panelRect.h));
    Gdiplus::Pen line(toColor(theme::BORDER), 1.0f);
    g.DrawLine(&line, static_cast<Gdiplus::REAL>(panelRect.x), static_cast<Gdiplus::REAL>(panelRect.y),
               static_cast<Gdiplus::REAL>(panelRect.x), static_cast<Gdiplus::REAL>(panelRect.bottom()));

    Gdiplus::Region saved;
    g.GetClip(&saved);
    g.SetClip(Gdiplus::RectF(static_cast<Gdiplus::REAL>(panelRect.x), static_cast<Gdiplus::REAL>(panelRect.y),
                             static_cast<Gdiplus::REAL>(panelRect.w), static_cast<Gdiplus::REAL>(panelRect.h)),
              Gdiplus::CombineModeIntersect);

    const double x = panelRect.x + PADDING;
    const double w = panelRect.w - 2 * PADDING;
    double y = panelRect.y + PADDING - panelScroll;

    // --- сообщение об ошибке ввода ---
    if (!state.errorMessage.empty()) {
        const double h1 = measureWrapped(g, toWide(state.errorMessage), *fonts.ui(13, true), w - 20);
        const double h2 = state.errorHint.empty()
            ? 0 : measureWrapped(g, toWide(state.errorHint), *fonts.ui(11.5), w - 20) + 8;
        fillRoundRect(g, x, y, w, h1 + h2 + 20, 8, toColor(theme::BAD, 0.10));
        strokeRoundRect(g, x, y, w, h1 + h2 + 20, 8, toColor(theme::BAD, 0.35));
        drawWrapped(g, toWide(state.errorMessage), *fonts.ui(13, true), toColor(theme::BAD),
                    x + 10, y + 10, w - 20);
        if (!state.errorHint.empty()) {
            drawWrapped(g, toWide(state.errorHint), *fonts.ui(11.5), toColor(theme::TEXT_DIM),
                        x + 10, y + 10 + h1 + 8, w - 20);
        }
        y += h1 + h2 + 30;
    }

    if (state.analysis == nullptr) {
        panelContent = y + panelScroll - panelRect.y;
        g.SetClip(&saved, Gdiplus::CombineModeReplace);
        return;
    }
    const chem::Analysis& a = *state.analysis;

    // --- сообщение об изомерах ---
    if (!state.notice.empty()) {
        const double h = measureWrapped(g, toWide(state.notice), *fonts.ui(11.5), w - 20);
        fillRoundRect(g, x, y, w, h + 20, 8, toColor(theme::ACCENT, 0.08));
        strokeRoundRect(g, x, y, w, h + 20, 8, toColor(theme::ACCENT, 0.3));
        drawWrapped(g, toWide(state.notice), *fonts.ui(11.5), toColor(theme::TEXT_DIM),
                    x + 10, y + 10, w - 20);
        y += h + 30;
    }

    // --- шапка ---
    const std::string name = state.entry != nullptr ? state.entry->name : "Молекула из SMILES";
    drawText(g, toWide(name), *fonts.ui(16, true), toColor(theme::TEXT), x, y);
    y += 24;
    if (state.entry == nullptr && !state.customSmiles.empty()) {
        drawText(g, toWide(state.customSmiles), *fonts.mono(11), toColor(theme::TEXT_FAINT), x, y);
        y += 16;
    }
    drawText(g, toWide(a.formulaHtml), *fonts.ui(13), toColor(theme::ACCENT), x, y);
    y += 24;

    // чипы с краткой сводкой
    double chipX = x;
    double chipY = y;
    const auto chip = [&](const std::string& text, const Rgb& color) {
        const double cw = measure(g, toWide(text), *fonts.ui(11)).Width + 16;
        if (chipX + cw > x + w) { chipX = x; chipY += 25; }
        drawChip(g, fonts, text, chipX, chipY, color);
        chipX += cw + 6;
    };
    if (state.entry != nullptr) chip(state.entry->cat, theme::TEXT_DIM);
    chip("M = " + number(a.mass, 2) + " г/моль", theme::TEXT_DIM);
    for (const auto& item : a.hybridSummary) {
        chip(item.first + " × " + std::to_string(item.second), theme::ACCENT_VIOLET);
    }
    if (!a.rings.empty()) {
        int aromatic = 0;
        for (const chem::Ring& r : a.rings) if (r.aromatic) aromatic++;
        std::string text = std::to_string(a.rings.size()) + " "
                         + plural(static_cast<int>(a.rings.size()), "цикл", "цикла", "циклов");
        if (aromatic > 0) {
            text += ", " + std::to_string(aromatic) + " "
                  + plural(aromatic, "ароматический", "ароматических", "ароматических");
        }
        chip(text, theme::ANGLE_ARC_RING);
    }
    chip(a.polar ? "полярная" : "неполярная", a.polar ? theme::WARN : theme::ACCENT);
    y = chipY + 34;

    // --- пояснение из справочника ---
    if (state.entry != nullptr && !state.entry->note.empty()) {
        const double h = measureWrapped(g, toWide(state.entry->note), *fonts.ui(11.5), w - 22);
        fillRoundRect(g, x, y, w, h + 20, 8, toColor(theme::ACCENT, 0.06));
        Gdiplus::SolidBrush edge(toColor(theme::ACCENT, 0.55));
        g.FillRectangle(&edge, static_cast<Gdiplus::REAL>(x), static_cast<Gdiplus::REAL>(y + 4),
                        2.0f, static_cast<Gdiplus::REAL>(h + 12));
        drawWrapped(g, toWide(state.entry->note), *fonts.ui(11.5), toColor(theme::TEXT_DIM),
                    x + 12, y + 10, w - 22);
        y += h + 32;
    }

    // --- атомы ---
    std::vector<const chem::AtomAnalysis*> shown;
    for (const chem::AtomAnalysis& atom : a.atoms) {
        if (state.selectedAtom != -1 ? atom.id == state.selectedAtom
                                     : (atom.isCentral && atom.el != "H")) {
            shown.push_back(&atom);
        }
    }
    if (!shown.empty()) {
        drawText(g, state.selectedAtom != -1 ? L"Выбранный атом" : L"Центральные атомы",
                 *fonts.ui(10.5, true), toColor(theme::TEXT_FAINT), x, y);
        y += 20;
        const std::vector<const chem::AtomAnalysis*> groups =
            state.selectedAtom != -1 ? shown : dedupeByShape(shown, a);
        for (const chem::AtomAnalysis* atom : groups) {
            y += drawAtomCard(g, fonts, state, *atom, x, y, w) + 10;
        }
        if (state.selectedAtom != -1) {
            const RectD rect{x, y, 190, 28};
            drawButton(g, fonts, "← ко всем центрам", rect, false, rect.contains(mouseX, mouseY));
            strokeRoundRect(g, rect.x, rect.y, rect.w, rect.h, 7, toColor(theme::BORDER));
            addHotspot(rect, Action::SelectAtom, -1);
            y += 38;
        }
    }

    // --- углы ---
    y += drawAnglesTable(g, fonts, state, x, y, w) + 14;

    // --- полярность ---
    y += drawPolarity(g, fonts, state, x, y, w) + 14;

    // --- изомеры ---
    if (state.entry != nullptr) {
        const std::vector<const chem::DbEntry*> isomers = chem::isomersOf(*state.entry);
        if (!isomers.empty()) {
            drawText(g, L"Изомеры — та же формула, другое строение", *fonts.ui(10.5, true),
                     toColor(theme::TEXT_FAINT), x, y);
            y += 20;
            double bx = x;
            for (const chem::DbEntry* iso : isomers) {
                const double bw = measure(g, toWide(iso->name), *fonts.ui(12)).Width + 20;
                if (bx + bw > x + w) { bx = x; y += 32; }
                const RectD rect{bx, y, bw, 28};
                drawButton(g, fonts, iso->name, rect, false, rect.contains(mouseX, mouseY));
                strokeRoundRect(g, rect.x, rect.y, rect.w, rect.h, 7, toColor(theme::BORDER));
                addHotspot(rect, Action::LoadMolecule, 0, iso->id);
                bx += bw + 6;
            }
            y += 40;
        }
    }

    panelContent = y + panelScroll - panelRect.y;
    g.SetClip(&saved, Gdiplus::CombineModeReplace);
}

double Panels::drawAtomCard(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                            const chem::AtomAnalysis& atom, double x, double y, double w) {
    const chem::Analysis& a = *state.analysis;
    const chem::ElementInfo& info = chem::element(atom.el);
    const chem::GeometryInfo& geometry = chem::geometryFor(atom.sigma, atom.lonePairs);
    const bool active = state.selectedAtom == atom.id || state.hoveredAtom == atom.id;

    // строки таблицы «свойство — значение»
    std::vector<std::pair<std::string, std::string>> rows = {
        {"σ-связей", std::to_string(atom.sigma)},
        {"Неподелённых пар", std::to_string(atom.lonePairs)},
        {"Стерическое число", std::to_string(atom.steric)},
        {"Тип по Гиллеспи", atom.axe},
        {"Геометрия электронных пар", geometry.electronGeom},
        {"Форма молекулы", geometry.molecularGeom},
    };
    if (atom.hasIdealAngle) rows.emplace_back("Идеальный угол", number(atom.idealAngle, 1) + "°");
    if (atom.hasPredictedAngle) rows.emplace_back("Предсказанный угол", number(atom.predictedAngle, 1) + "°");

    // формула стерического числа
    const std::string formulaLine =
        "СЧ = " + std::to_string(atom.sigma) + " σ-связ" + bondEnding(atom.sigma)
        + " + " + std::to_string(atom.lonePairs) + " пар" + pairEnding(atom.lonePairs)
        + " = " + std::to_string(atom.steric) + " → " + atom.hybrid
        + (atom.steric >= 2 && atom.steric <= 7 ? ": " + chem::hybridExplanation(atom.steric) : "");

    // кратности связей, отличных от одинарной
    std::vector<std::string> orders;
    for (const chem::Bond& b : a.molecule.bonds) {
        if ((b.a != atom.id && b.b != atom.id) || b.order == 1) continue;
        const int other = b.a == atom.id ? b.b : b.a;
        const std::string text = atom.el + "–" + a.molecule.atoms[other].el + " = "
                               + chem::formatBondOrder(b.order);
        if (std::find(orders.begin(), orders.end(), text) == orders.end()) orders.push_back(text);
    }
    std::string orderLine;
    for (std::size_t k = 0; k < orders.size(); k++) {
        orderLine += (k > 0 ? ", " : "") + orders[k];
    }

    // высота карточки считается до отрисовки: рамка рисуется первой
    const double inner = w - 24;
    const double formulaHeight = measureWrapped(g, toWide(formulaLine), *fonts.ui(11), inner);
    const double orderHeight = orderLine.empty()
        ? 0 : measureWrapped(g, toWide("Порядок связей: " + orderLine), *fonts.ui(11), inner) + 6;
    // примеры «из учебника» — не про разбираемую молекулу, поэтому отдельной строкой;
    // у SF₆, XeF₂ и подобных единственный пример — сама эта молекула, тогда строки нет
    const std::string exampleLine = geometry.examples.empty() || geometry.examples == a.formulaHtml
        ? std::string() : "Такое же окружение атома: " + geometry.examples + ".";
    const double exampleHeight = exampleLine.empty()
        ? 0 : measureWrapped(g, toWide(exampleLine), *fonts.ui(11), inner) + 6;
    const double warnHeight = atom.warning.empty()
        ? 0 : measureWrapped(g, toWide(atom.warning), *fonts.ui(11), inner - 16) + 18;
    // подпись под названием атома переносится по словам: в вебе это делает браузер,
    // здесь высоту шапки приходится считать самим
    const std::string subLine = "атом №" + std::to_string(atom.id + 1) + " · " + geometry.hint;
    const double subWidth = w - 52 - 52;
    const double headHeight = std::max(52.0,
        30.0 + measureWrapped(g, toWide(subLine), *fonts.ui(10.5), subWidth) + 10.0);
    const double height = headHeight + static_cast<double>(rows.size()) * 19 + 10
                        + formulaHeight + orderHeight + exampleHeight + warnHeight + 16;

    fillRoundRect(g, x, y, w, height, 10, toColor(theme::SURFACE_RAISED));
    strokeRoundRect(g, x, y, w, height, 10,
                    toColor(active ? theme::ACCENT : theme::BORDER, active ? 0.6 : 1.0));
    addHotspot(RectD{x, y, w, height}, Action::SelectAtom, atom.id);

    // шапка карточки: значок элемента, название, гибридизация
    const Rgb badge = render::hexToRgb(info.color);
    fillRoundRect(g, x + 12, y + 12, 30, 30, 8, toColor(badge));
    drawTextCentered(g, toWide(atom.el), *fonts.ui(14, true),
                     toColor(render::luminance(badge) > 0.55 ? Rgb{10, 14, 26} : Rgb{255, 255, 255}),
                     x + 27, y + 27);

    drawText(g, toWide(capitalize(info.name) + chargeSuffix(atom.charge)), *fonts.ui(13, true),
             toColor(theme::TEXT), x + 52, y + 13);
    drawWrapped(g, toWide(subLine), *fonts.ui(10.5), toColor(theme::TEXT_FAINT),
                x + 52, y + 30, subWidth);
    drawText(g, toWide(atom.hybrid), *fonts.ui(14, true), toColor(theme::ACCENT_VIOLET),
             x + w - 12, y + 16, Align::Right);

    double ry = y + headHeight;
    for (const auto& row : rows) {
        drawText(g, toWide(row.first), *fonts.ui(11), toColor(theme::TEXT_DIM), x + 12, ry);
        drawText(g, toWide(row.second), *fonts.ui(11, true), toColor(theme::TEXT), x + w - 12, ry, Align::Right);
        ry += 19;
    }

    ry += 8;
    drawWrapped(g, toWide(formulaLine), *fonts.ui(11), toColor(theme::TEXT_DIM), x + 12, ry, inner);
    ry += formulaHeight;
    if (!orderLine.empty()) {
        ry += 6;
        drawWrapped(g, toWide("Порядок связей: " + orderLine), *fonts.ui(11),
                    toColor(theme::TEXT_DIM), x + 12, ry, inner);
        ry += orderHeight - 6;
    }
    if (!exampleLine.empty()) {
        ry += 6;
        drawWrapped(g, toWide(exampleLine), *fonts.ui(11), toColor(theme::TEXT_FAINT),
                    x + 12, ry, inner);
        ry += exampleHeight - 6;
    }
    if (!atom.warning.empty()) {
        ry += 8;
        const double h = warnHeight - 18;
        fillRoundRect(g, x + 12, ry, inner, h + 12, 6, toColor(theme::WARN, 0.10));
        drawWrapped(g, toWide(atom.warning), *fonts.ui(11), toColor(theme::WARN, 0.95),
                    x + 20, ry + 6, inner - 16);
    }
    return height;
}

double Panels::drawAnglesTable(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                               double x, double y, double w) {
    const chem::Analysis& a = *state.analysis;
    const double startY = y;

    std::vector<const chem::AngleRecord*> records;
    for (const chem::AngleRecord& r : a.angles) {
        if (state.selectedAtom != -1 && r.center != state.selectedAtom) continue;
        records.push_back(&r);
    }
    if (!state.showHydrogenAngles) {
        std::vector<const chem::AngleRecord*> heavyOnly;
        for (const chem::AngleRecord* r : records) {
            if (a.molecule.atoms[r->i].el != "H" && a.molecule.atoms[r->j].el != "H") heavyOnly.push_back(r);
        }
        if (!heavyOnly.empty()) records = heavyOnly;
    }

    drawText(g, state.selectedAtom != -1 ? L"Валентные углы выбранного атома" : L"Валентные углы",
             *fonts.ui(10.5, true), toColor(theme::TEXT_FAINT), x, y);
    y += 20;

    if (records.empty()) {
        y += drawWrapped(g, L"В этой частице нет ни одного атома с двумя связями — "
                            L"валентного угла не существует.",
                         *fonts.ui(11.5), toColor(theme::TEXT_DIM), x, y, w);
        return y - startY;
    }

    // одинаковые по смыслу углы (например, шесть углов H-C-H в этане)
    // сворачиваем в одну строку с указанием кратности
    std::stable_sort(records.begin(), records.end(),
                     [](const chem::AngleRecord* p, const chem::AngleRecord* q) {
                         if (p->center != q->center) return p->center < q->center;
                         return p->actual < q->actual;
                     });
    std::vector<std::pair<const chem::AngleRecord*, int>> grouped;
    std::vector<std::string> keys;
    for (const chem::AngleRecord* r : records) {
        const std::string key = a.molecule.atoms[r->i].el + "|" + std::to_string(r->center) + "|"
                              + a.molecule.atoms[r->j].el + "|" + number(r->actual, 1) + "|"
                              + (r->hasPredicted ? number(r->predicted, 1) : "-") + "|"
                              + (r->hasExperimental ? number(r->experimental, 1) : "-");
        const auto found = std::find(keys.begin(), keys.end(), key);
        if (found == keys.end()) {
            keys.push_back(key);
            grouped.emplace_back(r, 1);
        } else {
            grouped[static_cast<std::size_t>(found - keys.begin())].second++;
        }
    }

    bool hasExperiment = false;
    for (const chem::AngleRecord* r : records) if (r->hasExperimental) hasExperiment = true;

    // колонки
    const double colAngle = x;
    const double colVsepr = x + w * (hasExperiment ? 0.42 : 0.55);
    const double colModel = x + w * (hasExperiment ? 0.58 : 0.78);
    const double colExp = x + w * 0.76;
    const double colDelta = x + w;

    const Font* head = fonts.ui(10);
    drawText(g, L"Угол", *head, toColor(theme::TEXT_FAINT), colAngle, y);
    drawText(g, L"ОЭПВО", *head, toColor(theme::TEXT_FAINT), colVsepr, y, Align::Right);
    drawText(g, L"модель", *head, toColor(theme::TEXT_FAINT), colModel, y, Align::Right);
    if (hasExperiment) {
        drawText(g, L"опыт", *head, toColor(theme::TEXT_FAINT), colExp, y, Align::Right);
        drawText(g, L"Δ", *head, toColor(theme::TEXT_FAINT), colDelta, y, Align::Right);
    }
    y += 16;
    Gdiplus::Pen rule(toColor(theme::BORDER), 1.0f);
    g.DrawLine(&rule, static_cast<Gdiplus::REAL>(x), static_cast<Gdiplus::REAL>(y),
               static_cast<Gdiplus::REAL>(x + w), static_cast<Gdiplus::REAL>(y));
    y += 4;

    const Font* cell = fonts.ui(11);
    const Font* cellMono = fonts.mono(11);
    for (const auto& item : grouped) {
        const chem::AngleRecord& r = *item.first;
        std::string label = a.molecule.atoms[r.i].el + "–" + a.molecule.atoms[r.center].el
                          + std::to_string(r.center + 1) + "–" + a.molecule.atoms[r.j].el;
        if (item.second > 1) label += " ×" + std::to_string(item.second);

        addHotspot(RectD{x, y - 2, w, 18}, Action::SelectAtom, r.center);
        drawText(g, toWide(label), *cell, toColor(r.inRing ? theme::ANGLE_ARC_RING : theme::TEXT),
                 colAngle, y);
        drawText(g, toWide(r.hasPredicted ? number(r.predicted, 1) + "°" : "—"), *cellMono,
                 toColor(theme::TEXT_DIM), colVsepr, y, Align::Right);
        drawText(g, toWide(number(r.actual, 1) + "°"), *cellMono, toColor(theme::TEXT),
                 colModel, y, Align::Right);
        if (r.hasExperimental) {
            const double delta = r.actual - r.experimental;
            const Rgb color = std::fabs(delta) < 3 ? theme::GOOD
                            : (std::fabs(delta) < 8 ? theme::WARN : theme::BAD);
            drawText(g, toWide(number(r.experimental, 1) + "°"), *cellMono, toColor(theme::TEXT_DIM),
                     colExp, y, Align::Right);
            drawText(g, toWide((delta >= 0 ? "+" : "") + number(delta, 1)), *cellMono, toColor(color),
                     colDelta, y, Align::Right);
        } else if (hasExperiment) {
            drawText(g, L"—", *cellMono, toColor(theme::TEXT_FAINT), colExp, y, Align::Right);
            drawText(g, L"—", *cellMono, toColor(theme::TEXT_FAINT), colDelta, y, Align::Right);
        }
        y += 18;
    }

    y += 6;
    const std::string note = std::string("«ОЭПВО» — предсказание теории для изолированного атома, ")
        + "«модель» — угол, измеренный по построенной трёхмерной структуре"
        + (hasExperiment ? ", «опыт» — справочные данные." : ".");
    y += drawWrapped(g, toWide(note), *fonts.ui(10.5), toColor(theme::TEXT_FAINT), x, y, w);
    return y - startY;
}

double Panels::drawPolarity(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                            double x, double y, double w) {
    const chem::Analysis& a = *state.analysis;
    const double startY = y;

    drawText(g, L"Полярность", *fonts.ui(10.5, true), toColor(theme::TEXT_FAINT), x, y);
    y += 20;

    const std::string text = a.polar
        ? "Диполи связей и неподелённых пар не компенсируют друг друга — молекула полярна."
        : "Диполи связей расположены симметрично и полностью компенсируются — молекула неполярна.";
    y += drawWrapped(g, toWide(text), *fonts.ui(11.5), toColor(theme::TEXT_DIM), x, y, w) + 10;

    // Расчёт надёжно отвечает на вопрос «полярна или нет» (это следствие
    // симметрии), но не даёт значения в дебаях. Поэтому шкала — безразмерная,
    // а число приводится по справочнику.
    const double level = a.dipoleValue < 0.06 ? 0 : std::min(1.0, a.dipoleValue / 3.2);
    const std::string wording = a.dipoleValue < 0.06 ? "диполя нет"
        : (a.dipoleValue < 0.9 ? "слабо полярная"
        : (a.dipoleValue < 2.0 ? "умеренно полярная" : "сильно полярная"));

    const double labelW = 74;
    const double valueW = 118;
    const double barW = std::max(40.0, w - labelW - valueW - 16);

    const auto bar = [&](const std::string& caption, double fraction, const std::string& value,
                         const Rgb& color) {
        drawText(g, toWide(caption), *fonts.ui(11), toColor(theme::TEXT_FAINT), x, y + 1);
        fillRoundRect(g, x + labelW, y + 4, barW, 8, 4, toColor(theme::BORDER, 0.8));
        if (fraction > 0) {
            fillRoundRect(g, x + labelW, y + 4, std::max(4.0, barW * fraction), 8, 4, toColor(color));
        }
        drawText(g, toWide(value), *fonts.ui(11), toColor(theme::TEXT_DIM),
                 x + labelW + barW + 8, y + 1);
        y += 22;
    };

    bar("расчёт", level, wording, theme::ACCENT);
    if (state.entry != nullptr && state.entry->hasDipole) {
        bar("эксперимент", std::min(1.0, state.entry->dipole / 4.0),
            "μ = " + number(state.entry->dipole, 2) + " Д", theme::WARN);
    }

    y += 4;
    const std::string explain =
        "Программа складывает векторы диполей связей (по разности электроотрицательностей) "
        "и вкладов неподелённых пар. Такой расчёт надёжно отвечает на вопрос «полярна молекула "
        "или нет» — это следствие симметрии, — но не даёт значения в дебаях: точный расчёт "
        "требует квантовохимических методов.";
    y += drawWrapped(g, toWide(explain), *fonts.ui(10.5), toColor(theme::TEXT_FAINT), x, y, w);
    return y - startY;
}

// ---------------------------------------------------------------------------
// Подсказки поиска и справка
// ---------------------------------------------------------------------------

void Panels::drawSuggestions(Gdiplus::Graphics& g, Fonts& fonts, const AppState& state,
                             double mouseX, double mouseY) {
    const double rowHeight = 26;
    const double h = rowHeight * static_cast<double>(state.suggestions.size()) + 8;
    const double x = searchRect.x - 6;
    const double y = searchRect.bottom() + 8;
    const double w = searchRect.w + 12;

    fillRoundRect(g, x, y, w, h, 8, toColor(theme::SURFACE_RAISED, 0.98));
    strokeRoundRect(g, x, y, w, h, 8, toColor(theme::BORDER));

    double ry = y + 4;
    for (std::size_t k = 0; k < state.suggestions.size(); k++) {
        const chem::DbEntry* entry = state.suggestions[k];
        const RectD rect{x + 4, ry, w - 8, rowHeight};
        const bool active = static_cast<int>(k) == state.suggestIndex || rect.contains(mouseX, mouseY);
        if (active) fillRoundRect(g, rect.x, rect.y, rect.w, rect.h, 6, toColor(theme::ACCENT, 0.12));
        drawText(g, toWide(entry->name), *fonts.ui(12), toColor(active ? theme::ACCENT : theme::TEXT),
                 rect.x + 8, rect.y + 6);
        drawText(g, toWide(entry->formulaPretty), *fonts.ui(11), toColor(theme::TEXT_FAINT),
                 rect.right() - 10, rect.y + 7, Align::Right);
        addHotspot(rect, Action::Suggestion, static_cast<int>(k));
        ry += rowHeight;
    }
}

void Panels::drawHelp(Gdiplus::Graphics& g, Fonts& fonts, double mouseX, double mouseY) {
    Gdiplus::SolidBrush shade(Gdiplus::Color(190, 4, 6, 12));
    g.FillRectangle(&shade, 0.0f, 0.0f, static_cast<Gdiplus::REAL>(windowRect.w),
                    static_cast<Gdiplus::REAL>(windowRect.h));
    addHotspot(windowRect, Action::CloseHelp);

    const double w = std::min(620.0, windowRect.w - 60);
    const double x = (windowRect.w - w) / 2;

    struct Section { const char* title; const char* body; };
    const Section sections[] = {
        {"Что вводить",
         "Название по-русски («вода», «бензол»), брутто-формулу (H2O, C2H6O) или строение "
         "в виде SMILES (CCO, c1ccccc1). Формула проверяется раньше SMILES: на запрос C2H6O "
         "программа показывает оба изомера, а не молча разбирает строку как SMILES."},
        {"Строение в виде SMILES",
         "Атомы записываются символами элементов: CCO — этанол. Кратные связи: = двойная, "
         "# тройная. Ветвление в круглых скобках: CC(=O)O — уксусная кислота. Цикл замыкается "
         "парой одинаковых цифр: C1CCCCC1 — циклогексан. Строчные буквы означают ароматический "
         "атом: c1ccccc1 — бензол. Заряд и явные водороды — в квадратных скобках: [NH4+]."},
        {"Управление сценой",
         "Перетаскивание мышью — поворот, колесо — масштаб, двойной щелчок — вернуть вид. "
         "Щелчок по атому открывает его разбор и оставляет только его валентные углы. "
         "Клавиша «/» переводит курсор в строку поиска, Esc закрывает это окно."},
        {"Орбитали",
         "Кнопка «Орбитали» показывает гибридные орбитали центрального атома: лепестки "
         "связей и лепестки неподелённых пар покрашены по-разному, маленький обратный лепесток "
         "— третьим цветом. Это посчитанные поверхности уровня волновой функции в приближении "
         "Слейтера — так орбитали рисуют в учебниках. «Как считается» переключает на точные "
         "функции атома водорода с радиальным узлом: положительная часть прячется в ядре, а "
         "отрицательная торчит наружу конусами. Считать так правильно, а смотреть — не на что; "
         "потому в учебниках и рисуют иначе."},
        {"Откуда берётся геометрия",
         "Координаты не хранятся в базе, а вычисляются: электронные группы расставляются "
         "по сфере минимизацией энергии их отталкивания, затем молекула наращивается обходом "
         "графа и разглаживается релаксацией. В базе лежит только связность в виде SMILES."},
    };

    // высота окна считается по содержимому: пустой низ выглядел бы неряшливо
    double contentHeight = 34;
    for (const Section& section : sections) {
        contentHeight += 20 + measureWrapped(g, toWide(section.body), *fonts.ui(11.5), w - 48) + 14;
    }
    const double h = std::min(windowRect.h - 60, contentHeight + 20 + 46);
    const double y = (windowRect.h - h) / 2;

    fillRoundRect(g, x, y, w, h, 14, toColor(theme::SURFACE));
    strokeRoundRect(g, x, y, w, h, 14, toColor(theme::BORDER));
    // щелчок по самому окну справки его не закрывает
    addHotspot(RectD{x, y, w, h}, Action::None);

    double ty = y + 20;
    drawText(g, L"Как пользоваться", *fonts.ui(18, true), toColor(theme::TEXT), x + 24, ty);
    ty += 34;

    for (const Section& section : sections) {
        drawText(g, toWide(section.title), *fonts.ui(12.5, true), toColor(theme::ACCENT), x + 24, ty);
        ty += 20;
        ty += drawWrapped(g, toWide(section.body), *fonts.ui(11.5), toColor(theme::TEXT_DIM),
                          x + 24, ty, w - 48) + 14;
    }

    const RectD close{x + w - 104, y + h - 46, 80, 30};
    drawButton(g, fonts, "Закрыть", close, false, close.contains(mouseX, mouseY));
    strokeRoundRect(g, close.x, close.y, close.w, close.h, 7, toColor(theme::BORDER));
    addHotspot(close, Action::CloseHelp);
}

}  // namespace ui
