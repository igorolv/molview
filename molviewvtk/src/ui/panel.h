#pragma once

// Каталог слева, панель разбора справа, панель инструментов снизу и окно
// справки. Всё рисуется самостоятельно через GDI+; единственный нативный
// контрол в программе — поле ввода, и оно живёт в window.cpp.
//
// Приём «непосредственного режима»: при отрисовке панели попутно
// запоминают прямоугольники всего, по чему можно щёлкнуть. Проверка
// попадания мыши потом просто просматривает этот список — не нужно
// поддерживать отдельное дерево элементов, которое всегда расходится
// с тем, что нарисовано.

#include <string>
#include <vector>

#include "chem/database.h"
#include "chem/types.h"
#include "render/draw.h"
#include "render/renderer.h"

namespace ui {

using render::RectD;

/** Что произойдёт по щелчку. */
enum class Action {
    None,
    LoadMolecule,      // id молекулы в поле id
    SelectAtom,        // номер атома в value, −1 — снять выбор
    SetStyle,          // value — render::Style
    ToggleOption,      // value — Option
    SetAngleMode,      // value — render::AngleMode
    ResetView,
    SaveImage,         // снимок сцены в PNG
    ShowHelp,
    CloseHelp,
    Random,
    Suggestion,        // value — номер подсказки
};

enum class Option { LonePairs, Orbitals, OrbitalView, Labels, Dipole, AutoRotate };

struct Hotspot {
    RectD rect;
    Action action = Action::None;
    int value = 0;
    std::string id;
};

/** Всё, что панелям нужно знать о текущем состоянии программы. */
struct AppState {
    const chem::Analysis* analysis = nullptr;
    /** Запись справочника или nullptr, если молекула построена по SMILES. */
    const chem::DbEntry* entry = nullptr;
    std::string customSmiles;
    std::string errorMessage;
    std::string errorHint;
    /** Сообщение об изомерах после поиска по брутто-формуле. */
    std::string notice;
    /** Показывать ли в таблице углы с участием водорода. */
    bool showHydrogenAngles = false;
    int selectedAtom = -1;
    int hoveredAtom = -1;
    bool helpVisible = false;
    /** Подсказки автодополнения; пусто — список скрыт. */
    std::vector<const chem::DbEntry*> suggestions;
    int suggestIndex = -1;
};

class Panels {
public:
    void layout(double width, double height);

    const RectD& scene() const { return sceneRect; }
    const RectD& catalog() const { return catalogRect; }
    const RectD& panel() const { return panelRect; }
    const RectD& search() const { return searchRect; }

    /** Отрисовка всего интерфейса, кроме сцены. */
    void draw(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
              const render::ViewOptions& options, double mouseX, double mouseY);

    const Hotspot* hit(double x, double y) const;

    void scrollCatalog(double delta);
    void scrollPanel(double delta);
    void resetPanelScroll() { panelScroll = 0; }
    /** Прокрутить каталог так, чтобы активная запись оказалась видна. */
    void revealActive();

private:
    void drawTopBar(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                    double mouseX, double mouseY);
    void drawCatalog(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                     double mouseX, double mouseY);
    void drawAnalysisPanel(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                           double mouseX, double mouseY);
    void drawToolbar(Gdiplus::Graphics& g, render::Fonts& fonts,
                     const render::ViewOptions& options, double mouseX, double mouseY);
    void drawHeading(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state);
    void drawLegend(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                    const render::ViewOptions& options);
    void drawSuggestions(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                         double mouseX, double mouseY);
    void drawHelp(Gdiplus::Graphics& g, render::Fonts& fonts, double mouseX, double mouseY);

    /** Карточка одного атома. Возвращает занятую высоту. */
    double drawAtomCard(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                        const chem::AtomAnalysis& atom, double x, double y, double w);
    double drawAnglesTable(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                           double x, double y, double w);
    double drawPolarity(Gdiplus::Graphics& g, render::Fonts& fonts, const AppState& state,
                        double x, double y, double w);

    void addHotspot(const RectD& rect, Action action, int value = 0, const std::string& id = {});

    RectD windowRect{0, 0, 0, 0};
    RectD topBarRect{0, 0, 0, 0};
    RectD catalogRect{0, 0, 0, 0};
    RectD sceneRect{0, 0, 0, 0};
    RectD panelRect{0, 0, 0, 0};
    RectD toolbarRect{0, 0, 0, 0};
    RectD searchRect{0, 0, 0, 0};

    std::vector<Hotspot> hotspots;
    double catalogScroll = 0;
    double catalogContent = 0;
    double panelScroll = 0;
    double panelContent = 0;
    /** Положение активной записи каталога — для прокрутки к ней. */
    double activeCatalogY = -1;
    bool revealRequested = false;
};

}  // namespace ui
