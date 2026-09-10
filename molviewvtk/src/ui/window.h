#pragma once

// Окно программы: класс окна, цикл сообщений, мышь и клавиатура.
// Здесь же живёт единственный нативный контрол — поле ввода EDIT для строки
// поиска. Писать курсор, выделение, Ctrl+V и работу с раскладкой вручную —
// яма времени без всякой пользы для проекта.

#include "render/draw.h"

#include <string>

#include "chem/analyze.h"
#include "chem/database.h"
#include "render/renderer.h"
#include "ui/panel.h"

namespace ui {

class Window {
public:
    bool create(HINSTANCE instance, int showCommand);
    /** Цикл сообщений. Возвращает код завершения программы. */
    int run();

private:
    static LRESULT CALLBACK windowProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam);
    static LRESULT CALLBACK editProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam,
                                     UINT_PTR id, DWORD_PTR data);
    LRESULT handle(UINT message, WPARAM wParam, LPARAM lParam);

    void onCreate();
    void onSize(int width, int height);
    void onPaint();
    void onMouseMove(double x, double y);
    void onLeftDown(double x, double y);
    void onLeftUp(double x, double y);
    void onWheel(double x, double y, double delta);
    void onKeyDown(WPARAM key);
    void setHelpVisible(bool visible);
    /** Реакция на щелчок по элементу интерфейса. */
    void perform(const Hotspot& spot);

    void loadEntry(const chem::DbEntry& entry);
    void loadSmiles(const std::string& smiles);
    void submitQuery();
    void updateSuggestions();
    void applyAnalysis(chem::Analysis&& value);
    void showError(const std::string& message, const std::string& hint);
    /** Текущее время в миллисекундах от запуска — для анимации. */
    double now() const;

    void ensureBackBuffer(int width, int height);
    void releaseBackBuffer();

    HWND hwnd = nullptr;
    HWND edit = nullptr;
    HFONT editFont = nullptr;
    HBRUSH editBrush = nullptr;

    // задний буфер: рисуем в DIB-секцию, затем выводим одним BitBlt
    HDC backDc = nullptr;
    HBITMAP backBitmap = nullptr;
    HBITMAP backOld = nullptr;
    int backWidth = 0;
    int backHeight = 0;

    render::Fonts fonts;
    render::MoleculeRenderer renderer;
    Panels panels;
    AppState state;

    /** Разбор текущей молекулы. Панели и рендер смотрят на него по указателю. */
    chem::Analysis analysis;
    bool hasAnalysis = false;
    std::string customSmiles;

    double mouseX = 0;
    double mouseY = 0;
    bool draggingScene = false;
    double lastFrameTime = 0;
    long long startTicks = 0;
    long long ticksPerSecond = 1;
};

}  // namespace ui
