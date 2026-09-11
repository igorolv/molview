#include "ui/window.h"

#include <commctrl.h>
#include <windowsx.h>

#include <algorithm>
#include <cstdlib>

#include "chem/smiles.h"
#include "ui/panel.h"

namespace ui {

using render::toUtf8;
using render::toWide;
namespace theme = render::theme;

namespace {

const wchar_t* const WINDOW_CLASS = L"MolviewMainWindow";
constexpr UINT_PTR TIMER_ID = 1;
constexpr UINT_PTR EDIT_ID = 100;
/** Примерно 60 кадров в секунду. */
constexpr UINT FRAME_INTERVAL = 16;

/** Случайная запись базы для кнопки «Наугад». */
const chem::DbEntry& randomEntry() {
    const std::vector<chem::DbEntry>& db = chem::database();
    return db[static_cast<std::size_t>(std::rand()) % db.size()];
}

}  // namespace

// ---------------------------------------------------------------------------
// Создание окна
// ---------------------------------------------------------------------------

bool Window::create(HINSTANCE instance, int showCommand) {
    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(wc);
    wc.style = CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS;
    wc.lpfnWndProc = &Window::windowProc;
    wc.hInstance = instance;
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hbrBackground = nullptr;  // фон рисуем сами, иначе будет мерцание
    wc.lpszClassName = WINDOW_CLASS;
    if (RegisterClassExW(&wc) == 0) return false;

    hwnd = CreateWindowExW(
        0, WINDOW_CLASS, L"Строение молекул — валентные углы и гибридизация",
        // WS_CLIPCHILDREN обязателен: без него BitBlt заднего буфера затирает
        // поле ввода, и текст в нём то виден, то пропадает до следующей
        // перерисовки самого контрола
        WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN, CW_USEDEFAULT, CW_USEDEFAULT, 1280, 820,
        nullptr, nullptr, instance, this);
    if (hwnd == nullptr) return false;

    ShowWindow(hwnd, showCommand);
    UpdateWindow(hwnd);
    return true;
}

int Window::run() {
    MSG message;
    while (GetMessageW(&message, nullptr, 0, 0) > 0) {
        TranslateMessage(&message);
        DispatchMessageW(&message);
    }
    return static_cast<int>(message.wParam);
}

LRESULT CALLBACK Window::windowProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
    Window* self = nullptr;
    if (message == WM_NCCREATE) {
        const CREATESTRUCTW* create = reinterpret_cast<CREATESTRUCTW*>(lParam);
        self = static_cast<Window*>(create->lpCreateParams);
        self->hwnd = hwnd;
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    } else {
        self = reinterpret_cast<Window*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    }
    if (self != nullptr) return self->handle(message, wParam, lParam);
    return DefWindowProcW(hwnd, message, wParam, lParam);
}

/**
 * Поле ввода перехватывает Enter и стрелки: без этого Windows съедает их сама,
 * и выбрать подсказку с клавиатуры было бы нельзя.
 */
LRESULT CALLBACK Window::editProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam,
                                  UINT_PTR id, DWORD_PTR data) {
    (void)id;
    Window* self = reinterpret_cast<Window*>(data);
    if (message == WM_KEYDOWN && self != nullptr) {
        if (wParam == VK_RETURN) {
            if (self->state.suggestIndex >= 0
                && self->state.suggestIndex < static_cast<int>(self->state.suggestions.size())) {
                self->loadEntry(*self->state.suggestions[static_cast<std::size_t>(self->state.suggestIndex)]);
                self->state.suggestions.clear();
                self->state.suggestIndex = -1;
            } else {
                self->submitQuery();
            }
            return 0;
        }
        if (wParam == VK_DOWN || wParam == VK_UP) {
            const int delta = wParam == VK_DOWN ? 1 : -1;
            const int count = static_cast<int>(self->state.suggestions.size());
            self->state.suggestIndex = std::max(-1, std::min(count - 1, self->state.suggestIndex + delta));
            return 0;
        }
        if (wParam == VK_ESCAPE) {
            self->state.suggestions.clear();
            self->state.suggestIndex = -1;
            SetFocus(self->hwnd);
            return 0;
        }
    }
    return DefSubclassProc(hwnd, message, wParam, lParam);
}

// ---------------------------------------------------------------------------
// Обработка сообщений
// ---------------------------------------------------------------------------

LRESULT Window::handle(UINT message, WPARAM wParam, LPARAM lParam) {
    switch (message) {
        case WM_CREATE:
            onCreate();
            return 0;

        case WM_SIZE:
            onSize(LOWORD(lParam), HIWORD(lParam));
            return 0;

        case WM_ERASEBKGND:
            return 1;  // фон целиком рисует WM_PAINT

        case WM_PAINT:
            onPaint();
            return 0;

        case WM_TIMER:
            if (wParam == TIMER_ID) {
                const double time = now();
                renderer.tick(time - lastFrameTime);
                lastFrameTime = time;
                InvalidateRect(hwnd, nullptr, FALSE);
            }
            return 0;

        case WM_MOUSEMOVE:
            onMouseMove(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
            return 0;

        case WM_LBUTTONDOWN:
            onLeftDown(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
            return 0;

        case WM_LBUTTONUP:
            onLeftUp(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam));
            return 0;

        case WM_LBUTTONDBLCLK: {
            const double x = GET_X_LPARAM(lParam);
            const double y = GET_Y_LPARAM(lParam);
            if (panels.scene().contains(x, y)) renderer.resetView();
            return 0;
        }

        case WM_MOUSEWHEEL: {
            POINT point{GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)};
            ScreenToClient(hwnd, &point);
            onWheel(point.x, point.y, GET_WHEEL_DELTA_WPARAM(wParam));
            return 0;
        }

        case WM_KEYDOWN:
            onKeyDown(wParam);
            return 0;

        case WM_COMMAND:
            if (LOWORD(wParam) == EDIT_ID && HIWORD(wParam) == EN_CHANGE) updateSuggestions();
            return 0;

        case WM_CTLCOLOREDIT: {
            // тёмное оформление поля ввода: светлый текст на фоне панели
            HDC dc = reinterpret_cast<HDC>(wParam);
            SetTextColor(dc, RGB(static_cast<BYTE>(theme::TEXT.r), static_cast<BYTE>(theme::TEXT.g),
                                 static_cast<BYTE>(theme::TEXT.b)));
            SetBkColor(dc, RGB(static_cast<BYTE>(theme::SURFACE_RAISED.r),
                               static_cast<BYTE>(theme::SURFACE_RAISED.g),
                               static_cast<BYTE>(theme::SURFACE_RAISED.b)));
            return reinterpret_cast<LRESULT>(editBrush);
        }

        case WM_DESTROY:
            KillTimer(hwnd, TIMER_ID);
            releaseBackBuffer();
            if (editFont != nullptr) DeleteObject(editFont);
            if (editBrush != nullptr) DeleteObject(editBrush);
            PostQuitMessage(0);
            return 0;

        default:
            break;
    }
    return DefWindowProcW(hwnd, message, wParam, lParam);
}

void Window::onCreate() {
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    ticksPerSecond = frequency.QuadPart;
    startTicks = counter.QuadPart;

    editBrush = CreateSolidBrush(RGB(static_cast<BYTE>(theme::SURFACE_RAISED.r),
                                     static_cast<BYTE>(theme::SURFACE_RAISED.g),
                                     static_cast<BYTE>(theme::SURFACE_RAISED.b)));
    editFont = CreateFontW(-17, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET,
                           OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                           DEFAULT_PITCH | FF_DONTCARE, L"Segoe UI");

    edit = CreateWindowExW(0, L"EDIT", L"",
                           WS_CHILD | WS_VISIBLE | ES_AUTOHSCROLL,
                           0, 0, 10, 10, hwnd, reinterpret_cast<HMENU>(EDIT_ID),
                           reinterpret_cast<HINSTANCE>(GetWindowLongPtrW(hwnd, GWLP_HINSTANCE)),
                           nullptr);
    SendMessageW(edit, WM_SETFONT, reinterpret_cast<WPARAM>(editFont), TRUE);
    SetWindowSubclass(edit, &Window::editProc, 1, reinterpret_cast<DWORD_PTR>(this));

    SetTimer(hwnd, TIMER_ID, FRAME_INTERVAL, nullptr);
    lastFrameTime = now();

    // стартовая молекула — вода, классический пример AX₂E₂
    const chem::DbEntry* water = chem::entryById("water");
    if (water != nullptr) loadEntry(*water);
}

void Window::onSize(int width, int height) {
    if (width <= 0 || height <= 0) return;
    panels.layout(width, height);
    renderer.setViewport(panels.scene());
    ensureBackBuffer(width, height);

    const RectD& search = panels.search();
    if (edit != nullptr) {
        MoveWindow(edit, static_cast<int>(search.x), static_cast<int>(search.y),
                   static_cast<int>(search.w), static_cast<int>(search.h), TRUE);
    }
}

// ---------------------------------------------------------------------------
// Отрисовка
// ---------------------------------------------------------------------------

void Window::ensureBackBuffer(int width, int height) {
    if (backDc != nullptr && width == backWidth && height == backHeight) return;
    releaseBackBuffer();

    HDC windowDc = GetDC(hwnd);
    backDc = CreateCompatibleDC(windowDc);

    BITMAPINFO info = {};
    info.bmiHeader.biSize = sizeof(info.bmiHeader);
    info.bmiHeader.biWidth = width;
    info.bmiHeader.biHeight = -height;  // сверху вниз
    info.bmiHeader.biPlanes = 1;
    info.bmiHeader.biBitCount = 32;
    info.bmiHeader.biCompression = BI_RGB;

    void* bits = nullptr;
    backBitmap = CreateDIBSection(windowDc, &info, DIB_RGB_COLORS, &bits, nullptr, 0);
    backOld = static_cast<HBITMAP>(SelectObject(backDc, backBitmap));
    ReleaseDC(hwnd, windowDc);

    backWidth = width;
    backHeight = height;
}

void Window::releaseBackBuffer() {
    if (backDc == nullptr) return;
    SelectObject(backDc, backOld);
    DeleteObject(backBitmap);
    DeleteDC(backDc);
    backDc = nullptr;
    backBitmap = nullptr;
    backOld = nullptr;
}

void Window::onPaint() {
    PAINTSTRUCT ps;
    HDC dc = BeginPaint(hwnd, &ps);
    if (backDc != nullptr) {
        Gdiplus::Graphics g(backDc);
        g.SetSmoothingMode(Gdiplus::SmoothingModeAntiAlias);
        // На тёмном фоне ClearType даёт цветную бахрому вокруг букв,
        // поэтому именно AntiAliasGridFit, а не ClearTypeGridFit.
        g.SetTextRenderingHint(Gdiplus::TextRenderingHintAntiAliasGridFit);
        g.SetInterpolationMode(Gdiplus::InterpolationModeHighQualityBicubic);

        // общий фон окна — под панелями и сценой
        Gdiplus::SolidBrush background(render::toColor(theme::BACKGROUND));
        g.FillRectangle(&background, 0, 0, backWidth, backHeight);

        const double time = now();

        // Кадр собирается в три слоя: фон и сетка от GDI+, сцена от VTK
        // поверх них, а сверху — то, что VTK пока не рисует. Между слоями
        // нужен Flush: GDI+ и обычный GDI пишут в один растр, и порядок
        // операций должен быть тем, в каком они выписаны.
        renderer.renderBackground(g, time);
        g.Flush(Gdiplus::FlushIntentionSync);
        scene.draw(backDc, renderer);
        renderer.renderOverlays(g, fonts, time);

        state.analysis = hasAnalysis ? &analysis : nullptr;
        state.selectedAtom = renderer.selected();
        state.hoveredAtom = renderer.hovered();
        panels.draw(g, fonts, state, renderer.options(), mouseX, mouseY);

        // готовый кадр выводится одним переносом — так нет мерцания
        BitBlt(dc, 0, 0, backWidth, backHeight, backDc, 0, 0, SRCCOPY);
    }
    EndPaint(hwnd, &ps);
}

double Window::now() const {
    LARGE_INTEGER counter;
    QueryPerformanceCounter(&counter);
    return static_cast<double>(counter.QuadPart - startTicks) * 1000.0
         / static_cast<double>(ticksPerSecond);
}

// ---------------------------------------------------------------------------
// Мышь и клавиатура
// ---------------------------------------------------------------------------

void Window::onMouseMove(double x, double y) {
    mouseX = x;
    mouseY = y;

    if (draggingScene) {
        renderer.dragTo(x, y);
        return;
    }
    if (panels.scene().contains(x, y) && !state.helpVisible) {
        renderer.setHovered(renderer.hitTest(x, y));
        SetCursor(LoadCursor(nullptr, renderer.hovered() != -1 ? IDC_HAND : IDC_ARROW));
    } else {
        renderer.setHovered(-1);
    }
}

void Window::onLeftDown(double x, double y) {
    mouseX = x;
    mouseY = y;

    const Hotspot* spot = panels.hit(x, y);
    if (spot != nullptr) {
        perform(*spot);
        return;
    }
    if (state.helpVisible) return;

    if (panels.scene().contains(x, y)) {
        draggingScene = true;
        renderer.beginDrag(x, y);
        SetCapture(hwnd);
    }
}

void Window::onLeftUp(double x, double y) {
    if (!draggingScene) return;
    ReleaseCapture();
    draggingScene = false;

    if (renderer.endDrag()) {
        // щелчок без перетаскивания — выбор атома
        const int id = renderer.hitTest(x, y);
        renderer.setHovered(id);
        renderer.setSelected(id);
        panels.resetPanelScroll();
    }
}

void Window::onWheel(double x, double y, double delta) {
    if (state.helpVisible) return;
    if (panels.catalog().contains(x, y)) {
        panels.scrollCatalog(delta * 0.5);
    } else if (panels.panel().contains(x, y)) {
        panels.scrollPanel(delta * 0.5);
    } else if (panels.scene().contains(x, y)) {
        // знак противоположен браузерному deltaY, отсюда минус
        renderer.wheel(-delta);
    }
}

/**
 * Окно справки затемняет всё вокруг, но нативный контрол рисуется поверх
 * любой нашей графики. Поэтому на время справки поле ввода прячется.
 */
void Window::setHelpVisible(bool visible) {
    state.helpVisible = visible;
    if (edit != nullptr) ShowWindow(edit, visible ? SW_HIDE : SW_SHOW);
}

void Window::onKeyDown(WPARAM key) {
    if (key == VK_ESCAPE) {
        if (state.helpVisible) setHelpVisible(false);
        else if (renderer.selected() != -1) renderer.setSelected(-1);
        return;
    }
    // «/» переводит курсор в строку поиска — как в веб-версии
    if (key == VK_OEM_2 || key == VK_DIVIDE) {
        SetFocus(edit);
        SendMessageW(edit, EM_SETSEL, 0, -1);
        return;
    }
    if (key == VK_F1) setHelpVisible(true);
}

void Window::perform(const Hotspot& spot) {
    switch (spot.action) {
        case Action::None:
            break;

        case Action::LoadMolecule: {
            const chem::DbEntry* entry = chem::entryById(spot.id);
            if (entry != nullptr) loadEntry(*entry);
            break;
        }

        case Action::SelectAtom:
            renderer.setSelected(spot.value);
            panels.resetPanelScroll();
            break;

        case Action::SetStyle:
            renderer.options().style = static_cast<render::Style>(spot.value);
            break;

        case Action::ToggleOption: {
            render::ViewOptions& options = renderer.options();
            switch (static_cast<Option>(spot.value)) {
                case Option::LonePairs: options.showLonePairs = !options.showLonePairs; break;
                case Option::Orbitals:  options.showOrbitals = !options.showOrbitals; break;
                // Режим «как считается» без орбиталей на экране показать нечего,
                // поэтому его включение заодно включает и их.
                case Option::OrbitalView:
                    options.orbitalView = options.orbitalView == render::OrbitalView::Textbook
                                              ? render::OrbitalView::Computed
                                              : render::OrbitalView::Textbook;
                    if (options.orbitalView == render::OrbitalView::Computed) options.showOrbitals = true;
                    break;
                case Option::Labels:    options.showLabels = !options.showLabels; break;
                case Option::Dipole:    options.showDipole = !options.showDipole; break;
                case Option::AutoRotate: options.autoRotate = !options.autoRotate; break;
            }
            break;
        }

        case Action::SetAngleMode:
            renderer.options().showAngles = static_cast<render::AngleMode>(spot.value);
            break;

        case Action::ResetView:
            renderer.resetView();
            break;

        case Action::SaveImage: {
            // Действие редкое и явное, поэтому окно с подтверждением уместно:
            // иначе непонятно, произошло ли вообще что-нибудь и куда смотреть.
            const std::wstring path = scene.saveImage(renderer);
            if (path.empty()) {
                MessageBoxW(hwnd, L"Не удалось сохранить снимок.", L"Строение молекул",
                            MB_ICONWARNING);
            } else {
                MessageBoxW(hwnd, (L"Снимок сохранён:\n" + path).c_str(), L"Строение молекул",
                            MB_ICONINFORMATION);
            }
            break;
        }

        case Action::ShowHelp:
            setHelpVisible(true);
            break;

        case Action::CloseHelp:
            setHelpVisible(false);
            break;

        case Action::Random:
            loadEntry(randomEntry());
            break;

        case Action::Suggestion:
            if (spot.value >= 0 && spot.value < static_cast<int>(state.suggestions.size())) {
                loadEntry(*state.suggestions[static_cast<std::size_t>(spot.value)]);
            }
            break;
    }
}

// ---------------------------------------------------------------------------
// Загрузка молекулы
// ---------------------------------------------------------------------------

void Window::applyAnalysis(chem::Analysis&& value) {
    analysis = std::move(value);
    hasAnalysis = true;

    // у маленьких частиц углы с водородом интересны, у больших только мешают
    int heavy = 0;
    for (const chem::Atom& a : analysis.molecule.atoms) {
        if (a.el != "H") heavy++;
    }
    state.showHydrogenAngles = heavy <= 4;
    state.analysis = &analysis;
    state.errorMessage.clear();
    state.errorHint.clear();

    renderer.setAnalysis(&analysis, now());
    scene.setAnalysis(&analysis);
    panels.resetPanelScroll();
}

void Window::loadEntry(const chem::DbEntry& entry) {
    state.entry = &entry;
    state.customSmiles.clear();
    customSmiles.clear();
    state.notice.clear();
    applyAnalysis(chem::analyze(chem::moleculeFromEntry(entry)));

    SetWindowTextW(edit, toWide(entry.name).c_str());
    state.suggestions.clear();
    state.suggestIndex = -1;
    panels.revealActive();
}

void Window::loadSmiles(const std::string& smiles) {
    try {
        chem::Molecule molecule = chem::parseSmiles(smiles);
        state.entry = nullptr;
        customSmiles = smiles;
        state.customSmiles = smiles;
        state.notice.clear();
        applyAnalysis(chem::analyze(std::move(molecule)));
    } catch (const std::exception& err) {
        showError("Не удалось построить молекулу", err.what());
    }
}

void Window::showError(const std::string& message, const std::string& hint) {
    state.errorMessage = message;
    state.errorHint = hint;
    panels.resetPanelScroll();
}

void Window::submitQuery() {
    wchar_t buffer[512];
    GetWindowTextW(edit, buffer, 512);
    const std::string query = toUtf8(buffer);

    state.suggestions.clear();
    state.suggestIndex = -1;
    if (query.empty()) return;

    const chem::Resolution result = chem::resolveQuery(query);
    switch (result.kind) {
        case chem::Resolution::Kind::Entries: {
            loadEntry(*result.entries[0]);
            if (result.entries.size() > 1 && result.reason == chem::Resolution::Reason::Formula) {
                // брутто-формула не определяет строение однозначно — это и есть изомерия
                state.notice = "Формуле " + result.entries[0]->formulaPretty + " отвечает "
                             + std::to_string(result.entries.size())
                             + " вещества. Брутто-формула не определяет строение однозначно — "
                               "это и есть структурная изомерия. Остальные изомеры перечислены ниже.";
            }
            break;
        }
        case chem::Resolution::Kind::Smiles:
            loadSmiles(result.smiles);
            break;
        case chem::Resolution::Kind::Error:
            showError(result.message, result.hint);
            break;
    }
}

void Window::updateSuggestions() {
    wchar_t buffer[512];
    GetWindowTextW(edit, buffer, 512);
    const std::string query = toUtf8(buffer);

    state.suggestions = query.empty() ? std::vector<const chem::DbEntry*>()
                                      : chem::suggestions(query);
    state.suggestIndex = -1;
}

}  // namespace ui
