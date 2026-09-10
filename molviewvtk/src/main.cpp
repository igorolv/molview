// Точка входа оконной программы.
//
// Здесь запускается и останавливается GDI+, а вся работа идёт в ui::Window.
// Ядро расчёта (src/chem) про Windows не знает ничего и собирается отдельно —
// его проверяют консольные тесты.

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#include <windows.h>
#include <objidl.h>
#include <gdiplus.h>

#include <ctime>
#include <exception>

#include "chem/json.h"
#include "render/draw.h"
#include "ui/window.h"

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int showCommand) {
    // Программа сама масштабирует раскладку, поэтому просит систему
    // не растягивать окно за неё — иначе на мониторах с крупным
    // масштабированием картинка получится размытой.
    SetProcessDPIAware();
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    // Каталог data ищется рядом с .exe, а не относительно рабочего каталога:
    // у ярлыка он может быть каким угодно.
    wchar_t exePath[MAX_PATH] = {};
    if (GetModuleFileNameW(nullptr, exePath, MAX_PATH) > 0) {
        std::wstring path(exePath);
        const std::size_t slash = path.find_last_of(L"\\/");
        if (slash != std::wstring::npos) {
            chem::setProgramDirectory(render::toUtf8(path.substr(0, slash)));
        }
    }

    Gdiplus::GdiplusStartupInput startupInput;
    ULONG_PTR gdiplusToken = 0;
    if (Gdiplus::GdiplusStartup(&gdiplusToken, &startupInput, nullptr) != Gdiplus::Ok) {
        MessageBoxW(nullptr, L"Не удалось запустить GDI+.", L"Строение молекул", MB_ICONERROR);
        return 1;
    }

    int code = 1;
    try {
        ui::Window window;
        if (window.create(instance, showCommand)) {
            code = window.run();
        } else {
            MessageBoxW(nullptr, L"Не удалось создать окно.", L"Строение молекул", MB_ICONERROR);
        }
    } catch (const std::exception& err) {
        // чаще всего сюда попадает отсутствие файлов из каталога data
        const std::wstring text = L"Сбой при запуске:\n" + render::toWide(err.what())
            + L"\n\nПроверьте, что рядом с программой лежит каталог data "
              L"с файлами elements.json и molecules.json.";
        MessageBoxW(nullptr, text.c_str(), L"Строение молекул", MB_ICONERROR);
    }

    Gdiplus::GdiplusShutdown(gdiplusToken);
    return code;
}
