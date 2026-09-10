// ПРОБА. Запускать первой, до всякой архитектуры.
//
// Она проверяет разом всё, что может не заработать, и ради чего не стоит
// писать неделю кода вслепую:
//
//   1. VTK из пакета MSYS2 подключается и линкуется тем же g++;
//   2. модуль DomainsChemistry на месте — vtkMolecule и vtkMoleculeMapper
//      собирают шаростержневую модель;
//   3. vtk_module_autoinit отработал: без него программа соберётся,
//      запустится и покажет ЧЁРНЫЙ ПРЯМОУГОЛЬНИК без единой ошибки,
//      потому что фабрики объектов OpenGL не зарегистрируются;
//   4. закадровый рендер (SetShowWindow(false) + буферы в FBO) работает
//      на этой видеокарте;
//   5. кадр удаётся забрать из VTK и вывести в обычное окно Win32 через
//      StretchDIBits — это и есть схема, на которой держится вся версия:
//      VTK рисует сцену, GDI+ рисует интерфейс ПОВЕРХ неё;
//   6. и всё это укладывается в бюджет кадра. Замеры печатаются в заголовок
//      окна: рендер, чтение кадра, вывод.
//
// Ожидается: окно с вращающимся метаном, затенение в углублениях между
// атомами (SSAO), суммарное время кадра заметно меньше 16 мс.
//
// Если пункт 4 не заработает — отступать к дочернему окну: убрать закадровый
// режим и позвать SetParentId(hwnd) у vtkWin32OpenGLRenderWindow. Тогда VTK
// заведёт собственное дочернее окно, но рисовать поверх него уже не выйдет,
// и панель инструментов придётся выносить из области сцены.

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#include <windows.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkLightKit.h>
#include <vtkMolecule.h>
#include <vtkMoleculeMapper.h>
#include <vtkNew.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>

#include <cwchar>
#include <vector>

namespace {

vtkNew<vtkRenderer> renderer;
vtkNew<vtkRenderWindow> renderWindow;
vtkNew<vtkUnsignedCharArray> pixels;

/** Кадр в том виде, в каком его понимает StretchDIBits: 32 бита, BGRA. */
std::vector<unsigned char> frame;
int frameWidth = 0;
int frameHeight = 0;

double milliseconds(LARGE_INTEGER from, LARGE_INTEGER to) {
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return double(to.QuadPart - from.QuadPart) * 1000.0 / double(freq.QuadPart);
}

/** Метан: один углерод и четыре водорода по вершинам тетраэдра. */
void buildScene() {
    vtkNew<vtkMolecule> molecule;
    const vtkIdType c = molecule->AppendAtom(6, 0.0, 0.0, 0.0).GetId();
    const double d = 0.63;  // 1.09 Å по диагонали куба
    const double corners[4][3] = {{d, d, d}, {d, -d, -d}, {-d, d, -d}, {-d, -d, d}};
    for (const auto& p : corners) {
        const vtkIdType h = molecule->AppendAtom(1, p[0], p[1], p[2]).GetId();
        molecule->AppendBond(c, h, 1);
    }

    vtkNew<vtkMoleculeMapper> mapper;
    mapper->SetInputData(molecule);
    mapper->UseBallAndStickSettings();

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    renderer->AddActor(actor);

    // Тот самый набор, ради которого всё и затевается: объёмное затенение
    // в щелях между атомами и трёхточечная схема света вместо одной лампы.
    renderer->SetUseSSAO(true);
    renderer->SetSSAORadius(0.5);
    renderer->SetSSAOBias(0.01);
    renderer->SetSSAOKernelSize(32);
    renderer->SetSSAOBlur(true);
    renderer->SetUseFXAA(true);
    renderer->SetBackground(0.043, 0.055, 0.086);  // тот же тёмный фон, что в оригинале

    vtkNew<vtkLightKit> lights;
    lights->AddLightsToRenderer(renderer);

    renderWindow->AddRenderer(renderer);
    // Ключевая пара: окна VTK на экране нет, картинка живёт в кадровом буфере.
    renderWindow->SetShowWindow(false);
    renderWindow->SetUseOffScreenBuffers(true);
    renderWindow->SetMultiSamples(0);  // сглаживание берём с FXAA, MSAA поверх FBO капризнее
}

/** Рендер, чтение кадра и вывод. Возвращает замеры в миллисекундах. */
void drawFrame(HWND window, HDC dc, int width, int height) {
    if (width <= 0 || height <= 0) return;

    LARGE_INTEGER t0, t1, t2, t3;
    QueryPerformanceCounter(&t0);

    renderWindow->SetSize(width, height);

    // Камеру наводим ПОСЛЕ первого рендера. До него маппер геометрию ещё
    // не построил, границы сцены пустые, и ResetCamera наводится в никуда —
    // молекула выходит крошечной посреди пустого кадра.
    static bool cameraReady = false;
    if (!cameraReady) {
        renderWindow->Render();
        renderer->ResetCamera();
        renderer->GetActiveCamera()->Zoom(1.5);
        cameraReady = true;
    }

    renderer->GetActiveCamera()->Azimuth(0.6);
    renderWindow->Render();
    QueryPerformanceCounter(&t1);

    // front = 0: читаем задний буфер, тот, в который только что нарисовали.
    renderWindow->GetRGBACharPixelData(0, 0, width - 1, height - 1, 0, pixels);
    QueryPerformanceCounter(&t2);

    if (frameWidth != width || frameHeight != height) {
        frame.assign(std::size_t(width) * std::size_t(height) * 4, 0);
        frameWidth = width;
        frameHeight = height;
    }

    // OpenGL отдаёт RGBA снизу вверх; DIB с положительной высотой тоже идёт
    // снизу вверх, поэтому переворачивать ничего не нужно — только поменять
    // местами красный и синий.
    const unsigned char* src = pixels->GetPointer(0);
    unsigned char* dst = frame.data();
    const std::size_t count = std::size_t(width) * std::size_t(height);
    for (std::size_t i = 0; i < count; i++) {
        dst[i * 4 + 0] = src[i * 4 + 2];
        dst[i * 4 + 1] = src[i * 4 + 1];
        dst[i * 4 + 2] = src[i * 4 + 0];
        dst[i * 4 + 3] = 255;
    }

    BITMAPINFO info = {};
    info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    info.bmiHeader.biWidth = width;
    info.bmiHeader.biHeight = height;  // положительная — строки снизу вверх
    info.bmiHeader.biPlanes = 1;
    info.bmiHeader.biBitCount = 32;
    info.bmiHeader.biCompression = BI_RGB;
    StretchDIBits(dc, 0, 0, width, height, 0, 0, width, height, frame.data(), &info,
                  DIB_RGB_COLORS, SRCCOPY);
    QueryPerformanceCounter(&t3);

    wchar_t title[160];
    std::swprintf(title, 160,
                  L"Проба VTK — рендер %.1f мс · чтение %.1f мс · вывод %.1f мс · всего %.1f мс",
                  milliseconds(t0, t1), milliseconds(t1, t2), milliseconds(t2, t3),
                  milliseconds(t0, t3));
    SetWindowTextW(window, title);
}

LRESULT CALLBACK windowProc(HWND window, UINT message, WPARAM wParam, LPARAM lParam) {
    switch (message) {
        case WM_PAINT: {
            PAINTSTRUCT ps;
            const HDC dc = BeginPaint(window, &ps);
            RECT client;
            GetClientRect(window, &client);
            drawFrame(window, dc, client.right, client.bottom);
            EndPaint(window, &ps);
            return 0;
        }
        case WM_TIMER:
            InvalidateRect(window, nullptr, FALSE);
            return 0;
        case WM_ERASEBKGND:
            return 1;  // фон закрывается кадром целиком, стирать нечего
        case WM_DESTROY:
            PostQuitMessage(0);
            return 0;
        default:
            return DefWindowProcW(window, message, wParam, lParam);
    }
}

}  // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int showCommand) {
    SetProcessDPIAware();
    buildScene();

    WNDCLASSEXW cls = {};
    cls.cbSize = sizeof(cls);
    cls.lpfnWndProc = windowProc;
    cls.hInstance = instance;
    cls.hCursor = LoadCursor(nullptr, IDC_ARROW);
    cls.lpszClassName = L"MolviewVtkProbe";
    RegisterClassExW(&cls);

    const HWND window = CreateWindowExW(0, cls.lpszClassName, L"Проба VTK",
                                        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
                                        900, 700, nullptr, nullptr, instance, nullptr);
    if (!window) return 1;

    ShowWindow(window, showCommand);
    SetTimer(window, 1, 16, nullptr);  // ~60 кадров в секунду

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }
    return 0;
}
