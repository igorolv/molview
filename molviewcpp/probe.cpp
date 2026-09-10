// Проба из раздела 3 задания: собирается ли GDI+ этим компилятором.
// Если здесь получится 0 — второй этап (окно и рендер) технически возможен.
#define NOMINMAX
#define UNICODE
#include <windows.h>
#include <objidl.h>
#include <gdiplus.h>
#include <cstdio>

int main() {
    Gdiplus::GdiplusStartupInput in;
    ULONG_PTR token;
    Gdiplus::Status s = Gdiplus::GdiplusStartup(&token, &in, nullptr);
    std::printf("GdiplusStartup -> %d\n", (int)s);
    Gdiplus::GdiplusShutdown(token);
    return 0;
}
