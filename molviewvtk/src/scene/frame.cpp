#include "scene/frame.h"

#include <algorithm>

namespace scene {

Frame::Frame() {
    window->AddRenderer(ren);
    // Окна на экране нет: картинка живёт в кадровом буфере и забирается оттуда.
    window->SetShowWindow(false);
    window->SetUseOffScreenBuffers(true);
    // Без плоскостей альфы фон вернётся непрозрачным и закроет сияние,
    // которое рисует GDI+ под сценой.
    window->SetAlphaBitPlanes(1);
    window->SetMultiSamples(0);

    ren->SetBackground(0.0, 0.0, 0.0);
    ren->SetBackgroundAlpha(0.0);
}

Frame::~Frame() { releaseSurface(); }

void Frame::releaseSurface() {
    if (surfaceDc == nullptr) return;
    SelectObject(surfaceDc, surfaceOld);
    DeleteObject(surface);
    DeleteDC(surfaceDc);
    surfaceDc = nullptr;
    surface = nullptr;
    surfaceOld = nullptr;
    surfaceBits = nullptr;
}

void Frame::resize(int width, int height) {
    width = std::max(1, width);
    height = std::max(1, height);
    if (width == frameWidth && height == frameHeight && surfaceDc != nullptr) return;

    releaseSurface();
    frameWidth = width;
    frameHeight = height;
    ready = false;

    BITMAPINFO info = {};
    info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    info.bmiHeader.biWidth = frameWidth;
    // Положительная высота — строки снизу вверх, ровно как их отдаёт OpenGL.
    // Переворачивать кадр не нужно; это совпадение, а не случайность в коде,
    // и «поправить» его нельзя.
    info.bmiHeader.biHeight = frameHeight;
    info.bmiHeader.biPlanes = 1;
    info.bmiHeader.biBitCount = 32;
    info.bmiHeader.biCompression = BI_RGB;

    surfaceDc = CreateCompatibleDC(nullptr);
    void* bits = nullptr;
    surface = CreateDIBSection(surfaceDc, &info, DIB_RGB_COLORS, &bits, nullptr, 0);
    surfaceBits = static_cast<unsigned char*>(bits);
    surfaceOld = static_cast<HBITMAP>(SelectObject(surfaceDc, surface));

    window->SetSize(frameWidth, frameHeight);
}

void Frame::render() {
    if (surfaceBits == nullptr) return;
    window->SetSize(frameWidth, frameHeight);
    window->Render();

    // front = 0: читаем задний буфер, тот, в который только что нарисовали.
    window->GetRGBACharPixelData(0, 0, frameWidth - 1, frameHeight - 1, 0, pixels);
    if (pixels->GetNumberOfTuples() < static_cast<vtkIdType>(frameWidth) * frameHeight) return;

    const unsigned char* src = pixels->GetPointer(0);
    for (int row = 0; row < frameHeight; row++) {
        const unsigned char* s = src + static_cast<std::size_t>(row) * frameWidth * 4;
        unsigned char* d = surfaceBits + static_cast<std::size_t>(row) * frameWidth * 4;
        for (int x = 0; x < frameWidth; x++) {
            // ЗЕРКАЛО ПО ГОРИЗОНТАЛИ. Собственный рендер кладёт мировой +x
            // ВПРАВО, а правосторонняя камера при взгляде с −z кладёт его
            // влево: right = forward × up = z × y = −x. Расхождение в одну
            // ось, то есть зеркало. Разворачиваем здесь, чтобы картинка
            // совпала с наложениями GDI+ и с прежними версиями программы,
            // и чтобы не вносить отражение в матрицу — оно вывернуло бы
            // нормали и испортило освещение.
            const unsigned char* p = s + static_cast<std::size_t>(frameWidth - 1 - x) * 4;
            unsigned char* q = d + static_cast<std::size_t>(x) * 4;
            const unsigned int alpha = p[3];
            // AC_SRC_ALPHA требует предумноженных цветов, иначе полупрозрачные
            // края светятся.
            q[0] = static_cast<unsigned char>(p[2] * alpha / 255);
            q[1] = static_cast<unsigned char>(p[1] * alpha / 255);
            q[2] = static_cast<unsigned char>(p[0] * alpha / 255);
            q[3] = static_cast<unsigned char>(alpha);
        }
    }
    ready = true;
}

void Frame::blendTo(HDC target, int x, int y) {
    if (!ready || surfaceDc == nullptr) return;
    BLENDFUNCTION blend = {};
    blend.BlendOp = AC_SRC_OVER;
    blend.SourceConstantAlpha = 255;
    blend.AlphaFormat = AC_SRC_ALPHA;
    AlphaBlend(target, x, y, frameWidth, frameHeight, surfaceDc, 0, 0, frameWidth, frameHeight,
               blend);
}

}  // namespace scene
