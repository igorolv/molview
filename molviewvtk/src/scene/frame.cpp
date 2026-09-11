#include "scene/frame.h"

#include <vtkPNGWriter.h>
#include <vtkWindowToImageFilter.h>

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

    // Затенение в углублениях между атомами. Главный источник ощущения
    // объёма и единственное, чего принципиально не мог собственный рендер:
    // там объём изображался нарисованным бликом, не зависящим от соседей.
    // Радиус выборки задаётся в ангстремах и подбирается под стиль модели,
    // см. Scene::Impl::buildRadii.
    ren->SetUseSSAO(true);
    ren->SetSSAOBias(0.02);
    ren->SetSSAOKernelSize(64);
    ren->SetSSAOBlur(true);
    // Края шаров без сглаживания заметно лестничные, а вторая версия
    // сглаживала их средствами GDI+ — без этого третья выглядела бы хуже
    // предшественницы.
    ren->SetUseFXAA(true);

    // Честная прозрачность. На этапе 3 её включать было незачем: в сцене не
    // было ни одного прозрачного предмета. Теперь их два — заливка сектора
    // валентного угла и изоповерхности орбиталей, причём лепестки одного
    // атома пересекаются между собой, и без послойного разбора видно то,
    // что нарисовалось последним, а не то, что ближе. Слоёв восемь — столько
    // и названо в разделе 7 задания. У четырёх пересекающихся лепестков sp³
    // луч проходит и больше поверхностей, но разницы уже не видно, а каждый
    // лишний слой — это ещё одна отрисовка всей сцены.
    ren->SetUseDepthPeeling(true);
    ren->SetMaximumNumberOfPeels(8);
    // Доля пикселей, на которой разбор останавливается досрочно. По
    // умолчанию ноль, то есть все восемь слоёв считаются всегда; сотая
    // доля отсекает хвост, которого на экране всё равно не различить.
    ren->SetOcclusionRatio(0.01);
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

bool Frame::saveImage(const std::string& utf8Path, int scale, unsigned char r, unsigned char g,
                      unsigned char b) {
    if (scale < 1) scale = 1;
    const int keepWidth = frameWidth;
    const int keepHeight = frameHeight;

    // Прозрачный фон годится для наложения на сияние от GDI+, но не для файла.
    ren->SetBackground(r / 255.0, g / 255.0, b / 255.0);
    ren->SetBackgroundAlpha(1.0);
    window->SetSize(keepWidth * scale, keepHeight * scale);
    window->Render();

    vtkNew<vtkWindowToImageFilter> grab;
    grab->SetInput(window);
    grab->SetInputBufferTypeToRGB();
    grab->ReadFrontBufferOff();
    grab->Update();

    vtkNew<vtkPNGWriter> writer;
    writer->SetFileName(utf8Path.c_str());
    writer->SetInputConnection(grab->GetOutputPort());
    writer->Write();
    const bool ok = writer->GetErrorCode() == 0;

    ren->SetBackground(0.0, 0.0, 0.0);
    ren->SetBackgroundAlpha(0.0);
    window->SetSize(keepWidth, keepHeight);
    // Кадр на экране остался от увеличенного рендера — перерисовываем.
    render();
    return ok;
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
