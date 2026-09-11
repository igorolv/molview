#pragma once

/**
 * ЗАКАДРОВЫЙ КАДР VTK И ЕГО ПЕРЕНОС В РАСТР GDI
 *
 * VTK рисует не в окно, а в кадровый буфер; готовый кадр накладывается на тот
 * же растр, в который рисует GDI+. Почему так, а не дочерним окном — раздел 5
 * задания: поверх дочернего окна Windows не даст ничего нарисовать, а панель
 * инструментов и заголовок молекулы лежат поверх сцены.
 *
 * Фон сцены (сияние и сетка) остаётся за GDI+, поэтому VTK рисует на
 * ПРОЗРАЧНОМ фоне, а кадр накладывается через AlphaBlend. Отсюда требование
 * предумноженной альфы — так устроен AC_SRC_ALPHA.
 */

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#include <windows.h>

#include <vtkNew.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>

namespace scene {

class Frame {
public:
    Frame();
    ~Frame();
    Frame(const Frame&) = delete;
    Frame& operator=(const Frame&) = delete;

    /** Рендерер, в который складывается сцена. */
    vtkRenderer* renderer() { return ren; }

    /** Размер кадра в пикселях. Повторный вызов с тем же размером бесплатен. */
    void resize(int width, int height);

    /** Нарисовать сцену в закадровый буфер и забрать её в растр. */
    void render();

    /** Наложить последний кадр на растр приёмника в точке (x, y). */
    void blendTo(HDC target, int x, int y);

    /**
     * Снимок для печати: кадр рисуется увеличенным в scale раз и на
     * НЕПРОЗРАЧНОМ фоне — иначе в файле вместо сияния будет дыра.
     * Размер и прозрачность восстанавливаются, путь в UTF-8.
     */
    bool saveImage(const std::string& utf8Path, int scale, unsigned char r, unsigned char g,
                   unsigned char b);

    int width() const { return frameWidth; }
    int height() const { return frameHeight; }

private:
    void releaseSurface();

    vtkNew<vtkRenderer> ren;
    vtkNew<vtkRenderWindow> window;
    /** Буфер чтения из OpenGL: RGBA, строки снизу вверх. */
    vtkNew<vtkUnsignedCharArray> pixels;

    // Растр, к которому применяется AlphaBlend. Пишем прямо в его пиксели,
    // поэтому промежуточный буфер не нужен.
    HDC surfaceDc = nullptr;
    HBITMAP surface = nullptr;
    HBITMAP surfaceOld = nullptr;
    unsigned char* surfaceBits = nullptr;

    int frameWidth = 0;
    int frameHeight = 0;
    bool ready = false;
};

}  // namespace scene
