#pragma once

/**
 * СЦЕНА НА VTK: шары и связи.
 *
 * Здесь живёт всё, что знает про VTK, кроме переноса кадра — тот в frame.h.
 * Заголовок нарочно не тянет заголовки VTK: интерфейсу про неё знать нечего,
 * а компиляция ui/ от этого не разбухает вдвое.
 *
 * Камера не выбирается самостоятельно, а берётся у собственного рендера
 * (render::CameraParams). Иначе шары разъедутся с дугами углов и подписями,
 * которые пока рисует GDI+ поверх кадра.
 */

#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef UNICODE
#define UNICODE
#endif
#include <windows.h>

#include <memory>
#include <string>

#include "chem/types.h"
#include "render/renderer.h"

namespace scene {

class Scene {
public:
    Scene();
    ~Scene();
    Scene(const Scene&) = delete;
    Scene& operator=(const Scene&) = delete;

    /** Новая молекула. nullptr — сцена пустеет. */
    void setAnalysis(const chem::Analysis* analysis);

    /**
     * Нарисовать сцену и наложить её на растр приёмника.
     * Камеру, стиль модели и радиусы шаров берёт у рендера — там они
     * единственный источник истины.
     */
    void draw(HDC target, const render::MoleculeRenderer& view);

    /**
     * Снимок сцены в PNG для записки и презентации. Путь выбирается сам —
     * каталог «Изображения» и имя по названию молекулы; он же возвращается,
     * чтобы программе было что показать. Пусто — не получилось.
     */
    std::wstring saveImage(const render::MoleculeRenderer& view);

    /** Замеры последнего кадра в миллисекундах — для отладки. */
    double lastRenderMs() const;
    double lastTransferMs() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

}  // namespace scene
