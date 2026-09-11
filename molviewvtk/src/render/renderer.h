#pragma once

/**
 * СОБСТВЕННЫЙ ТРЁХМЕРНЫЙ РЕНДЕР НА GDI+
 *
 * Никаких сторонних библиотек и никакого OpenGL: перспективная проекция
 * считается вручную, а объём создаётся сортировкой примитивов по глубине
 * (алгоритм художника), радиальными градиентами на шарах и продольными
 * градиентами на связях. От GDI+ нужны только заливки, градиенты, пути и текст.
 */

#include <functional>
#include <string>
#include <vector>

#include "chem/orbital.h"
#include "chem/types.h"
#include "chem/vec.h"
#include "render/draw.h"
#include "render/occlusion.h"
#include "render/palette.h"

namespace render {

enum class AngleMode { None, Selected, All };
enum class Style { BallStick, SpaceFill, Wire };

/**
 * Какую радиальную часть показывать у орбиталей.
 *
 * Textbook — приближение Слейтера: лепесток вперёд, хвостик назад, четыре
 * sp³ читаются как тетраэдр. Это картинка из учебника, и она стоит по
 * умолчанию: программа учебная, и понятность здесь важнее точности.
 * Computed — водородоподобная 2s с радиальным узлом, как записано в задании.
 * Посчитана верно, а смотреть на неё нельзя: положительная часть прячется
 * в ядре, отрицательная торчит конусами (раздел 11 задания). Режим оставлен
 * нарочно — чтобы показать, ПОЧЕМУ в учебниках рисуют не это.
 */
enum class OrbitalView { Textbook, Computed };

struct ViewOptions {
    bool showLonePairs = true;
    bool showOrbitals = false;
    OrbitalView orbitalView = OrbitalView::Textbook;
    bool showLabels = true;
    bool showDipole = false;
    AngleMode showAngles = AngleMode::Selected;
    bool autoRotate = true;
    Style style = Style::BallStick;
};

/**
 * Камера, которой пользуется собственный рендер. Отдаётся наружу, чтобы сцена
 * на VTK смотрела ровно туда же: иначе шары разъедутся с дугами углов и
 * подписями, которые по-прежнему рисует GDI+.
 */
struct CameraParams {
    /** Поворот мира в систему камеры. */
    chem::Mat3 rotation = chem::IDENTITY;
    /** Расстояние от глаза до центра сцены, Å. */
    double distance = 0;
    /** Фокусное расстояние в пикселях. */
    double focal = 0;
    /** Радиус сцены, Å — по нему считаются плоскости отсечения. */
    double sceneRadius = 0;
    RectD viewport;
};

/**
 * Дуга валентного угла. Геометрию считает вид, рисует её сцена на VTK — там
 * дуга правильно уходит за шар центрального атома, чего при рисовании поверх
 * кадра быть не может. Подпись остаётся за GDI+, поэтому расчёт общий.
 */
struct AngleArc {
    const chem::AngleRecord* record = nullptr;
    chem::Vec3 center;
    /** Единичные направления на соседей. */
    chem::Vec3 from, to;
    /** Точки дуги в мировых координатах — по ним сцена строит трубку. */
    std::vector<chem::Vec3> points;
    double radius = 0;
    /** 180°: дуги не существует, вместо неё отрезок-диаметр. */
    bool linear = false;
    /** Угол при выбранном или подсвеченном атоме — рисуется ярче. */
    bool emphasised = false;
    Rgb color;
    double alpha = 1;
};

/**
 * Гибридные орбитали одного атома — задание для сцены: она считает по ним
 * волновые функции и проводит изоповерхности. Отбор атомов делает вид, там
 * же, где он делается для дуг углов.
 *
 * Направления собраны в один набор не для удобства, а ради скорости:
 * chem::orbitalGrids считает радиальную часть один раз на все направления
 * атома, и это втрое быстрее, чем считать каждый гибрид отдельно.
 */
struct OrbitalSet {
    /** Центр атома, Å. */
    chem::Vec3 center;
    /** Тип гибрида — от него зависят коэффициенты при s и p. */
    chem::Orbital kind = chem::Orbital::Sp3;
    /** Единичные направления гибридов: сначала связывающие, за ними — занятые парами. */
    std::vector<chem::Vec3> dirs;
    /** Сколько первых направлений в dirs — связи; остальные — неподелённые пары. */
    int bonding = 0;
    /** Размер сетки, заряд ядра, вид радиальной части. */
    chem::OrbitalOptions options;
};

/** Облако неподелённой пары. */
struct LonePairCloud {
    chem::Vec3 position;
    /** Направление от атома — вдоль него разносятся два электрона. */
    chem::Vec3 direction;
    double radius = 0;
    double alpha = 1;
};

class MoleculeRenderer {
public:
    /** Разбор принадлежит приложению; рендер только читает его. */
    void setAnalysis(const chem::Analysis* analysis, double timeMs, bool animate = true);
    void setViewport(const RectD& rect) { viewport = rect; }
    const RectD& area() const { return viewport; }

    ViewOptions& options() { return opts; }
    const ViewOptions& options() const { return opts; }

    void resetView();
    /** Поворот по таймеру: dt в миллисекундах. */
    void tick(double dtMs);

    void render(Gdiplus::Graphics& g, Fonts& fonts, double timeMs);

    /**
     * Кадр, разрезанный надвое: между этими двумя вызовами сцену рисует VTK.
     * Первый готовит фон и пересчитывает проекции атомов, второй кладёт
     * поверх всё, что VTK пока не рисует.
     */
    void renderBackground(Gdiplus::Graphics& g, double timeMs);
    void renderOverlays(Gdiplus::Graphics& g, Fonts& fonts, double timeMs);

    /** Камера для внешнего рендера. */
    CameraParams camera() const;

    /** Дуги углов для сцены. Считаются по режиму показа и текущему выбору. */
    std::vector<AngleArc> angleArcs() const;

    /**
     * Гибридные орбитали для сцены — по одному набору на атом.
     *
     * Здесь же выбираются два параметра ПОКАЗА, которых нет в задании и не
     * должно быть в ядре: эффективный заряд ядра (размер орбитали) и вид
     * радиальной части 2s (её форма). Оба разобраны у chem::OrbitalOptions
     * и chem::Radial; коротко — заряд подбирается так, чтобы радиус
     * наибольшей плотности p-части совпал с ковалентным радиусом элемента,
     * а радиальная часть берётся слейтеровская, потому что водородоподобная
     * даёт верную, но нерисуемую картинку.
     */
    std::vector<OrbitalSet> orbitals() const;

    /**
     * Радиус шара, Å. Наружу нужен затем, чтобы VTK строила шары ровно того же
     * размера: подписи и дуги углов рассчитаны на эти радиусы.
     */
    double atomRadius(const std::string& el) const;

    // --- мышь ---
    void beginDrag(double x, double y);
    void dragTo(double x, double y);
    /** Возвращает true, если это был щелчок без перетаскивания. */
    bool endDrag();
    void wheel(double delta);
    /** Атом под курсором или -1. */
    int hitTest(double px, double py) const;

    void setHovered(int id) { hoveredAtom = id; }
    void setSelected(int id) { selectedAtom = id; }
    int hovered() const { return hoveredAtom; }
    int selected() const { return selectedAtom; }
    bool isDragging() const { return dragging; }

private:
    struct Projected {
        double x = 0, y = 0, z = 0;
        /** Пикселей на ангстрем на этой глубине. */
        double scale = 0;
        bool visible = false;
    };

    /**
     * Элемент сцены. Собирается в первом проходе, рисуется во втором —
     * после сортировки по глубине. Замыкание хранит всё, что нужно для
     * отрисовки, ровно как в оригинале.
     */
    struct Primitive {
        double depth = 0;
        std::function<void(Gdiplus::Graphics&)> draw;
    };

    /** Рисуются ли у этого атома гибридные орбитали. */
    bool hasOrbitals(const chem::AtomAnalysis& info) const;

    double sceneRadius() const;
    double cameraDistance() const;
    double focal() const;
    Projected project(const chem::Vec3& p) const;
    Rgb fog(const Rgb& color, double z) const;

    void drawBackground(Gdiplus::Graphics& g) const;
    void collectAtoms(std::vector<Primitive>& out, double progress, double timeMs) const;
    void collectBonds(std::vector<Primitive>& out, double progress) const;
    void collectLonePairs(std::vector<Primitive>& out, double progress) const;

    void drawAngles(Gdiplus::Graphics& g, Fonts& fonts, double progress);
    void drawLinearAngle(Gdiplus::Graphics& g, Fonts& fonts, const chem::AngleRecord& record,
                         const chem::Vec3& center, const chem::Vec3& u, const chem::Vec3& w,
                         double radius, double fade);
    void drawDipole(Gdiplus::Graphics& g, Fonts& fonts, double progress) const;
    void drawLabels(Gdiplus::Graphics& g, Fonts& fonts, double progress) const;

    bool labelCollides(const RectD& r) const;

    /** Пересчитать положение глаза, оси камеры и заслоняющие шары. */
    void updateOcclusion();

    const chem::Analysis* analysis = nullptr;
    ViewOptions opts;

    // Видимость подписей. Пересчитывается раз в кадр вместе с проекциями.
    Occluder occluder;
    chem::Vec3 eyePos, camRight, camUp;

    chem::Mat3 rotation = chem::IDENTITY;
    double zoom = 1;
    RectD viewport{0, 0, 800, 600};

    double entryStart = 0;
    /** Ход анимации появления в текущем кадре. Ставится в renderBackground. */
    double frameProgress = 1;
    std::vector<Projected> projected;
    /** Прямоугольники уже нарисованных подписей — чтобы они не налезали друг на друга. */
    std::vector<RectD> labelRects;

    int hoveredAtom = -1;
    int selectedAtom = -1;

    bool dragging = false;
    double lastPointerX = 0, lastPointerY = 0;
    bool pointerMoved = false;
};

}  // namespace render
