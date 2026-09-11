#include "scene/scene.h"

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkPointData.h>
#include <vtkCamera.h>
#include <vtkAppendPolyData.h>
#include <vtkBlueObeliskData.h>
#include <vtkFloatArray.h>
#include <vtkFlyingEdges3D.h>
#include <vtkImageData.h>
#include <vtkLightKit.h>
#include <vtkMolecule.h>
#include <vtkMoleculeMapper.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPeriodicTable.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>
#include <vtkUnsignedCharArray.h>

#include <shlobj.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "chem/orbital.h"
#include "chem/periodic.h"
#include "chem/vec.h"
#include "render/palette.h"
#include "render/draw.h"
#include "scene/frame.h"

namespace scene {

namespace {

constexpr double BOND_RADIUS = 0.115;  // как в собственном рендере

// --- орбитали ---------------------------------------------------------------
// Уровни двух изоповерхностей — в долях от наибольшего и от наименьшего
// значения ψ на этой же сетке, каждый от своего.
//
// Доли, а не абсолютные величины: у sp, sp² и sp³ разные коэффициенты при s,
// а значит и разная амплитуда. Порознь — потому что обратный лепесток слабее
// переднего вдвое у sp³ и вчетверо у sp: на общем уровне выше половины его
// нет вовсе, а ниже половины передние лепестки раздуваются в яйцо шире
// собственной длины.
//
// Числа подобраны под ПОНЯТНОСТЬ, а не под плотность: программа учебная.
// Передний лепесток при 0,7 — капля длиной 1,05 ковалентного радиуса и
// полушириной 0,55, кончик касается шара соседа. Обратный при 0,8 — хвостик
// длиной 1,0 радиуса, из шара атома (0,5 радиуса) выглядывает наполовину.
// При прежних 0,45 и 0,55 хвостик был почти с передний лепесток (1,26 против
// 1,42), сам передний — яйцо шире своей длины, и четыре sp³ сливались вокруг
// атома в сплошную кляксу.
constexpr double ORBITAL_ISO_PLUS = 0.70;
constexpr double ORBITAL_ISO_MINUS = 0.80;
// Обратный лепесток и полупрозрачным читается хорошо, а передний на нём же
// теряется — он крупнее и его больше перекрывают соседние.
constexpr double ORBITAL_OPACITY_PLUS = 0.55;
constexpr double ORBITAL_OPACITY_MINUS = 0.45;

double milliseconds(LARGE_INTEGER from, LARGE_INTEGER to) {
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return double(to.QuadPart - from.QuadPart) * 1000.0 / double(freq.QuadPart);
}

}  // namespace

struct Scene::Impl {
    Frame frame;
    vtkNew<vtkMolecule> molecule;
    vtkNew<vtkMoleculeMapper> mapper;
    vtkNew<vtkActor> actor;
    vtkNew<vtkLightKit> lights;

    // Дуги валентных углов. Собираются заново каждый кадр: их состав зависит
    // от того, какой атом выбран и над каким висит мышь.
    vtkNew<vtkPolyData> arcLines;
    vtkNew<vtkTubeFilter> arcTubes;
    vtkNew<vtkPolyDataMapper> arcMapper;
    vtkNew<vtkActor> arcActor;
    /** Заливка сектора между сторонами угла. */
    vtkNew<vtkPolyData> wedgeMesh;
    vtkNew<vtkPolyDataMapper> wedgeMapper;
    vtkNew<vtkActor> wedgeActor;

    /**
     * Изоповерхности гибридных орбиталей ОДНОГО атома.
     *
     * Три актёра на атом: передние лепестки связывающих орбиталей, передние
     * лепестки орбиталей с неподелёнными парами и обратные лепестки всех
     * разом. Цвет переднего лепестка — по назначению, а не по знаку ψ: для
     * школьника главное в картинке воды — «эти две смотрят на водороды, а
     * эти две заняты парами, потому и угол 104,5°, а не 109,5°». Знак при
     * этом не теряется: обратный лепесток всегда отдельного цвета.
     *
     * Гибриды одного цвета сшиваются в одного актёра через vtkAppendPolyData —
     * и дело не в опрятности. Прозрачность рисуется послойно, и каждый слой
     * заново обходит ВСЕХ прозрачных актёров: четыре орбитали воды — это
     * восемь актёров против трёх, то есть вдвое-втрое больше вызовов
     * отрисовки на тот же кадр.
     */
    struct SurfacePart {
        vtkNew<vtkAppendPolyData> mesh;
        vtkNew<vtkPolyDataMapper> mapper;
        vtkNew<vtkActor> actor;
    };
    struct OrbitalSurface {
        // Поля и контуры живут ровно затем, чтобы конвейер до них дотянулся;
        // vtkNew в вектор не кладётся, поэтому умные указатели.
        std::vector<vtkSmartPointer<vtkImageData>> fields;
        std::vector<vtkSmartPointer<vtkFlyingEdges3D>> contours;
        SurfacePart bond;
        SurfacePart pair;
        SurfacePart back;
        /** Сколько орбиталей попало в каждую сшивку — пустой актёр VTK не любит. */
        int bonds = 0;
        int pairs = 0;
    };
    // unique_ptr, а не сами объекты: vtkNew внутри не копируется и не
    // перемещается, поэтому в vector структура целиком не влезает.
    std::vector<std::unique_ptr<OrbitalSurface>> orbitals;
    /** Набор, под который построены поверхности — чтобы не считать их заново. */
    std::vector<render::OrbitalSet> builtOrbitals;

    const chem::Analysis* analysis = nullptr;
    /** Стиль, под который построены радиусы. -1 — ещё ни под какой. */
    int builtStyle = -1;
    bool builtMolecule = false;

    double renderMs = 0;
    double transferMs = 0;

    Impl() {
        mapper->SetInputData(molecule);
        mapper->SetRenderAtoms(true);
        mapper->SetRenderBonds(true);
        mapper->SetBondColorModeToDiscreteByAtom();
        mapper->SetUseMultiCylindersForBonds(true);
        // Радиусы задаются поштучно массивом «radii»: они должны совпасть
        // с теми, по которым собственный рендер размещает подписи и дуги.
        mapper->SetAtomicRadiusTypeToCustomArrayRadius();

        actor->SetMapper(mapper);
        frame.renderer()->AddActor(actor);
        lights->AddLightsToRenderer(frame.renderer());

        // Дуга — светящаяся линия, а не поверхность: свет ей ни к чему,
        // поэтому весь цвет идёт через ambient.
        arcTubes->SetInputData(arcLines);
        arcTubes->SetNumberOfSides(8);
        arcTubes->CappingOn();
        arcMapper->SetInputConnection(arcTubes->GetOutputPort());
        arcMapper->SetColorModeToDirectScalars();
        arcActor->SetMapper(arcMapper);
        arcActor->GetProperty()->SetAmbient(1.0);
        arcActor->GetProperty()->SetDiffuse(0.0);
        arcActor->GetProperty()->SetSpecular(0.0);
        frame.renderer()->AddActor(arcActor);

        wedgeMapper->SetInputData(wedgeMesh);
        wedgeMapper->SetColorModeToDirectScalars();
        wedgeActor->SetMapper(wedgeMapper);
        wedgeActor->GetProperty()->SetAmbient(1.0);
        wedgeActor->GetProperty()->SetDiffuse(0.0);
        wedgeActor->GetProperty()->SetSpecular(0.0);
        frame.renderer()->AddActor(wedgeActor);
    }

    /** Перенос графа молекулы в vtkMolecule. */
    void buildMolecule() {
        molecule->Initialize();
        builtStyle = -1;  // радиусы тоже придётся переложить
        if (analysis == nullptr) return;

        applyPalette();

        std::vector<vtkIdType> ids;
        ids.reserve(analysis->molecule.atoms.size());
        for (const chem::Atom& atom : analysis->molecule.atoms) {
            const int z = chem::element(atom.el).z;
            ids.push_back(molecule
                              ->AppendAtom(static_cast<unsigned short>(z > 0 ? z : 1), atom.pos.x,
                                           atom.pos.y, atom.pos.z)
                              .GetId());
        }
        for (const chem::Bond& bond : analysis->molecule.bonds) {
            if (bond.a < 0 || bond.b < 0) continue;
            if (static_cast<std::size_t>(bond.a) >= ids.size()) continue;
            if (static_cast<std::size_t>(bond.b) >= ids.size()) continue;
            // Кратность в ядре вещественная: 1,5 у ароматической связи,
            // 1⅓ в нитрат-ионе. VTK принимает целое и рисует только 1, 2, 3,
            // поэтому дробные округляются вниз, а признак delocalized
            // отрабатывается отдельной отрисовкой — это этап 4.
            const int order = std::max(1, std::min(3, static_cast<int>(bond.order)));
            molecule->AppendBond(ids[bond.a], ids[bond.b], static_cast<unsigned short>(order));
        }
        builtMolecule = true;
    }

    /**
     * Своя палитра вместо стандартной CPK.
     *
     * Подменить таблицу цветов у маппера нельзя — она внутренняя. Но
     * vtkPeriodicTable отдаёт наружу сам справочник vtkBlueObeliskData, а в
     * нём массив DefaultColors изменяемый. Пишем туда цвета из elements.json:
     * они те же CPK, но подогнанные под тёмный фон, и ими же панель разбора
     * красит значок элемента — расхождение было бы заметно.
     */
    void applyPalette() {
        if (analysis == nullptr) return;
        vtkFloatArray* colors = mapper->GetPeriodicTable()->GetBlueObeliskData()->GetDefaultColors();
        if (colors == nullptr) return;
        for (const chem::Atom& atom : analysis->molecule.atoms) {
            const int z = chem::element(atom.el).z;
            if (z <= 0 || z >= colors->GetNumberOfTuples()) continue;
            const render::Rgb c = render::hexToRgb(chem::element(atom.el).color);
            colors->SetTuple3(z, c.r / 255.0, c.g / 255.0, c.b / 255.0);
        }
        colors->Modified();
    }

    /** Радиусы шаров под текущий стиль. */
    void buildRadii(const render::MoleculeRenderer& view) {
        if (analysis == nullptr) return;
        vtkNew<vtkFloatArray> radii;
        radii->SetName("radii");
        radii->SetNumberOfComponents(1);
        radii->SetNumberOfTuples(static_cast<vtkIdType>(analysis->molecule.atoms.size()));
        for (std::size_t i = 0; i < analysis->molecule.atoms.size(); i++) {
            radii->SetValue(static_cast<vtkIdType>(i),
                            static_cast<float>(view.atomRadius(analysis->molecule.atoms[i].el)));
        }
        molecule->GetAtomData()->AddArray(radii);
        // конвейер VTK смотрит на MTime ВХОДА; без этого маппер
        // не заметит, что массив появился
        molecule->Modified();

        const render::Style style = view.options().style;
        mapper->SetBondRadius(static_cast<float>(style == render::Style::Wire ? 0.05
                                                                             : BOND_RADIUS));
        // В объёмной модели связи всё равно целиком внутри шаров — считать их
        // незачем.
        mapper->SetRenderBonds(style != render::Style::SpaceFill);
        mapper->Modified();
        builtStyle = static_cast<int>(style);

        // Радиус выборки затенения. Должен доставать от одного атома до
        // соседнего, иначе в объёмной модели щели не потемнеют, но и не
        // больше того: на слишком большом радиусе выборка становится
        // разреженной и по шарам идёт зерно.
        double biggest = 0.3;
        for (const chem::Atom& atom : analysis->molecule.atoms) {
            biggest = std::max(biggest, view.atomRadius(atom.el));
        }
        frame.renderer()->SetSSAORadius(biggest * 0.8);
    }

    /**
     * Дуги углов и заливка секторов.
     *
     * Толщина трубки задаётся в ангстремах и растёт вместе с радиусом дуги:
     * в пикселях это даёт примерно ту же линию, что рисовал GDI+, на любой
     * молекуле и любом масштабе.
     */
    void buildArcs(const render::MoleculeRenderer& view) {
        const std::vector<render::AngleArc> arcs = view.angleArcs();

        vtkNew<vtkPoints> arcPoints;
        vtkNew<vtkCellArray> arcCells;
        vtkNew<vtkUnsignedCharArray> arcColors;
        arcColors->SetNumberOfComponents(4);
        arcColors->SetName("colors");

        vtkNew<vtkPoints> wedgePoints;
        vtkNew<vtkCellArray> wedgeCells;
        vtkNew<vtkUnsignedCharArray> wedgeColors;
        wedgeColors->SetNumberOfComponents(4);
        wedgeColors->SetName("colors");

        double widest = 0.02;
        for (const render::AngleArc& arc : arcs) {
            if (arc.linear || arc.points.size() < 2) continue;
            widest = std::max(widest, arc.radius * (arc.emphasised ? 0.055 : 0.038));

            const unsigned char r = static_cast<unsigned char>(arc.color.r);
            const unsigned char g = static_cast<unsigned char>(arc.color.g);
            const unsigned char b = static_cast<unsigned char>(arc.color.b);
            const unsigned char a = static_cast<unsigned char>(
                std::max(0.0, std::min(1.0, arc.alpha)) * 255);

            const vtkIdType first = arcPoints->GetNumberOfPoints();
            arcCells->InsertNextCell(static_cast<vtkIdType>(arc.points.size()));
            for (const chem::Vec3& p : arc.points) {
                arcCells->InsertCellPoint(arcPoints->InsertNextPoint(p.x, p.y, p.z));
                arcColors->InsertNextTuple4(r, g, b, a);
            }
            (void)first;

            // сектор: веер треугольников из вершины угла
            const vtkIdType apex = wedgePoints->InsertNextPoint(arc.center.x, arc.center.y,
                                                                arc.center.z);
            const unsigned char wedgeAlpha = static_cast<unsigned char>(
                std::max(0.0, std::min(1.0, arc.alpha)) * 0.10 * 255);
            wedgeColors->InsertNextTuple4(r, g, b, wedgeAlpha);
            std::vector<vtkIdType> rim;
            for (const chem::Vec3& p : arc.points) {
                rim.push_back(wedgePoints->InsertNextPoint(p.x, p.y, p.z));
                wedgeColors->InsertNextTuple4(r, g, b, wedgeAlpha);
            }
            for (std::size_t k = 1; k < rim.size(); k++) {
                wedgeCells->InsertNextCell(3);
                wedgeCells->InsertCellPoint(apex);
                wedgeCells->InsertCellPoint(rim[k - 1]);
                wedgeCells->InsertCellPoint(rim[k]);
            }
        }

        arcLines->SetPoints(arcPoints);
        arcLines->SetLines(arcCells);
        arcLines->GetPointData()->SetScalars(arcColors);
        arcTubes->SetRadius(widest);
        arcLines->Modified();
        arcActor->SetVisibility(arcPoints->GetNumberOfPoints() > 0);

        wedgeMesh->SetPoints(wedgePoints);
        wedgeMesh->SetPolys(wedgeCells);
        wedgeMesh->GetPointData()->SetScalars(wedgeColors);
        wedgeMesh->Modified();
        wedgeActor->SetVisibility(wedgePoints->GetNumberOfPoints() > 0);
    }

    /**
     * Изоповерхности гибридных орбиталей.
     *
     * Во второй версии лепестки рисовались кривыми Безье «на глаз». Здесь
     * ядро считает волновую функцию в узлах куба (chem/orbital.cpp), сетка
     * кладётся в vtkImageData, и vtkFlyingEdges3D проводит по ней две
     * поверхности: ψ = +iso и ψ = −iso. Разные знаки — разного цвета, как
     * орбитали и рисуют в учебниках, только посчитанные.
     *
     * Считать это каждый кадр незачем: от поворота камеры сетка не меняется.
     * Пересчёт идёт, только когда сменился набор орбиталей — другая молекула,
     * другой выбранный атом, включили или выключили показ.
     */
    void buildOrbitals(const render::MoleculeRenderer& view) {
        const std::vector<render::OrbitalSet> sets = view.orbitals();
        if (sameOrbitals(sets, builtOrbitals)) return;
        builtOrbitals = sets;

        for (const std::unique_ptr<OrbitalSurface>& old : orbitals) {
            frame.renderer()->RemoveActor(old->bond.actor);
            frame.renderer()->RemoveActor(old->pair.actor);
            frame.renderer()->RemoveActor(old->back.actor);
        }
        orbitals.clear();

        for (const render::OrbitalSet& set : sets) {
            const std::vector<chem::OrbitalGrid> grids =
                chem::orbitalGrids(set.kind, set.dirs, set.center, set.options);
            auto surface = std::make_unique<OrbitalSurface>();

            for (std::size_t k = 0; k < grids.size(); k++) {
                const chem::OrbitalGrid& grid = grids[k];
                const double plus = ORBITAL_ISO_PLUS * grid.highest();
                const double minus = ORBITAL_ISO_MINUS * grid.lowest();
                if (!(plus > 0)) continue;

                vtkNew<vtkFloatArray> scalars;
                scalars->SetName("psi");
                scalars->SetNumberOfComponents(1);
                scalars->SetNumberOfTuples(static_cast<vtkIdType>(grid.values.size()));
                for (std::size_t n = 0; n < grid.values.size(); n++) {
                    scalars->SetValue(static_cast<vtkIdType>(n), static_cast<float>(grid.values[n]));
                }
                auto field = vtkSmartPointer<vtkImageData>::New();
                field->SetDimensions(grid.size, grid.size, grid.size);
                field->SetOrigin(grid.center.x - grid.half, grid.center.y - grid.half,
                                 grid.center.z - grid.half);
                field->SetSpacing(grid.step(), grid.step(), grid.step());
                field->GetPointData()->SetScalars(scalars);
                surface->fields.push_back(field);

                const bool isBond = static_cast<int>(k) < set.bonding;
                addContour(*surface, field, plus, *(isBond ? surface->bond : surface->pair).mesh);
                addContour(*surface, field, minus, *surface->back.mesh);
                (isBond ? surface->bonds : surface->pairs)++;
            }
            if (surface->fields.empty()) continue;

            if (surface->bonds > 0) {
                setupSurface(surface->bond, render::theme::ORBITAL, ORBITAL_OPACITY_PLUS);
            }
            if (surface->pairs > 0) {
                setupSurface(surface->pair, render::theme::LONE_PAIR, ORBITAL_OPACITY_PLUS);
            }
            setupSurface(surface->back, render::theme::ORBITAL_MINUS, ORBITAL_OPACITY_MINUS);
            orbitals.push_back(std::move(surface));
        }
    }

    /** Изоповерхность одного уровня, приложенная к общей сшивке своего знака. */
    void addContour(OrbitalSurface& surface, vtkImageData* field, double level,
                    vtkAppendPolyData& into) {
        auto edges = vtkSmartPointer<vtkFlyingEdges3D>::New();
        edges->SetInputData(field);
        edges->SetNumberOfContours(1);
        edges->SetValue(0, level);
        // Нормали от градиента самой функции, а не от треугольников: с ними
        // поверхность выглядит гладкой при куда более грубой сетке.
        edges->ComputeNormalsOn();
        edges->ComputeGradientsOff();
        edges->ComputeScalarsOff();
        into.AddInputConnection(edges->GetOutputPort());
        surface.contours.push_back(edges);
    }

    /** Актёр одной сшивки: цвет, прозрачность, свет. */
    void setupSurface(SurfacePart& part, const render::Rgb& color, double opacity) {
        vtkActor& actor = *part.actor;
        part.mapper->SetInputConnection(part.mesh->GetOutputPort());
        part.mapper->ScalarVisibilityOff();
        actor.SetMapper(part.mapper);
        actor.GetProperty()->SetColor(color.r / 255.0, color.g / 255.0, color.b / 255.0);
        actor.GetProperty()->SetOpacity(opacity);
        // Изнутри лепестка видна его обратная сторона, и без освещения
        // обратных граней она выглядит чёрной дырой.
        actor.GetProperty()->BackfaceCullingOff();
        actor.GetProperty()->SetAmbient(0.25);
        actor.GetProperty()->SetDiffuse(0.75);
        actor.GetProperty()->SetSpecular(0.2);
        frame.renderer()->AddActor(&actor);
    }

    /** Совпадает ли набор орбиталей с тем, под который уже построены сетки. */
    static bool sameOrbitals(const std::vector<render::OrbitalSet>& a,
                             const std::vector<render::OrbitalSet>& b) {
        if (a.size() != b.size()) return false;
        for (std::size_t i = 0; i < a.size(); i++) {
            if (a[i].kind != b[i].kind || a[i].dirs.size() != b[i].dirs.size()) return false;
            if (a[i].bonding != b[i].bonding) return false;
            if (a[i].options.radial != b[i].options.radial) return false;
            if (std::abs(a[i].options.zeff - b[i].options.zeff) > 1e-9) return false;
            if (std::abs(a[i].options.half - b[i].options.half) > 1e-9) return false;
            if (chem::dist(a[i].center, b[i].center) > 1e-9) return false;
            for (std::size_t k = 0; k < a[i].dirs.size(); k++) {
                if (chem::dist(a[i].dirs[k], b[i].dirs[k]) > 1e-9) return false;
            }
        }
        return true;
    }

    /**
     * Камера VTK по параметрам собственного рендера.
     *
     * Собственная проекция устроена так: c = R·p — точка в системе камеры,
     * глубина z = c.z + d, экранные координаты (cx + c.x·f/z, cy − c.y·f/z).
     * Значит глаз стоит в системе камеры в точке (0, 0, −d), смотрит в +z,
     * верх у него +y. Переводим это в мировые координаты обратным поворотом.
     *
     * Зеркальность (мировой +x у собственного рендера уходит вправо, у
     * правосторонней камеры — влево) снимается разворотом кадра в Frame.
     */
    void applyCamera(const render::CameraParams& cam) {
        const chem::Mat3 inverse = chem::transpose(cam.rotation);
        const chem::Vec3 eye = chem::applyMat(inverse, chem::v3(0, 0, -cam.distance));
        const chem::Vec3 up = chem::applyMat(inverse, chem::v3(0, 1, 0));

        vtkCamera* camera = frame.renderer()->GetActiveCamera();
        camera->SetPosition(eye.x, eye.y, eye.z);
        camera->SetFocalPoint(0, 0, 0);
        camera->SetViewUp(up.x, up.y, up.z);
        // Угол обзора по вертикали: половина высоты кадра делится на фокусное.
        const double halfHeight = std::max(1.0, cam.viewport.h) / 2;
        const double angle = 2 * std::atan(halfHeight / std::max(1.0, cam.focal));
        camera->SetViewAngle(angle * 180.0 / chem::PI);
        // Плоскости отсечения с запасом: иначе ближние шары срежет.
        // имена near/far заняты макросами из windows.h
        const double nearPlane = std::max(0.05, cam.distance - cam.sceneRadius - 1.0);
        camera->SetClippingRange(nearPlane, cam.distance + cam.sceneRadius + 1.0);
    }
};

Scene::Scene() : impl(std::make_unique<Impl>()) {}
Scene::~Scene() = default;

void Scene::setAnalysis(const chem::Analysis* analysis) {
    impl->analysis = analysis;
    impl->builtMolecule = false;
    impl->builtOrbitals.clear();
}

void Scene::draw(HDC target, const render::MoleculeRenderer& view) {
    const render::CameraParams cam = view.camera();
    const int width = static_cast<int>(cam.viewport.w + 0.5);
    const int height = static_cast<int>(cam.viewport.h + 0.5);
    if (width <= 0 || height <= 0) return;

    if (!impl->builtMolecule) impl->buildMolecule();
    if (impl->analysis == nullptr) return;
    if (impl->builtStyle != static_cast<int>(view.options().style)) impl->buildRadii(view);

    LARGE_INTEGER t0, t1, t2;
    QueryPerformanceCounter(&t0);
    impl->buildArcs(view);
    impl->buildOrbitals(view);
    impl->frame.resize(width, height);
    impl->applyCamera(cam);
    impl->frame.render();
    QueryPerformanceCounter(&t1);

    impl->frame.blendTo(target, static_cast<int>(cam.viewport.x + 0.5),
                        static_cast<int>(cam.viewport.y + 0.5));
    QueryPerformanceCounter(&t2);

    impl->renderMs = milliseconds(t0, t1);
    impl->transferMs = milliseconds(t1, t2);

}

std::wstring Scene::saveImage(const render::MoleculeRenderer& view) {
    if (impl->analysis == nullptr) return L"";
    const render::CameraParams cam = view.camera();
    const int width = static_cast<int>(cam.viewport.w + 0.5);
    if (width <= 0) return L"";

    // Снимок должен годиться для печати: не меньше 3200 пикселей по ширине.
    // Больше шестикратного увеличения не берём — смысла нет, а памяти
    // на кадр уходит уже под сотню мегабайт.
    const int scale = std::max(2, std::min(6, (3200 + width - 1) / width));

    wchar_t folder[MAX_PATH] = {};
    if (SHGetFolderPathW(nullptr, CSIDL_MYPICTURES, nullptr, 0, folder) != S_OK) return L"";

    // Имя латиницей: путь уходит в VTK, и с кириллицей в имени файла она
    // ведёт себя по-разному в зависимости от сборки.
    std::string id = "molecule";
    if (impl->analysis->molecule.meta != nullptr && !impl->analysis->molecule.meta->id.empty()) {
        id = impl->analysis->molecule.meta->id;
    }
    const std::wstring path = std::wstring(folder) + L"\\molview-" + render::toWide(id) + L".png";

    const render::Rgb back = render::theme::BACKGROUND;
    if (!impl->frame.saveImage(render::toUtf8(path), scale, static_cast<unsigned char>(back.r),
                               static_cast<unsigned char>(back.g),
                               static_cast<unsigned char>(back.b))) {
        return L"";
    }
    return path;
}

double Scene::lastRenderMs() const { return impl->renderMs; }
double Scene::lastTransferMs() const { return impl->transferMs; }

}  // namespace scene
