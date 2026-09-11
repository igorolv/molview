#include "scene/scene.h"

#include <vtkActor.h>
#include <vtkDataSetAttributes.h>
#include <vtkCamera.h>
#include <vtkBlueObeliskData.h>
#include <vtkFloatArray.h>
#include <vtkLightKit.h>
#include <vtkMolecule.h>
#include <vtkMoleculeMapper.h>
#include <vtkNew.h>
#include <vtkPeriodicTable.h>
#include <vtkProperty.h>

#include <shlobj.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "chem/periodic.h"
#include "chem/vec.h"
#include "render/palette.h"
#include "render/draw.h"
#include "scene/frame.h"

namespace scene {

namespace {

constexpr double BOND_RADIUS = 0.115;  // как в собственном рендере

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
