#include "render/renderer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "chem/periodic.h"

namespace render {

using chem::add;
using chem::cross;
using chem::dot;
using chem::len;
using chem::Mat3;
using chem::mul;
using chem::norm;
using chem::sub;
using chem::v3;
using chem::Vec3;

namespace {

constexpr double ATOM_RADIUS_FACTOR = 0.34;
constexpr double BOND_RADIUS = 0.115;
constexpr double ENTRY_DURATION = 750;

double clamp01(double v) { return std::max(0.0, std::min(1.0, v)); }

/** Появление атома при загрузке молекулы: волной от центра, с лёгким «перелётом». */
double atomAppearance(double progress, std::size_t index, std::size_t total) {
    if (progress >= 1) return 1;
    const double delay = (static_cast<double>(index) / std::max<std::size_t>(1, total)) * 0.45;
    const double t = clamp01((progress - delay) / 0.55);
    const double eased = 1 - std::pow(1 - t, 3);
    return eased * (1 + 0.12 * std::sin(t * chem::PI));
}

/** Интерполяция направлений по дуге большого круга. */
Vec3 slerp(const Vec3& a, const Vec3& b, double t) {
    const double cosine = std::max(-1.0, std::min(1.0, dot(a, b)));
    const double omega = std::acos(cosine);
    if (omega < 1e-4) return a;
    if (chem::PI - omega < 1e-4) {
        // векторы противоположны — идём через произвольную перпендикулярную ось
        const Vec3 axis = norm(cross(a, std::fabs(a.x) < 0.9 ? v3(1, 0, 0) : v3(0, 1, 0)));
        return chem::applyMat(chem::rotationAxis(axis, omega * t), a);
    }
    const double s = std::sin(omega);
    return add(mul(a, std::sin((1 - t) * omega) / s), mul(b, std::sin(t * omega) / s));
}

std::wstring degrees(double value) {
    wchar_t buffer[32];
    std::swprintf(buffer, 32, L"%.1f°", value);
    return buffer;
}

}  // namespace

// ---------------------------------------------------------------------------
// Состояние
// ---------------------------------------------------------------------------

void MoleculeRenderer::setAnalysis(const chem::Analysis* value, double timeMs, bool animate) {
    analysis = value;
    selectedAtom = -1;
    hoveredAtom = -1;
    // проекции прошлой молекулы больше не действительны
    projected.clear();

    rotation = chem::IDENTITY;
    zoom = 1;
    entryStart = animate ? timeMs : 0;
}

void MoleculeRenderer::resetView() {
    rotation = chem::IDENTITY;
    zoom = 1;
}

void MoleculeRenderer::tick(double dtMs) {
    if (!opts.autoRotate || dragging) return;
    const double dt = std::min(60.0, dtMs);
    rotation = chem::matMul(chem::rotationAxis(v3(0, 1, 0), dt * 0.00022), rotation);
}

void MoleculeRenderer::beginDrag(double x, double y) {
    dragging = true;
    pointerMoved = false;
    lastPointerX = x;
    lastPointerY = y;
}

void MoleculeRenderer::dragTo(double x, double y) {
    if (!dragging) return;
    const double dx = x - lastPointerX;
    const double dy = y - lastPointerY;
    if (std::fabs(dx) + std::fabs(dy) > 2) pointerMoved = true;
    lastPointerX = x;
    lastPointerY = y;

    const double speed = 0.0085;
    const Mat3 spin = chem::matMul(chem::rotationAxis(v3(0, 1, 0), dx * speed),
                                   chem::rotationAxis(v3(1, 0, 0), dy * speed));
    rotation = chem::matMul(spin, rotation);
}

bool MoleculeRenderer::endDrag() {
    const bool wasClick = dragging && !pointerMoved;
    dragging = false;
    return wasClick;
}

void MoleculeRenderer::wheel(double delta) {
    zoom = std::max(0.35, std::min(4.0, zoom * std::exp(-delta * 0.0012)));
}

int MoleculeRenderer::hitTest(double px, double py) const {
    if (analysis == nullptr) return -1;
    int found = -1;
    double bestZ = 1e300;

    for (std::size_t i = 0; i < projected.size() && i < analysis->molecule.atoms.size(); i++) {
        const Projected& p = projected[i];
        if (!p.visible) continue;
        const double r = std::max(9.0, atomRadius(analysis->molecule.atoms[i].el) * p.scale);
        const double ddx = px - p.x;
        const double ddy = py - p.y;
        if (ddx * ddx + ddy * ddy <= r * r && p.z < bestZ) {
            bestZ = p.z;
            found = static_cast<int>(i);
        }
    }
    return found;
}

// ---------------------------------------------------------------------------
// Проекция
// ---------------------------------------------------------------------------

double MoleculeRenderer::atomRadius(const std::string& el) const {
    const chem::ElementInfo& info = chem::element(el);
    if (opts.style == Style::SpaceFill) return info.vdw * 0.62;
    if (opts.style == Style::Wire) return std::max(0.12, info.r * 0.18);
    return std::max(0.26, info.r * ATOM_RADIUS_FACTOR + 0.12);
}

double MoleculeRenderer::sceneRadius() const {
    if (analysis == nullptr) return 4;
    double r = 0.8;
    for (const chem::Atom& atom : analysis->molecule.atoms) {
        r = std::max(r, len(atom.pos) + atomRadius(atom.el));
    }
    if (opts.showOrbitals) r += 0.5;
    else if (opts.showLonePairs) r += 0.35;
    return r + 0.3;
}

double MoleculeRenderer::cameraDistance() const { return (sceneRadius() * 2.85) / zoom; }

double MoleculeRenderer::focal() const { return std::min(viewport.w, viewport.h) * 1.25; }

MoleculeRenderer::Projected MoleculeRenderer::project(const Vec3& p) const {
    const Vec3 camera = chem::applyMat(rotation, p);
    const double z = camera.z + cameraDistance();
    if (z < 0.15) return Projected{0, 0, z, 0, false};
    const double s = focal() / z;
    return Projected{
        viewport.x + viewport.w / 2 + camera.x * s,
        viewport.y + viewport.h / 2 - camera.y * s,
        z, s, true,
    };
}

Rgb MoleculeRenderer::fog(const Rgb& color, double z) const {
    const double radius = sceneRadius();
    // near и far заняты макросами из windows.h, поэтому имена другие
    const double nearZ = cameraDistance() - radius;
    const double farZ = cameraDistance() + radius;
    const double t = clamp01((z - nearZ) / std::max(0.001, farZ - nearZ));
    return mix(color, theme::FOG, t * 0.55);
}

// ---------------------------------------------------------------------------
// Кадр целиком
// ---------------------------------------------------------------------------

CameraParams MoleculeRenderer::camera() const {
    return CameraParams{rotation, cameraDistance(), focal(), sceneRadius(), viewport};
}

void MoleculeRenderer::render(Gdiplus::Graphics& g, Fonts& fonts, double timeMs) {
    renderBackground(g, timeMs);
    if (analysis == nullptr) return;

    const double progress = entryStart == 0 ? 1.0 : clamp01((timeMs - entryStart) / ENTRY_DURATION);

    std::vector<Primitive> primitives;
    collectBonds(primitives, progress);
    collectAtoms(primitives, progress, timeMs);

    // алгоритм художника: сначала дальние
    std::stable_sort(primitives.begin(), primitives.end(),
                     [](const Primitive& a, const Primitive& b) { return a.depth > b.depth; });
    for (const Primitive& p : primitives) p.draw(g);

    renderOverlays(g, fonts, timeMs);
}

void MoleculeRenderer::renderBackground(Gdiplus::Graphics& g, double timeMs) {
    frameProgress = entryStart == 0 ? 1.0 : clamp01((timeMs - entryStart) / ENTRY_DURATION);
    drawBackground(g);
    projected.clear();
    updateOcclusion();
    if (analysis == nullptr) return;
    // Проекции нужны и наложениям, и проверке попадания мыши, поэтому
    // считаются до сцены, а не вместе с шарами.
    for (const chem::Atom& a : analysis->molecule.atoms) projected.push_back(project(a.pos));
}

void MoleculeRenderer::updateOcclusion() {
    // Обращение поворота переводит оси камеры в мировые координаты. Глаз стоит
    // в системе камеры в точке (0, 0, −d) — см. комментарий к project().
    const Mat3 inverse = chem::transpose(rotation);
    eyePos = chem::applyMat(inverse, v3(0, 0, -cameraDistance()));
    camRight = chem::applyMat(inverse, v3(1, 0, 0));
    camUp = chem::applyMat(inverse, v3(0, 1, 0));

    occluder.clear();
    if (analysis == nullptr) return;
    for (const chem::Atom& atom : analysis->molecule.atoms) {
        occluder.add(atom.pos, atomRadius(atom.el));
    }
}

void MoleculeRenderer::renderOverlays(Gdiplus::Graphics& g, Fonts& fonts, double timeMs) {
    if (analysis == nullptr) return;
    const double progress = entryStart == 0 ? 1.0 : clamp01((timeMs - entryStart) / ENTRY_DURATION);

    // Неподелённые пары и орбитали пока рисуются здесь, а значит ложатся
    // ПОВЕРХ шаров, даже когда должны быть за ними. Сортировка по глубине
    // вместе с атомами вернётся на этапе 4, когда они переедут в VTK.
    std::vector<Primitive> primitives;
    if (opts.showLonePairs) collectLonePairs(primitives, progress);
    if (opts.showOrbitals) collectOrbitals(primitives, progress);
    std::stable_sort(primitives.begin(), primitives.end(),
                     [](const Primitive& a, const Primitive& b) { return a.depth > b.depth; });
    for (const Primitive& p : primitives) p.draw(g);

    if (opts.showDipole) drawDipole(g, fonts, progress);
    if (opts.showAngles != AngleMode::None) drawAngles(g, fonts, progress);
    if (opts.showLabels) drawLabels(g, fonts, progress);
}

void MoleculeRenderer::drawBackground(Gdiplus::Graphics& g) const {
    const double cx = viewport.x + viewport.w / 2;
    const double cy = viewport.y + viewport.h * 0.46;
    const double radius = std::max(viewport.w, viewport.h) * 0.78;

    // прямоугольник закрашивается целиком, потом поверх ложится сияние
    Gdiplus::SolidBrush base(toColor(theme::BACKGROUND));
    g.FillRectangle(&base, static_cast<Gdiplus::REAL>(viewport.x), static_cast<Gdiplus::REAL>(viewport.y),
                    static_cast<Gdiplus::REAL>(viewport.w), static_cast<Gdiplus::REAL>(viewport.h));

    Gdiplus::Region saved;
    g.GetClip(&saved);
    g.SetClip(Gdiplus::RectF(static_cast<Gdiplus::REAL>(viewport.x), static_cast<Gdiplus::REAL>(viewport.y),
                             static_cast<Gdiplus::REAL>(viewport.w), static_cast<Gdiplus::REAL>(viewport.h)),
              Gdiplus::CombineModeIntersect);

    fillRadial(g, cx, cy, radius, radius, cx, cy, {
        {0.0, toColor(theme::BACKGROUND_GLOW)},
        {0.55, toColor(theme::BACKGROUND_MID)},
        {1.0, toColor(theme::BACKGROUND)},
    });

    // едва заметная сетка «лабораторного стола»
    Gdiplus::Pen grid(Gdiplus::Color(12, 120, 170, 255), 1.0f);
    const double step = 46;
    for (double x = viewport.x + std::fmod(viewport.w / 2, step); x < viewport.right(); x += step) {
        g.DrawLine(&grid, static_cast<Gdiplus::REAL>(x), static_cast<Gdiplus::REAL>(viewport.y),
                   static_cast<Gdiplus::REAL>(x), static_cast<Gdiplus::REAL>(viewport.bottom()));
    }
    for (double y = viewport.y + std::fmod(viewport.h / 2, step); y < viewport.bottom(); y += step) {
        g.DrawLine(&grid, static_cast<Gdiplus::REAL>(viewport.x), static_cast<Gdiplus::REAL>(y),
                   static_cast<Gdiplus::REAL>(viewport.right()), static_cast<Gdiplus::REAL>(y));
    }

    g.SetClip(&saved, Gdiplus::CombineModeReplace);
}

// ---------------------------------------------------------------------------
// Атомы
// ---------------------------------------------------------------------------

void MoleculeRenderer::collectAtoms(std::vector<Primitive>& out, double progress, double timeMs) const {
    const chem::Molecule& mol = analysis->molecule;

    for (std::size_t index = 0; index < mol.atoms.size(); index++) {
        const Projected p = projected[index];
        if (!p.visible) continue;

        const double appear = atomAppearance(progress, index, mol.atoms.size());
        if (appear <= 0.001) continue;

        const double baseRadius = atomRadius(mol.atoms[index].el) * p.scale * appear;
        const Rgb color = fog(hexToRgb(chem::element(mol.atoms[index].el).color), p.z);
        const bool isHovered = hoveredAtom == static_cast<int>(index);
        const bool isSelected = selectedAtom == static_cast<int>(index);

        Primitive prim;
        prim.depth = p.z;
        prim.draw = [p, baseRadius, color, isHovered, isSelected, timeMs](Gdiplus::Graphics& g) {
            const double r = baseRadius * (isHovered ? 1.1 : 1.0);
            if (r < 0.4) return;

            if (isSelected || isHovered) {
                const double pulse = isSelected ? 0.5 + 0.5 * std::sin(timeMs * 0.005) : 0.55;
                fillRadial(g, p.x, p.y, r * 2.3, r * 2.3, p.x, p.y, {
                    {0.0, Gdiplus::Color(static_cast<BYTE>(0.36 * pulse * 255), 94, 234, 212)},
                    {1.0, Gdiplus::Color(0, 94, 234, 212)},
                });
            }

            // тело шара: блик смещён к верхнему левому краю
            fillRadial(g, p.x, p.y, r, r, p.x - r * 0.36, p.y - r * 0.42, {
                {0.0, toColor(lighten(color, 0.62))},
                {0.42, toColor(lighten(color, 0.1))},
                {1.0, toColor(darken(color, 0.52))},
            });

            // ободок — отделяет атом от соседей
            const double lineWidth = std::max(0.6, r * 0.07);
            const Gdiplus::Color edge = isSelected
                ? Gdiplus::Color(242, 94, 234, 212)
                : Gdiplus::Color(static_cast<BYTE>((0.16 + (isHovered ? 0.3 : 0.0)) * 255), 255, 255, 255);
            Gdiplus::Pen pen(edge, static_cast<Gdiplus::REAL>(lineWidth));
            g.DrawEllipse(&pen, static_cast<Gdiplus::REAL>(p.x - r), static_cast<Gdiplus::REAL>(p.y - r),
                          static_cast<Gdiplus::REAL>(r * 2), static_cast<Gdiplus::REAL>(r * 2));

            // блик
            const double hx = p.x - r * 0.34;
            const double hy = p.y - r * 0.42;
            fillRadial(g, hx, hy, r * 0.42, r * 0.42, hx, hy, {
                {0.0, Gdiplus::Color(140, 255, 255, 255)},
                {1.0, Gdiplus::Color(0, 255, 255, 255)},
            });
        };
        out.push_back(std::move(prim));
    }
}

// ---------------------------------------------------------------------------
// Связи
// ---------------------------------------------------------------------------

void MoleculeRenderer::collectBonds(std::vector<Primitive>& out, double progress) const {
    if (opts.style == Style::SpaceFill) return;
    const chem::Molecule& mol = analysis->molecule;

    for (const chem::Bond& bond : mol.bonds) {
        const Projected pa = projected[bond.a];
        const Projected pb = projected[bond.b];
        if (!pa.visible || !pb.visible) continue;

        const double grow = clamp01((progress - 0.12) / 0.55);
        if (grow <= 0.001) continue;

        const Projected mid = project(mul(add(mol.atoms[bond.a].pos, mol.atoms[bond.b].pos), 0.5));

        const Rgb colorA = fog(hexToRgb(chem::element(mol.atoms[bond.a].el).color), (pa.z + mid.z) / 2);
        const Rgb colorB = fog(hexToRgb(chem::element(mol.atoms[bond.b].el).color), (pb.z + mid.z) / 2);

        // перпендикуляр к оси связи в экранных координатах
        const double dx = pb.x - pa.x;
        const double dy = pb.y - pa.y;
        const double l = std::max(1e-6, std::hypot(dx, dy));
        const double perpX = -dy / l;
        const double perpY = dx / l;

        const int count = bond.order >= 2.9 ? 3 : (bond.order >= 1.4 ? 2 : 1);
        const double separation = BOND_RADIUS * 2.05;
        const bool dashedBond = bond.aromatic || bond.delocalized;

        for (int k = 0; k < count; k++) {
            const double offset = count == 1 ? 0 : (k - (count - 1) / 2.0) * separation;
            // вторая линия у ароматической связи рисуется пунктиром — так на схеме
            // видно, что порядок связи дробный
            const bool isDashedLine = dashedBond && k == count - 1;
            const double widthScale = isDashedLine ? 0.55 : 1.0;

            const auto half = [&](const Projected& from, const Rgb& color) {
                const double w1 = BOND_RADIUS * from.scale * widthScale;
                const double w2 = BOND_RADIUS * mid.scale * widthScale;
                const double o1 = offset * from.scale;
                const double o2 = offset * mid.scale;

                const double ax = from.x + perpX * o1;
                const double ay = from.y + perpY * o1;
                const double bx = ax + (mid.x + perpX * o2 - ax) * grow;
                const double by = ay + (mid.y + perpY * o2 - ay) * grow;

                Primitive prim;
                prim.depth = (from.z + mid.z) / 2;
                prim.draw = [=](Gdiplus::Graphics& g) {
                    if (isDashedLine) {
                        Gdiplus::Pen pen(toColor(lighten(color, 0.18), 0.8),
                                         static_cast<Gdiplus::REAL>(std::max(0.4, w1 * 1.7)));
                        const Gdiplus::REAL dashes[] = {2.2f, 1.8f};
                        pen.SetDashPattern(dashes, 2);
                        pen.SetStartCap(Gdiplus::LineCapRound);
                        pen.SetEndCap(Gdiplus::LineCapRound);
                        g.DrawLine(&pen, static_cast<Gdiplus::REAL>(ax), static_cast<Gdiplus::REAL>(ay),
                                   static_cast<Gdiplus::REAL>(bx), static_cast<Gdiplus::REAL>(by));
                        return;
                    }

                    // трапеция, имитирующая цилиндр: у ближнего конца шире
                    Gdiplus::GraphicsPath path;
                    Gdiplus::PointF points[4] = {
                        {static_cast<Gdiplus::REAL>(ax + perpX * w1), static_cast<Gdiplus::REAL>(ay + perpY * w1)},
                        {static_cast<Gdiplus::REAL>(bx + perpX * w2), static_cast<Gdiplus::REAL>(by + perpY * w2)},
                        {static_cast<Gdiplus::REAL>(bx - perpX * w2), static_cast<Gdiplus::REAL>(by - perpY * w2)},
                        {static_cast<Gdiplus::REAL>(ax - perpX * w1), static_cast<Gdiplus::REAL>(ay - perpY * w1)},
                    };
                    path.AddPolygon(points, 4);

                    const double gx = (ax + bx) / 2;
                    const double gy = (ay + by) / 2;
                    const double gw = std::max(0.3, (w1 + w2) / 2);
                    fillPathLinear(g, path,
                                   gx + perpX * gw, gy + perpY * gw,
                                   gx - perpX * gw, gy - perpY * gw, {
                                       {0.0, toColor(darken(color, 0.55))},
                                       {0.38, toColor(lighten(color, 0.3))},
                                       {1.0, toColor(darken(color, 0.55))},
                                   });
                };
                out.push_back(std::move(prim));
            };

            half(pa, colorA);
            half(pb, colorB);
        }
    }
}

// ---------------------------------------------------------------------------
// Неподелённые пары
// ---------------------------------------------------------------------------

void MoleculeRenderer::collectLonePairs(std::vector<Primitive>& out, double progress) const {
    if (progress < 0.5) return;
    const double fade = std::min(1.0, (progress - 0.5) / 0.35);

    // Пары концевых атомов (три пары на каждом атоме хлора и т.п.) сильно
    // засоряют картинку, а форму молекулы определяют пары ЦЕНТРАЛЬНЫХ атомов.
    // Поэтому по умолчанию показываем только их — плюс пары выбранного атома.
    bool hasCentral = false;
    for (const chem::AtomAnalysis& a : analysis->atoms) {
        if (a.isCentral && !a.lonePairDirs.empty()) { hasCentral = true; break; }
    }

    for (const chem::AtomAnalysis& info : analysis->atoms) {
        if (info.lonePairDirs.empty()) continue;
        if (hasCentral && !info.isCentral && selectedAtom != info.id) continue;

        const chem::Atom& atom = analysis->molecule.atoms[info.id];
        const double reach = atomRadius(atom.el) + 0.42;

        for (const Vec3& dir : info.lonePairDirs) {
            const Vec3 cloud = add(atom.pos, mul(dir, reach));
            const Projected p = project(cloud);
            if (!p.visible) continue;
            // Облако — мягкое пятно диаметром в треть ангстрема, то есть почти
            // точка: одной проверки хватает. Закрытое наполовину исчезает
            // целиком, а не срезается по краю шара, но при таких размерах это
            // заметно только если специально искать.
            if (!occluder.visible(eyePos, cloud)) continue;

            const Projected base = project(atom.pos);
            const double axisX = p.x - base.x;
            const double axisY = p.y - base.y;
            const double axisLen = std::hypot(axisX, axisY);
            const double angle = axisLen > 1 ? std::atan2(axisY, axisX) : 0;

            // сплющиваем облако, когда пара смотрит «на зрителя»
            const double foreshorten = std::max(0.32, std::min(1.0, axisLen / std::max(1e-6, reach * p.scale)));
            const double rx = 0.30 * p.scale;
            const double ry = rx * (0.55 + 0.45 * foreshorten);
            const bool dimmed = hoveredAtom != -1 && hoveredAtom != info.id;
            const double alpha = (dimmed ? 0.28 : 0.62) * fade;

            Primitive prim;
            prim.depth = p.z;
            prim.draw = [p, rx, ry, angle, alpha](Gdiplus::Graphics& g) {
                if (rx < 0.5) return;
                const Gdiplus::GraphicsState state = g.Save();
                g.TranslateTransform(static_cast<Gdiplus::REAL>(p.x), static_cast<Gdiplus::REAL>(p.y));
                g.RotateTransform(static_cast<Gdiplus::REAL>(angle * 180.0 / chem::PI));

                fillRadial(g, 0, 0, rx, ry, 0, 0, {
                    {0.0, Gdiplus::Color(static_cast<BYTE>(std::min(1.0, alpha) * 255), 147, 214, 255)},
                    {0.6, Gdiplus::Color(static_cast<BYTE>(std::min(1.0, alpha * 0.5) * 255), 125, 211, 252)},
                    {1.0, Gdiplus::Color(0, 125, 211, 252)},
                });

                // два электрона
                const double dotR = std::max(1.1, rx * 0.15);
                const double gap = ry * 0.5;
                Gdiplus::SolidBrush dots(Gdiplus::Color(
                    static_cast<BYTE>(std::min(1.0, alpha * 1.5) * 255), 224, 245, 255));
                for (int sign = -1; sign <= 1; sign += 2) {
                    g.FillEllipse(&dots, static_cast<Gdiplus::REAL>(-dotR),
                                  static_cast<Gdiplus::REAL>(sign * gap - dotR),
                                  static_cast<Gdiplus::REAL>(dotR * 2),
                                  static_cast<Gdiplus::REAL>(dotR * 2));
                }
                g.Restore(state);
            };
            out.push_back(std::move(prim));
        }
    }
}

// ---------------------------------------------------------------------------
// Гибридные орбитали
// ---------------------------------------------------------------------------

void MoleculeRenderer::collectOrbitals(std::vector<Primitive>& out, double progress) const {
    if (progress < 0.55) return;
    const double fade = std::min(1.0, (progress - 0.55) / 0.35);

    for (const chem::AtomAnalysis& info : analysis->atoms) {
        if (info.steric < 2 || info.el == "H") continue;
        if (selectedAtom != -1 && selectedAtom != info.id) continue;
        if (selectedAtom == -1 && !info.isCentral) continue;

        const chem::Atom& atom = analysis->molecule.atoms[info.id];
        const double lobeLength = 1.15;

        for (std::size_t k = 0; k < info.orbitalDirs.size(); k++) {
            const bool isLonePair = static_cast<int>(k) >= info.sigma;
            const Projected tip = project(add(atom.pos, mul(info.orbitalDirs[k], lobeLength)));
            const Projected base = project(atom.pos);
            if (!tip.visible || !base.visible) continue;

            const double dx = tip.x - base.x;
            const double dy = tip.y - base.y;
            const double l = std::hypot(dx, dy);
            const double ux = l > 0.001 ? dx / l : 1;
            const double uy = l > 0.001 ? dy / l : 0;
            const double px = -uy;
            const double py = ux;
            const double width = 0.34 * base.scale;
            const Rgb color = isLonePair ? theme::LONE_PAIR : theme::ORBITAL;
            const double alpha = 0.30 * fade;

            Primitive prim;
            prim.depth = (tip.z + base.z) / 2 + 0.01;
            prim.draw = [=](Gdiplus::Graphics& g) {
                const auto f = [](double v) { return static_cast<Gdiplus::REAL>(v); };

                // передний лепесток — две кривые Безье, сходящиеся в остриё
                Gdiplus::GraphicsPath path;
                path.AddBezier(f(base.x), f(base.y),
                               f(base.x + ux * l * 0.28 + px * width), f(base.y + uy * l * 0.28 + py * width),
                               f(base.x + ux * l * 0.82 + px * width * 0.85), f(base.y + uy * l * 0.82 + py * width * 0.85),
                               f(tip.x), f(tip.y));
                path.AddBezier(f(tip.x), f(tip.y),
                               f(base.x + ux * l * 0.82 - px * width * 0.85), f(base.y + uy * l * 0.82 - py * width * 0.85),
                               f(base.x + ux * l * 0.28 - px * width), f(base.y + uy * l * 0.28 - py * width),
                               f(base.x), f(base.y));
                path.CloseFigure();

                fillPathLinear(g, path, base.x, base.y, tip.x, tip.y, {
                    {0.0, toColor(color, alpha * 0.35)},
                    {0.55, toColor(color, alpha)},
                    {1.0, toColor(color, alpha * 0.15)},
                });
                Gdiplus::Pen pen(toColor(color, std::min(1.0, alpha * 1.3)), 1.0f);
                g.DrawPath(&pen, &path);

                // малый «задний» лепесток гибридной орбитали
                const double backLen = l * 0.26;
                Gdiplus::GraphicsPath back;
                back.AddBezier(f(base.x), f(base.y),
                               f(base.x - ux * backLen * 0.5 + px * width * 0.5), f(base.y - uy * backLen * 0.5 + py * width * 0.5),
                               f(base.x - ux * backLen + px * width * 0.2), f(base.y - uy * backLen + py * width * 0.2),
                               f(base.x - ux * backLen), f(base.y - uy * backLen));
                back.AddBezier(f(base.x - ux * backLen), f(base.y - uy * backLen),
                               f(base.x - ux * backLen - px * width * 0.2), f(base.y - uy * backLen - py * width * 0.2),
                               f(base.x - ux * backLen * 0.5 - px * width * 0.5), f(base.y - uy * backLen * 0.5 - py * width * 0.5),
                               f(base.x), f(base.y));
                back.CloseFigure();
                Gdiplus::SolidBrush brush(toColor(color, alpha * 0.35));
                g.FillPath(&brush, &back);
            };
            out.push_back(std::move(prim));
        }
    }
}

// ---------------------------------------------------------------------------
// Валентные углы
// ---------------------------------------------------------------------------

bool MoleculeRenderer::labelCollides(const RectD& r) const {
    const double gap = 3;
    for (const RectD& o : labelRects) {
        if (r.x < o.x + o.w + gap && r.x + r.w + gap > o.x
            && r.y < o.y + o.h + gap && r.y + r.h + gap > o.y) {
            return true;
        }
    }
    return false;
}

std::vector<AngleArc> MoleculeRenderer::angleArcs() const {
    std::vector<AngleArc> arcs;
    if (analysis == nullptr || opts.showAngles == AngleMode::None) return arcs;
    if (frameProgress < 0.6) return arcs;
    const double fade = std::min(1.0, (frameProgress - 0.6) / 0.3);
    const int focus = selectedAtom != -1 ? selectedAtom : hoveredAtom;

    const chem::Molecule& mol = analysis->molecule;
    std::vector<const chem::AngleRecord*> records;
    for (const chem::AngleRecord& r : analysis->angles) records.push_back(&r);

    // В режиме «все» на большой молекуле углы с участием водорода дают десятки
    // почти одинаковых дуг и превращают картинку в кашу. Оставляем скелет.
    if (opts.showAngles == AngleMode::All && mol.atoms.size() > 8) {
        std::vector<const chem::AngleRecord*> heavy;
        for (const chem::AngleRecord* r : records) {
            if (mol.atoms[r->i].el != "H" && mol.atoms[r->j].el != "H") heavy.push_back(r);
        }
        records = heavy;
    }

    if (opts.showAngles == AngleMode::Selected) {
        int center = focus;
        if (center == -1) {
            // ничего не выбрано — показываем углы самого «главного» центра
            const chem::AtomAnalysis* best = nullptr;
            for (const chem::AtomAnalysis& a : analysis->atoms) {
                if (!a.isCentral || a.el == "H") continue;
                if (best == nullptr || a.sigma > best->sigma) best = &a;
            }
            if (best == nullptr) return arcs;
            center = best->id;
        }
        std::vector<const chem::AngleRecord*> filtered;
        for (const chem::AngleRecord* r : records) {
            if (r->center == center) filtered.push_back(r);
        }
        records = filtered;
    }
    if (records.empty()) return arcs;

    for (const chem::AngleRecord* record : records) {
        const Vec3 c = mol.atoms[record->center].pos;
        AngleArc arc;
        arc.record = record;
        arc.center = c;
        arc.from = norm(sub(mol.atoms[record->i].pos, c));
        arc.to = norm(sub(mol.atoms[record->j].pos, c));
        arc.radius = 0.46 * std::min(len(sub(mol.atoms[record->i].pos, c)),
                                     len(sub(mol.atoms[record->j].pos, c)));
        // Для линейного фрагмента дуги не существует: вместо неё отрезок.
        arc.linear = record->actual > 172;
        arc.emphasised = focus == record->center;
        arc.color = record->inRing ? theme::ANGLE_ARC : theme::ANGLE_ARC;
        if (record->inRing) arc.color = theme::ANGLE_ARC_RING;
        arc.alpha = fade * (arc.emphasised ? 1.0 : 0.72);
        if (!arc.linear) {
            const int steps = 36;
            arc.points.reserve(steps + 1);
            for (int k = 0; k <= steps; k++) {
                arc.points.push_back(
                    add(c, mul(slerp(arc.from, arc.to, static_cast<double>(k) / steps), arc.radius)));
            }
        }
        arcs.push_back(arc);
    }
    return arcs;
}

void MoleculeRenderer::drawAngles(Gdiplus::Graphics& g, Fonts& fonts, double progress) {
    (void)progress;
    labelRects.clear();
    const std::vector<AngleArc> arcs = angleArcs();
    if (arcs.empty()) return;

    const render::Font* font = fonts.mono(12.5, true);

    for (const AngleArc& arc : arcs) {
        const chem::AngleRecord* record = arc.record;
        const Vec3 c = arc.center;
        const Vec3 u = arc.from;
        const Vec3 w = arc.to;
        const double radius = arc.radius;
        const double alpha = arc.alpha;
        const Rgb color = arc.color;
        if (arc.linear) {
            drawLinearAngle(g, fonts, *record, c, u, w, radius, alpha);
            continue;
        }

        // подпись: отодвигаем её вдоль биссектрисы, пока она не перестанет
        // накладываться на уже нарисованные
        const Vec3 outward = norm(slerp(u, w, 0.5));
        const std::wstring text = degrees(record->actual);
        const double boxW = measure(g, text, *font).Width + 10;
        const double boxH = 18;

        bool placed = false;
        Projected spot{};
        for (int attempt = 0; attempt < 5; attempt++) {
            const Vec3 anchor = add(c, mul(outward, radius * (1.34 + attempt * 0.55)));
            const Projected point = project(anchor);
            if (!point.visible) break;
            // Подпись, заехавшую за шар, не показываем: она стоит в стороне от
            // атомов, поэтому хватает проверки одной точки.
            if (!occluder.visible(eyePos, anchor)) continue;
            const RectD rect{point.x - boxW / 2, point.y - boxH / 2, boxW, boxH};
            if (!labelCollides(rect)) { spot = point; labelRects.push_back(rect); placed = true; break; }
        }
        if (!placed) continue;

        fillRoundRect(g, spot.x - boxW / 2, spot.y - boxH / 2, boxW, boxH, 5,
                      Gdiplus::Color(static_cast<BYTE>(0.82 * alpha * 255), 8, 12, 24));
        strokeRoundRect(g, spot.x - boxW / 2, spot.y - boxH / 2, boxW, boxH, 5, toColor(color, alpha));
        drawTextCentered(g, text, *font, toColor(color, alpha), spot.x, spot.y, Rgb{8, 12, 24});
    }
}

void MoleculeRenderer::drawLinearAngle(Gdiplus::Graphics& g, Fonts& fonts,
                                       const chem::AngleRecord& record, const Vec3& center,
                                       const Vec3& u, const Vec3& w, double radius, double fade) {
    const Projected a = project(add(center, mul(u, radius)));
    const Projected b = project(add(center, mul(w, radius)));
    if (!a.visible || !b.visible) return;

    const Rgb color = record.inRing ? theme::ANGLE_ARC_RING : theme::ANGLE_ARC;
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double l = std::max(1e-6, std::hypot(dx, dy));
    const double px = -dy / l;
    const double py = dx / l;

    Gdiplus::Pen pen(toColor(color, fade), 1.6f);
    const Gdiplus::REAL dashes[] = {3.75f, 3.1f};
    pen.SetDashPattern(dashes, 2);
    g.DrawLine(&pen, static_cast<Gdiplus::REAL>(a.x), static_cast<Gdiplus::REAL>(a.y),
               static_cast<Gdiplus::REAL>(b.x), static_cast<Gdiplus::REAL>(b.y));

    const render::Font* font = fonts.mono(12.5, true);
    const std::wstring text = degrees(record.actual);
    const double boxW = measure(g, text, *font).Width + 10;

    double lx = 0, ly = 0;
    bool placed = false;
    for (int attempt = 0; attempt < 5; attempt++) {
        const double shift = 26 + attempt * 22;
        lx = (a.x + b.x) / 2 + px * shift;
        ly = (a.y + b.y) / 2 + py * shift;
        const RectD rect{lx - boxW / 2, ly - 9, boxW, 18};
        if (!labelCollides(rect)) { labelRects.push_back(rect); placed = true; break; }
    }
    if (!placed) return;

    fillRoundRect(g, lx - boxW / 2, ly - 9, boxW, 18, 5,
                  Gdiplus::Color(static_cast<BYTE>(0.82 * fade * 255), 8, 12, 24));
    strokeRoundRect(g, lx - boxW / 2, ly - 9, boxW, 18, 5, toColor(color, fade));
    drawTextCentered(g, text, *font, toColor(color, fade), lx, ly, Rgb{8, 12, 24});
}

// ---------------------------------------------------------------------------
// Дипольный момент
// ---------------------------------------------------------------------------

void MoleculeRenderer::drawDipole(Gdiplus::Graphics& g, Fonts& fonts, double progress) const {
    // Стрелка НАМЕРЕННО рисуется поверх всего и не проверяется на видимость.
    // Это не предмет сцены, а обозначение: она заведомо длиннее молекулы и
    // проходит сквозь неё насквозь. Ниже под неё кладётся тёмная подложка —
    // ровно затем, чтобы стрелка читалась и на светлых атомах. Отдать её
    // на откуп глубине значило бы спрятать половину.
    if (analysis->dipoleValue < 0.06) return;
    const double fade = clamp01((progress - 0.65) / 0.3);
    if (fade <= 0) return;

    const Vec3 dir = norm(analysis->dipoleVec);
    // Стрелка заведомо длиннее молекулы, иначе она теряется среди атомов.
    const double length = sceneRadius() * 2.1;
    const Projected a = project(mul(dir, -length * 0.5));
    const Projected b = project(mul(dir, length * 0.5));
    if (!a.visible || !b.visible) return;

    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double l = std::max(1e-6, std::hypot(dx, dy));
    const double ux = dx / l;
    const double uy = dy / l;
    const double headSize = std::max(10.0, 0.22 * b.scale);
    const auto f = [](double v) { return static_cast<Gdiplus::REAL>(v); };

    const double tipX = b.x - ux * headSize * 0.8;
    const double tipY = b.y - uy * headSize * 0.8;

    // тёмная подложка — стрелка остаётся читаемой поверх светлых атомов
    Gdiplus::Pen backing(Gdiplus::Color(static_cast<BYTE>(0.85 * fade * 255), 6, 9, 18), 7.0f);
    backing.SetStartCap(Gdiplus::LineCapRound);
    backing.SetEndCap(Gdiplus::LineCapRound);
    g.DrawLine(&backing, f(a.x), f(a.y), f(tipX), f(tipY));

    Gdiplus::GraphicsPath shaft;
    shaft.AddLine(f(a.x), f(a.y), f(tipX), f(tipY));
    glowPath(g, shaft, theme::DIPOLE, 3.0, fade);

    Gdiplus::Pen pen(toColor(theme::DIPOLE, fade), 3.0f);
    pen.SetStartCap(Gdiplus::LineCapRound);
    pen.SetEndCap(Gdiplus::LineCapRound);
    g.DrawLine(&pen, f(a.x), f(a.y), f(tipX), f(tipY));

    Gdiplus::PointF head[4] = {
        {f(b.x), f(b.y)},
        {f(b.x - ux * headSize - uy * headSize * 0.42), f(b.y - uy * headSize + ux * headSize * 0.42)},
        {f(b.x - ux * headSize * 0.7), f(b.y - uy * headSize * 0.7)},
        {f(b.x - ux * headSize + uy * headSize * 0.42), f(b.y - uy * headSize - ux * headSize * 0.42)},
    };
    Gdiplus::SolidBrush headBrush(toColor(theme::DIPOLE, fade));
    g.FillPolygon(&headBrush, head, 4);

    // поперечная чёрточка у хвоста — принятое в химии обозначение диполя
    Gdiplus::Pen tail(toColor(theme::DIPOLE, fade), 2.5f);
    g.DrawLine(&tail,
               f(a.x + uy * headSize * 0.42), f(a.y - ux * headSize * 0.42),
               f(a.x - uy * headSize * 0.42), f(a.y + ux * headSize * 0.42));

    // На стрелке подписываем ИЗМЕРЕННОЕ значение, если оно есть в справочнике:
    // расчёт даёт лишь направление и порядок величины.
    std::wstring caption = L"μ";
    if (analysis->molecule.meta != nullptr && analysis->molecule.meta->hasDipole) {
        wchar_t buffer[48];
        std::swprintf(buffer, 48, L"μ = %.2f Д", analysis->molecule.meta->dipole);
        caption = buffer;
    }
    const render::Font* font = fonts.ui(12, true);
    const double cx = (a.x + b.x) / 2 - uy * 20;
    const double cy = (a.y + b.y) / 2 - ux * 20;
    const double w = measure(g, caption, *font).Width + 12;
    fillRoundRect(g, cx - w / 2, cy - 9, w, 18, 5,
                  Gdiplus::Color(static_cast<BYTE>(0.85 * fade * 255), 8, 12, 24));
    drawTextCentered(g, caption, *font, toColor(theme::DIPOLE, fade), cx, cy, Rgb{8, 12, 24});
}

// ---------------------------------------------------------------------------
// Подписи атомов
// ---------------------------------------------------------------------------

void MoleculeRenderer::drawLabels(Gdiplus::Graphics& g, Fonts& fonts, double progress) const {
    const chem::Molecule& mol = analysis->molecule;
    const bool showHydrogens = mol.atoms.size() <= 14 || opts.style != Style::SpaceFill;

    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < mol.atoms.size(); i++) {
        if (projected[i].visible) order.push_back(i);
    }
    // от дальних к ближним: ближние подписи ложатся поверх
    std::stable_sort(order.begin(), order.end(), [this](std::size_t a, std::size_t b) {
        return projected[a].z > projected[b].z;
    });

    for (std::size_t index : order) {
        const chem::Atom& atom = mol.atoms[index];
        if (atom.el == "H" && !showHydrogens) continue;
        const double appear = atomAppearance(progress, index, mol.atoms.size());
        if (appear < 0.6) continue;

        const Projected& p = projected[index];
        const double r = atomRadius(atom.el) * p.scale;
        const double fontSize = std::max(9.0, std::min(22.0, r * 0.95));
        if (fontSize < 9.5) continue;

        // Ангстремов на пиксель на этой глубине — переводим размер плашки
        // из экранных единиц в мировые, чтобы проверить её видимость.
        const double perPixel = p.z / focal();
        const std::wstring symbol = toWide(atom.el);
        const double halfW = measure(g, symbol, *fonts.ui(fontSize, true)).Width * 0.5 * perPixel;
        const double halfH = fontSize * 0.5 * perPixel;
        if (!occluder.rectVisible(eyePos, atom.pos, halfW, halfH, camRight, camUp,
                                  static_cast<int>(index))) {
            continue;
        }

        const double alpha = (appear - 0.6) / 0.4;
        const Rgb elementColor = hexToRgb(chem::element(atom.el).color);
        const Gdiplus::Color textColor = luminance(elementColor) > 0.55
            ? Gdiplus::Color(static_cast<BYTE>(0.92 * alpha * 255), 10, 14, 26)
            : Gdiplus::Color(static_cast<BYTE>(0.95 * alpha * 255), 255, 255, 255);

        const Rgb sphere = fog(elementColor, p.z);
        drawTextCentered(g, symbol, *fonts.ui(fontSize, true), textColor, p.x, p.y, sphere);

        const chem::AtomAnalysis& info = analysis->atoms[index];
        if (std::fabs(info.charge) > 0.01) {
            const std::wstring sign = info.charge > 0 ? L"+" : L"−";
            const double magnitude = std::fabs(info.charge);
            std::wstring text;
            if (std::fabs(magnitude - 1) < 0.01) {
                text = sign;
            } else if (std::fabs(magnitude - std::floor(magnitude + 0.5)) < 1e-9) {
                text = std::to_wstring(static_cast<int>(magnitude + 0.5)) + sign;
            } else {
                text = L"δ" + sign;  // дробный заряд после усреднения — частичный
            }
            const Rgb chargeColor = info.charge > 0 ? theme::POSITIVE : theme::NEGATIVE;
            drawTextCentered(g, text, *fonts.ui(fontSize * 0.62, true), toColor(chargeColor, alpha),
                             p.x + r * 0.78, p.y - r * 0.72, theme::BACKGROUND);
        }
    }
}

}  // namespace render
