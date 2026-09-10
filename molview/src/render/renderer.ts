import type { Analysis } from '../chem/types';
import type { Vec3, Mat3 } from '../chem/vec';
import {
  v3, add, sub, mul, dot, cross, norm, len, applyMat, rotationAxis, matMul, IDENTITY,
} from '../chem/vec';
import { element } from '../chem/periodic';
import { THEME, rgb, rgbToCss, lighten, darken, mix, type Rgb } from './palette';

/**
 * СОБСТВЕННЫЙ ТРЁХМЕРНЫЙ РЕНДЕР НА CANVAS 2D
 *
 * Никаких сторонних библиотек: перспективная проекция считается вручную,
 * а объём создаётся сортировкой примитивов по глубине (алгоритм художника),
 * радиальными градиентами на шарах и продольными градиентами на связях.
 * Благодаря этому программа не зависит ни от WebGL, ни от интернета,
 * и легко переносится на любой язык, где есть двумерная графика.
 */

export type AngleMode = 'none' | 'selected' | 'all';
export type Style = 'ballstick' | 'spacefill' | 'wire';

export interface ViewOptions {
  showLonePairs: boolean;
  showOrbitals: boolean;
  showLabels: boolean;
  showDipole: boolean;
  showAngles: AngleMode;
  autoRotate: boolean;
  style: Style;
}

export const DEFAULT_OPTIONS: ViewOptions = {
  showLonePairs: true,
  showOrbitals: false,
  showLabels: true,
  showDipole: false,
  showAngles: 'selected',
  autoRotate: true,
  style: 'ballstick',
};

/** Элемент сцены, отсортированный по глубине. */
interface Primitive {
  depth: number;
  draw: () => void;
}

interface Projected {
  x: number;
  y: number;
  z: number;
  /** Пикселей на ангстрем на этой глубине. */
  scale: number;
  visible: boolean;
}

const ATOM_RADIUS_FACTOR = 0.34;
const BOND_RADIUS = 0.115;
const ENTRY_DURATION = 750;

export class MoleculeRenderer {
  private ctx: CanvasRenderingContext2D;
  private analysis: Analysis | null = null;
  private options: ViewOptions = { ...DEFAULT_OPTIONS };

  private rotation: Mat3 = IDENTITY.slice() as Mat3;
  private zoom = 1;
  private modelRadius = 4;
  private width = 0;
  private height = 0;
  private dpr = 1;

  private frame = 0;
  private entryStart = 0;
  private lastTime = 0;

  private projected: Projected[] = [];
  /** Прямоугольники уже нарисованных подписей — чтобы они не налезали друг на друга. */
  private labelRects: { x: number; y: number; w: number; h: number }[] = [];
  hoveredAtom: number | null = null;
  selectedAtom: number | null = null;
  onHoverChange: ((id: number | null) => void) | null = null;
  onSelectChange: ((id: number | null) => void) | null = null;

  private dragging = false;
  private lastPointer = { x: 0, y: 0 };
  private pointerMoved = false;

  constructor(private canvas: HTMLCanvasElement) {
    const ctx = canvas.getContext('2d');
    if (!ctx) throw new Error('Браузер не поддерживает Canvas 2D');
    this.ctx = ctx;
    this.attachEvents();
    this.resize();
  }

  // --- внешний интерфейс ----------------------------------------------------

  setAnalysis(analysis: Analysis, animate = true): void {
    this.analysis = analysis;
    this.selectedAtom = null;
    this.hoveredAtom = null;
    // проекции прошлой молекулы больше не действительны
    this.projected = [];

    let radius = 1;
    for (const atom of analysis.molecule.atoms) radius = Math.max(radius, len(atom.pos) + 0.9);
    this.modelRadius = radius;

    this.rotation = IDENTITY.slice() as Mat3;
    this.zoom = 1;
    this.entryStart = animate ? performance.now() : 0;
  }

  setOptions(options: Partial<ViewOptions>): void {
    this.options = { ...this.options, ...options };
  }

  getOptions(): ViewOptions {
    return { ...this.options };
  }

  resetView(): void {
    this.rotation = IDENTITY.slice() as Mat3;
    this.zoom = 1;
  }

  start(): void {
    if (this.frame !== 0) return;
    this.lastTime = performance.now();
    const loop = (time: number): void => {
      // Следующий кадр заказываем ДО отрисовки: если внутри случится сбой,
      // анимация не должна остановиться навсегда.
      this.frame = requestAnimationFrame(loop);

      const dt = Math.min(60, time - this.lastTime);
      this.lastTime = time;
      if (this.options.autoRotate && !this.dragging) {
        this.rotation = matMul(rotationAxis(v3(0, 1, 0), dt * 0.00022), this.rotation);
      }
      try {
        this.render(time);
      } catch (err) {
        console.error('Сбой при отрисовке кадра:', err);
      }
    };
    this.frame = requestAnimationFrame(loop);
  }

  stop(): void {
    if (this.frame !== 0) cancelAnimationFrame(this.frame);
    this.frame = 0;
  }

  resize(): void {
    const rect = this.canvas.getBoundingClientRect();
    this.dpr = Math.min(2.5, window.devicePixelRatio || 1);
    this.width = Math.max(1, rect.width);
    this.height = Math.max(1, rect.height);
    this.canvas.width = Math.round(this.width * this.dpr);
    this.canvas.height = Math.round(this.height * this.dpr);
  }

  /** Сохранение текущего кадра в PNG. */
  toDataUrl(): string {
    return this.canvas.toDataURL('image/png');
  }

  // --- управление мышью -----------------------------------------------------

  private attachEvents(): void {
    const c = this.canvas;

    c.addEventListener('pointerdown', (e) => {
      this.dragging = true;
      this.pointerMoved = false;
      this.lastPointer = { x: e.clientX, y: e.clientY };
      c.setPointerCapture(e.pointerId);
    });

    c.addEventListener('pointermove', (e) => {
      const rect = c.getBoundingClientRect();
      if (this.dragging) {
        const dx = e.clientX - this.lastPointer.x;
        const dy = e.clientY - this.lastPointer.y;
        if (Math.abs(dx) + Math.abs(dy) > 2) this.pointerMoved = true;
        this.lastPointer = { x: e.clientX, y: e.clientY };
        const speed = 0.0085;
        const spin = matMul(
          rotationAxis(v3(0, 1, 0), dx * speed),
          rotationAxis(v3(1, 0, 0), dy * speed),
        );
        this.rotation = matMul(spin, this.rotation);
      } else {
        this.updateHover(e.clientX - rect.left, e.clientY - rect.top);
      }
    });

    const release = (e: PointerEvent): void => {
      if (this.dragging && !this.pointerMoved) {
        const rect = c.getBoundingClientRect();
        this.updateHover(e.clientX - rect.left, e.clientY - rect.top);
        this.setSelected(this.hoveredAtom);
      }
      this.dragging = false;
    };
    c.addEventListener('pointerup', release);
    c.addEventListener('pointercancel', () => { this.dragging = false; });

    c.addEventListener('pointerleave', () => {
      if (this.hoveredAtom !== null) {
        this.hoveredAtom = null;
        this.onHoverChange?.(null);
      }
    });

    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      this.zoom = Math.max(0.35, Math.min(4, this.zoom * Math.exp(-e.deltaY * 0.0012)));
    }, { passive: false });

    c.addEventListener('dblclick', () => this.resetView());
  }

  setSelected(id: number | null): void {
    if (this.selectedAtom === id) return;
    this.selectedAtom = id;
    this.onSelectChange?.(id);
  }

  setHovered(id: number | null): void {
    if (this.hoveredAtom === id) return;
    this.hoveredAtom = id;
    this.onHoverChange?.(id);
  }

  private updateHover(px: number, py: number): void {
    if (!this.analysis) return;
    let found: number | null = null;
    let bestZ = Infinity;

    for (let i = 0; i < this.projected.length; i++) {
      const p = this.projected[i];
      if (!p.visible) continue;
      const atom = this.analysis.molecule.atoms[i];
      if (!atom) continue;
      const r = Math.max(9, this.atomRadius(atom.el) * p.scale);
      if ((px - p.x) ** 2 + (py - p.y) ** 2 <= r * r && p.z < bestZ) {
        bestZ = p.z;
        found = i;
      }
    }
    this.setHovered(found);
    this.canvas.style.cursor = found === null ? 'grab' : 'pointer';
  }

  // --- проекция -------------------------------------------------------------

  private atomRadius(el: string): number {
    const info = element(el);
    if (this.options.style === 'spacefill') return info.vdw * 0.62;
    if (this.options.style === 'wire') return Math.max(0.12, info.r * 0.18);
    return Math.max(0.26, info.r * ATOM_RADIUS_FACTOR + 0.12);
  }

  /**
   * Радиус сцены с учётом размеров самих шаров: в объёмной модели атомы
   * заметно крупнее, и камеру нужно отодвинуть, иначе молекула не поместится.
   */
  private sceneRadius(): number {
    if (!this.analysis) return this.modelRadius;
    let r = 0.8;
    for (const atom of this.analysis.molecule.atoms) {
      r = Math.max(r, len(atom.pos) + this.atomRadius(atom.el));
    }
    if (this.options.showOrbitals) r += 0.5;
    else if (this.options.showLonePairs) r += 0.35;
    return r + 0.3;
  }

  private cameraDistance(): number {
    return (this.sceneRadius() * 2.85) / this.zoom;
  }

  private focal(): number {
    return Math.min(this.width, this.height) * 1.25;
  }

  private project(p: Vec3): Projected {
    const camera = applyMat(this.rotation, p);
    const z = camera.z + this.cameraDistance();
    if (z < 0.15) return { x: 0, y: 0, z, scale: 0, visible: false };
    const s = this.focal() / z;
    return {
      x: this.width / 2 + camera.x * s,
      y: this.height / 2 - camera.y * s,
      z,
      scale: s,
      visible: true,
    };
  }

  /** Затуманивание дальних объектов — главный источник ощущения глубины. */
  private fog(color: Rgb, z: number): Rgb {
    const radius = this.sceneRadius();
    const near = this.cameraDistance() - radius;
    const far = this.cameraDistance() + radius;
    const t = Math.max(0, Math.min(1, (z - near) / Math.max(0.001, far - near)));
    return mix(color, rgb(THEME.fog), t * 0.55);
  }

  // --- отрисовка ------------------------------------------------------------

  private render(time: number): void {
    const ctx = this.ctx;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    this.drawBackground();

    if (!this.analysis) return;
    const mol = this.analysis.molecule;

    const progress = this.entryStart === 0
      ? 1
      : Math.max(0, Math.min(1, (time - this.entryStart) / ENTRY_DURATION));

    this.projected = mol.atoms.map((a) => this.project(a.pos));

    const primitives: Primitive[] = [];
    this.collectBonds(primitives, progress);
    this.collectAtoms(primitives, progress, time);
    if (this.options.showLonePairs) this.collectLonePairs(primitives, progress);
    if (this.options.showOrbitals) this.collectOrbitals(primitives, progress);

    primitives.sort((a, b) => b.depth - a.depth);
    for (const p of primitives) p.draw();

    if (this.options.showDipole) this.drawDipole(progress);
    if (this.options.showAngles !== 'none') this.drawAngles(progress);
    if (this.options.showLabels) this.drawLabels(progress);
  }

  private drawBackground(): void {
    const ctx = this.ctx;
    ctx.clearRect(0, 0, this.width, this.height);

    const cx = this.width / 2;
    const cy = this.height * 0.46;
    const radius = Math.max(this.width, this.height) * 0.78;
    const g = ctx.createRadialGradient(cx, cy, 0, cx, cy, radius);
    g.addColorStop(0, THEME.backgroundGlow);
    g.addColorStop(0.55, '#0b1020');
    g.addColorStop(1, THEME.background);
    ctx.fillStyle = g;
    ctx.fillRect(0, 0, this.width, this.height);

    // едва заметная сетка «лабораторного стола»
    ctx.save();
    ctx.strokeStyle = 'rgba(120, 170, 255, 0.045)';
    ctx.lineWidth = 1;
    const step = 46;
    ctx.beginPath();
    for (let x = (this.width / 2) % step; x < this.width; x += step) {
      ctx.moveTo(x, 0); ctx.lineTo(x, this.height);
    }
    for (let y = (this.height / 2) % step; y < this.height; y += step) {
      ctx.moveTo(0, y); ctx.lineTo(this.width, y);
    }
    ctx.stroke();
    ctx.restore();
  }

  // --- атомы ----------------------------------------------------------------

  private collectAtoms(out: Primitive[], progress: number, time: number): void {
    const mol = this.analysis!.molecule;
    const ctx = this.ctx;

    mol.atoms.forEach((atom, index) => {
      const p = this.projected[index];
      if (!p.visible) return;

      const appear = atomAppearance(progress, index, mol.atoms.length);
      if (appear <= 0.001) return;

      const baseRadius = this.atomRadius(atom.el) * p.scale * appear;
      const color = this.fog(rgb(element(atom.el).color), p.z);
      const isHovered = this.hoveredAtom === index;
      const isSelected = this.selectedAtom === index;

      out.push({
        depth: p.z,
        draw: () => {
          const r = baseRadius * (isHovered ? 1.1 : 1);

          if (isSelected || isHovered) {
            const pulse = isSelected ? 0.5 + 0.5 * Math.sin(time * 0.005) : 0.55;
            const glow = ctx.createRadialGradient(p.x, p.y, r * 0.7, p.x, p.y, r * 2.3);
            glow.addColorStop(0, `rgba(94, 234, 212, ${0.36 * pulse})`);
            glow.addColorStop(1, 'rgba(94, 234, 212, 0)');
            ctx.fillStyle = glow;
            ctx.beginPath();
            ctx.arc(p.x, p.y, r * 2.3, 0, Math.PI * 2);
            ctx.fill();
          }

          const g = ctx.createRadialGradient(
            p.x - r * 0.36, p.y - r * 0.42, r * 0.05,
            p.x, p.y, r * 1.02,
          );
          g.addColorStop(0, rgbToCss(lighten(color, 0.62)));
          g.addColorStop(0.42, rgbToCss(lighten(color, 0.1)));
          g.addColorStop(1, rgbToCss(darken(color, 0.52)));

          ctx.beginPath();
          ctx.arc(p.x, p.y, r, 0, Math.PI * 2);
          ctx.fillStyle = g;
          ctx.fill();

          // ободок — отделяет атом от соседей
          ctx.lineWidth = Math.max(0.6, r * 0.07);
          ctx.strokeStyle = isSelected
            ? 'rgba(94, 234, 212, 0.95)'
            : `rgba(255, 255, 255, ${0.16 + (isHovered ? 0.3 : 0)})`;
          ctx.stroke();

          // блик
          const hl = ctx.createRadialGradient(
            p.x - r * 0.34, p.y - r * 0.42, 0,
            p.x - r * 0.34, p.y - r * 0.42, r * 0.42,
          );
          hl.addColorStop(0, 'rgba(255,255,255,0.55)');
          hl.addColorStop(1, 'rgba(255,255,255,0)');
          ctx.fillStyle = hl;
          ctx.beginPath();
          ctx.arc(p.x - r * 0.34, p.y - r * 0.42, r * 0.42, 0, Math.PI * 2);
          ctx.fill();
        },
      });
    });
  }

  // --- связи ----------------------------------------------------------------

  private collectBonds(out: Primitive[], progress: number): void {
    const mol = this.analysis!.molecule;
    if (this.options.style === 'spacefill') return;

    for (const bond of mol.bonds) {
      const pa = this.projected[bond.a];
      const pb = this.projected[bond.b];
      if (!pa.visible || !pb.visible) continue;

      const grow = Math.max(0, Math.min(1, (progress - 0.12) / 0.55));
      if (grow <= 0.001) continue;

      const world = mul(add(mol.atoms[bond.a].pos, mol.atoms[bond.b].pos), 0.5);
      const mid = this.project(world);

      const colorA = this.fog(rgb(element(mol.atoms[bond.a].el).color), (pa.z + mid.z) / 2);
      const colorB = this.fog(rgb(element(mol.atoms[bond.b].el).color), (pb.z + mid.z) / 2);

      // перпендикуляр к оси связи в экранных координатах
      const dx = pb.x - pa.x;
      const dy = pb.y - pa.y;
      const l = Math.hypot(dx, dy) || 1;
      const perp = { x: -dy / l, y: dx / l };

      const count = bond.order >= 2.9 ? 3 : bond.order >= 1.4 ? 2 : 1;
      const separation = BOND_RADIUS * 2.05;
      const dashed = bond.aromatic || bond.delocalized;

      for (let k = 0; k < count; k++) {
        const offset = count === 1 ? 0 : (k - (count - 1) / 2) * separation;
        // вторая линия у ароматической связи рисуется пунктиром — так на схеме
        // видно, что порядок связи дробный
        const isDashedLine = dashed && k === count - 1;
        const widthScale = isDashedLine ? 0.55 : 1;

        this.pushHalfBond(out, pa, mid, perp, offset, colorA, grow, widthScale, isDashedLine);
        this.pushHalfBond(out, pb, mid, perp, offset, colorB, grow, widthScale, isDashedLine);
      }
    }
  }

  private pushHalfBond(
    out: Primitive[],
    from: Projected,
    to: Projected,
    perp: { x: number; y: number },
    offsetAngstrom: number,
    color: Rgb,
    grow: number,
    widthScale: number,
    dashed: boolean,
  ): void {
    const ctx = this.ctx;
    const w1 = BOND_RADIUS * from.scale * widthScale;
    const w2 = BOND_RADIUS * to.scale * widthScale;
    const o1 = offsetAngstrom * from.scale;
    const o2 = offsetAngstrom * to.scale;

    const ax = from.x + perp.x * o1;
    const ay = from.y + perp.y * o1;
    const bx = ax + (to.x + perp.x * o2 - ax) * grow;
    const by = ay + (to.y + perp.y * o2 - ay) * grow;

    out.push({
      depth: (from.z + to.z) / 2,
      draw: () => {
        if (dashed) {
          ctx.save();
          ctx.setLineDash([w1 * 2.2, w1 * 1.8]);
          ctx.lineWidth = w1 * 1.7;
          ctx.strokeStyle = rgbToCss(lighten(color, 0.18), 0.8);
          ctx.lineCap = 'round';
          ctx.beginPath();
          ctx.moveTo(ax, ay);
          ctx.lineTo(bx, by);
          ctx.stroke();
          ctx.restore();
          return;
        }

        ctx.beginPath();
        ctx.moveTo(ax + perp.x * w1, ay + perp.y * w1);
        ctx.lineTo(bx + perp.x * w2, by + perp.y * w2);
        ctx.lineTo(bx - perp.x * w2, by - perp.y * w2);
        ctx.lineTo(ax - perp.x * w1, ay - perp.y * w1);
        ctx.closePath();

        // поперечный градиент имитирует цилиндр
        const gx = (ax + bx) / 2;
        const gy = (ay + by) / 2;
        const gw = (w1 + w2) / 2;
        const g = ctx.createLinearGradient(
          gx + perp.x * gw, gy + perp.y * gw,
          gx - perp.x * gw, gy - perp.y * gw,
        );
        g.addColorStop(0, rgbToCss(darken(color, 0.55)));
        g.addColorStop(0.38, rgbToCss(lighten(color, 0.3)));
        g.addColorStop(1, rgbToCss(darken(color, 0.55)));
        ctx.fillStyle = g;
        ctx.fill();
      },
    });
  }

  // --- неподелённые пары ----------------------------------------------------

  private collectLonePairs(out: Primitive[], progress: number): void {
    if (progress < 0.5) return;
    const analysis = this.analysis!;
    const ctx = this.ctx;
    const fade = Math.min(1, (progress - 0.5) / 0.35);

    // Пары концевых атомов (три пары на каждом атоме хлора и т.п.) сильно
    // засоряют картинку, а форму молекулы определяют пары ЦЕНТРАЛЬНЫХ атомов.
    // Поэтому по умолчанию показываем только их — плюс пары выбранного атома.
    const hasCentral = analysis.atoms.some((a) => a.isCentral && a.lonePairDirs.length > 0);

    for (const info of analysis.atoms) {
      if (info.lonePairDirs.length === 0) continue;
      if (hasCentral && !info.isCentral && this.selectedAtom !== info.id) continue;
      const atom = analysis.molecule.atoms[info.id];
      const reach = this.atomRadius(atom.el) + 0.42;

      for (const dir of info.lonePairDirs) {
        const centerWorld = add(atom.pos, mul(dir, reach));
        const p = this.project(centerWorld);
        if (!p.visible) continue;

        const base = this.project(atom.pos);
        const axisX = p.x - base.x;
        const axisY = p.y - base.y;
        const axisLen = Math.hypot(axisX, axisY);
        const angle = axisLen > 1 ? Math.atan2(axisY, axisX) : 0;

        // сплющиваем облако, когда пара смотрит «на зрителя»
        const foreshorten = Math.max(0.32, Math.min(1, axisLen / (reach * p.scale)));
        const rx = 0.30 * p.scale;
        const ry = rx * (0.55 + 0.45 * foreshorten);
        const dimmed = this.hoveredAtom !== null && this.hoveredAtom !== info.id;

        out.push({
          depth: p.z,
          draw: () => {
            const alpha = (dimmed ? 0.28 : 0.62) * fade;
            ctx.save();
            ctx.translate(p.x, p.y);
            ctx.rotate(angle);

            const g = ctx.createRadialGradient(0, 0, 0, 0, 0, rx);
            g.addColorStop(0, `rgba(147, 214, 255, ${alpha})`);
            g.addColorStop(0.6, `rgba(125, 211, 252, ${alpha * 0.5})`);
            g.addColorStop(1, 'rgba(125, 211, 252, 0)');
            ctx.fillStyle = g;
            ctx.beginPath();
            ctx.ellipse(0, 0, rx, ry, 0, 0, Math.PI * 2);
            ctx.fill();

            // два электрона
            const dotR = Math.max(1.1, 0.045 * p.scale);
            const gap = ry * 0.5;
            ctx.fillStyle = `rgba(224, 245, 255, ${Math.min(1, alpha * 1.5)})`;
            for (const sign of [-1, 1]) {
              ctx.beginPath();
              ctx.arc(0, sign * gap, dotR, 0, Math.PI * 2);
              ctx.fill();
            }
            ctx.restore();
          },
        });
      }
    }
  }

  // --- гибридные орбитали ---------------------------------------------------

  private collectOrbitals(out: Primitive[], progress: number): void {
    if (progress < 0.55) return;
    const analysis = this.analysis!;
    const ctx = this.ctx;
    const fade = Math.min(1, (progress - 0.55) / 0.35);

    for (const info of analysis.atoms) {
      if (info.steric < 2 || info.el === 'H') continue;
      if (this.selectedAtom !== null && this.selectedAtom !== info.id) continue;
      if (this.selectedAtom === null && !info.isCentral) continue;

      const atom = analysis.molecule.atoms[info.id];
      const lobeLength = 1.15;

      info.orbitalDirs.forEach((dir, k) => {
        const isLonePair = k >= info.sigma;
        const tipWorld = add(atom.pos, mul(dir, lobeLength));
        const tip = this.project(tipWorld);
        const base = this.project(atom.pos);
        if (!tip.visible || !base.visible) return;

        const dx = tip.x - base.x;
        const dy = tip.y - base.y;
        const l = Math.hypot(dx, dy);
        const ux = l > 0.001 ? dx / l : 1;
        const uy = l > 0.001 ? dy / l : 0;
        const px = -uy;
        const py = ux;
        const width = 0.34 * base.scale;

        out.push({
          depth: (tip.z + base.z) / 2 + 0.01,
          draw: () => {
            const color = isLonePair ? '125, 211, 252' : '192, 132, 252';
            const alpha = 0.30 * fade;

            ctx.save();
            ctx.beginPath();
            ctx.moveTo(base.x, base.y);
            ctx.bezierCurveTo(
              base.x + ux * l * 0.28 + px * width, base.y + uy * l * 0.28 + py * width,
              base.x + ux * l * 0.82 + px * width * 0.85, base.y + uy * l * 0.82 + py * width * 0.85,
              tip.x, tip.y,
            );
            ctx.bezierCurveTo(
              base.x + ux * l * 0.82 - px * width * 0.85, base.y + uy * l * 0.82 - py * width * 0.85,
              base.x + ux * l * 0.28 - px * width, base.y + uy * l * 0.28 - py * width,
              base.x, base.y,
            );
            ctx.closePath();

            const g = ctx.createLinearGradient(base.x, base.y, tip.x, tip.y);
            g.addColorStop(0, `rgba(${color}, ${alpha * 0.35})`);
            g.addColorStop(0.55, `rgba(${color}, ${alpha})`);
            g.addColorStop(1, `rgba(${color}, ${alpha * 0.15})`);
            ctx.fillStyle = g;
            ctx.fill();
            ctx.strokeStyle = `rgba(${color}, ${alpha * 1.3})`;
            ctx.lineWidth = 1;
            ctx.stroke();

            // малый «задний» лепесток гибридной орбитали
            const backLen = l * 0.26;
            ctx.beginPath();
            ctx.moveTo(base.x, base.y);
            ctx.bezierCurveTo(
              base.x - ux * backLen * 0.5 + px * width * 0.5, base.y - uy * backLen * 0.5 + py * width * 0.5,
              base.x - ux * backLen + px * width * 0.2, base.y - uy * backLen + py * width * 0.2,
              base.x - ux * backLen, base.y - uy * backLen,
            );
            ctx.bezierCurveTo(
              base.x - ux * backLen - px * width * 0.2, base.y - uy * backLen - py * width * 0.2,
              base.x - ux * backLen * 0.5 - px * width * 0.5, base.y - uy * backLen * 0.5 - py * width * 0.5,
              base.x, base.y,
            );
            ctx.closePath();
            ctx.fillStyle = `rgba(${color}, ${alpha * 0.35})`;
            ctx.fill();
            ctx.restore();
          },
        });
      });
    }
  }

  // --- валентные углы -------------------------------------------------------

  private drawAngles(progress: number): void {
    if (progress < 0.6) return;
    const analysis = this.analysis!;
    const ctx = this.ctx;
    const fade = Math.min(1, (progress - 0.6) / 0.3);
    const focus = this.selectedAtom ?? this.hoveredAtom;
    this.labelRects = [];

    let records = analysis.angles;

    // В режиме «все» на большой молекуле углы с участием водорода дают десятки
    // почти одинаковых дуг и превращают картинку в кашу. Оставляем скелет.
    if (this.options.showAngles === 'all' && analysis.molecule.atoms.length > 8) {
      records = records.filter(
        (r) => analysis.molecule.atoms[r.i].el !== 'H' && analysis.molecule.atoms[r.j].el !== 'H',
      );
    }

    if (this.options.showAngles === 'selected') {
      if (focus === null) {
        // ничего не выбрано — показываем углы самого «главного» центра
        const central = analysis.atoms
          .filter((a) => a.isCentral && a.el !== 'H')
          .sort((a, b) => b.sigma - a.sigma)[0];
        if (!central) return;
        records = records.filter((r) => r.center === central.id);
      } else {
        records = records.filter((r) => r.center === focus);
      }
    }
    if (records.length === 0) return;

    const mol = analysis.molecule;

    for (const record of records) {
      const c = mol.atoms[record.center].pos;
      const u = norm(sub(mol.atoms[record.i].pos, c));
      const w = norm(sub(mol.atoms[record.j].pos, c));
      const radius = 0.46 * Math.min(
        len(sub(mol.atoms[record.i].pos, c)),
        len(sub(mol.atoms[record.j].pos, c)),
      );

      // Для линейного фрагмента дуги не существует: рисуем прямую-диаметр.
      if (record.actual > 172) {
        this.drawLinearAngle(record, c, u, w, radius, fade);
        continue;
      }

      const steps = 36;
      const points: Projected[] = [];
      for (let s = 0; s <= steps; s++) {
        points.push(this.project(add(c, mul(slerp(u, w, s / steps), radius))));
      }
      if (points.some((p) => !p.visible)) continue;

      const color = record.inRing ? THEME.angleArcRing : THEME.angleArc;
      const emphasised = focus === record.center;

      ctx.save();
      ctx.globalAlpha = fade * (emphasised ? 1 : 0.72);

      // заливка сектора
      ctx.beginPath();
      const centerPoint = this.project(c);
      ctx.moveTo(centerPoint.x, centerPoint.y);
      points.forEach((p) => ctx.lineTo(p.x, p.y));
      ctx.closePath();
      ctx.fillStyle = record.inRing ? 'rgba(251, 113, 133, 0.10)' : 'rgba(94, 234, 212, 0.10)';
      ctx.fill();

      // дуга
      ctx.beginPath();
      points.forEach((p, k) => (k === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
      ctx.strokeStyle = color;
      ctx.lineWidth = emphasised ? 2 : 1.4;
      ctx.shadowColor = color;
      ctx.shadowBlur = 8;
      ctx.stroke();
      ctx.shadowBlur = 0;

      // подпись: отодвигаем её вдоль биссектрисы, пока она не перестанет
      // накладываться на уже нарисованные
      const outward = norm(slerp(u, w, 0.5));
      const text = `${record.actual.toFixed(1)}°`;
      ctx.font = '600 12.5px ui-monospace, "SF Mono", Menlo, Consolas, monospace';
      const boxW = ctx.measureText(text).width + 10;
      const boxH = 18;

      let spot: { x: number; y: number } | null = null;
      for (let attempt = 0; attempt < 5; attempt++) {
        const point = this.project(add(c, mul(outward, radius * (1.34 + attempt * 0.55))));
        if (!point.visible) break;
        const rect = { x: point.x - boxW / 2, y: point.y - boxH / 2, w: boxW, h: boxH };
        if (!this.labelCollides(rect)) { spot = point; this.labelRects.push(rect); break; }
      }
      if (spot === null) { ctx.restore(); continue; }

      ctx.fillStyle = 'rgba(8, 12, 24, 0.82)';
      roundRect(ctx, spot.x - boxW / 2, spot.y - boxH / 2, boxW, boxH, 5);
      ctx.fill();
      ctx.strokeStyle = color;
      ctx.lineWidth = 1;
      ctx.stroke();

      ctx.fillStyle = color;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText(text, spot.x, spot.y + 0.5);
      ctx.restore();
    }
  }

  /** Пересекается ли прямоугольник подписи с уже нарисованными. */
  private labelCollides(r: { x: number; y: number; w: number; h: number }): boolean {
    const gap = 3;
    return this.labelRects.some((o) =>
      r.x < o.x + o.w + gap && r.x + r.w + gap > o.x &&
      r.y < o.y + o.h + gap && r.y + r.h + gap > o.y);
  }

  /** Указатель для угла, близкого к 180°: дуга выродилась бы в отрезок. */
  private drawLinearAngle(
    record: { actual: number; inRing: boolean },
    center: Vec3, u: Vec3, w: Vec3, radius: number, fade: number,
  ): void {
    const ctx = this.ctx;
    const a = this.project(add(center, mul(u, radius)));
    const b = this.project(add(center, mul(w, radius)));
    if (!a.visible || !b.visible) return;

    const color = record.inRing ? THEME.angleArcRing : THEME.angleArc;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const l = Math.hypot(dx, dy) || 1;
    const px = -dy / l;
    const py = dx / l;

    ctx.save();
    ctx.globalAlpha = fade;
    ctx.setLineDash([6, 5]);
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.6;
    ctx.shadowColor = color;
    ctx.shadowBlur = 8;
    ctx.beginPath();
    ctx.moveTo(a.x, a.y);
    ctx.lineTo(b.x, b.y);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.shadowBlur = 0;

    const text = `${record.actual.toFixed(1)}°`;
    ctx.font = '600 12.5px ui-monospace, "SF Mono", Menlo, Consolas, monospace';
    const boxW = ctx.measureText(text).width + 10;

    let lx = 0, ly = 0, placed = false;
    for (let attempt = 0; attempt < 5; attempt++) {
      const shift = 26 + attempt * 22;
      lx = (a.x + b.x) / 2 + px * shift;
      ly = (a.y + b.y) / 2 + py * shift;
      const rect = { x: lx - boxW / 2, y: ly - 9, w: boxW, h: 18 };
      if (!this.labelCollides(rect)) { this.labelRects.push(rect); placed = true; break; }
    }
    if (!placed) { ctx.restore(); return; }

    ctx.fillStyle = 'rgba(8, 12, 24, 0.82)';
    roundRect(ctx, lx - boxW / 2, ly - 9, boxW, 18, 5);
    ctx.fill();
    ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.stroke();
    ctx.fillStyle = color;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(text, lx, ly + 0.5);
    ctx.restore();
  }

  // --- дипольный момент -----------------------------------------------------

  private drawDipole(progress: number): void {
    const analysis = this.analysis!;
    if (analysis.dipoleValue < 0.06) return;
    const ctx = this.ctx;
    const fade = Math.max(0, Math.min(1, (progress - 0.65) / 0.3));
    if (fade <= 0) return;

    const dir = norm(analysis.dipoleVec);
    // Стрелка заведомо длиннее молекулы, иначе она теряется среди атомов.
    const length = this.sceneRadius() * 2.1;
    const tail = mul(dir, -length * 0.5);
    const head = mul(dir, length * 0.5);

    const a = this.project(tail);
    const b = this.project(head);
    if (!a.visible || !b.visible) return;

    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const l = Math.hypot(dx, dy) || 1;
    const ux = dx / l;
    const uy = dy / l;
    const headSize = Math.max(10, 0.22 * b.scale);

    ctx.save();
    ctx.globalAlpha = fade;
    ctx.lineCap = 'round';

    // тёмная подложка — стрелка остаётся читаемой поверх светлых атомов
    ctx.strokeStyle = 'rgba(6, 9, 18, 0.85)';
    ctx.lineWidth = 7;
    ctx.beginPath();
    ctx.moveTo(a.x, a.y);
    ctx.lineTo(b.x - ux * headSize * 0.8, b.y - uy * headSize * 0.8);
    ctx.stroke();

    ctx.strokeStyle = THEME.dipole;
    ctx.fillStyle = THEME.dipole;
    ctx.lineWidth = 3;
    ctx.shadowColor = THEME.dipole;
    ctx.shadowBlur = 12;

    ctx.beginPath();
    ctx.moveTo(a.x, a.y);
    ctx.lineTo(b.x - ux * headSize * 0.8, b.y - uy * headSize * 0.8);
    ctx.stroke();

    ctx.beginPath();
    ctx.moveTo(b.x, b.y);
    ctx.lineTo(b.x - ux * headSize + -uy * headSize * 0.42, b.y - uy * headSize + ux * headSize * 0.42);
    ctx.lineTo(b.x - ux * headSize * 0.7, b.y - uy * headSize * 0.7);
    ctx.lineTo(b.x - ux * headSize - -uy * headSize * 0.42, b.y - uy * headSize - ux * headSize * 0.42);
    ctx.closePath();
    ctx.fill();

    // поперечная чёрточка у хвоста — принятое в химии обозначение диполя
    ctx.shadowBlur = 0;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(a.x - -uy * headSize * 0.42, a.y - ux * headSize * 0.42);
    ctx.lineTo(a.x + -uy * headSize * 0.42, a.y + ux * headSize * 0.42);
    ctx.stroke();

    // На стрелке подписываем ИЗМЕРЕННОЕ значение, если оно есть в справочнике:
    // расчёт даёт лишь направление и порядок величины.
    const measured = analysis.molecule.meta?.dipole;
    const caption = measured !== undefined ? `μ = ${measured.toFixed(2)} Д` : 'μ';
    ctx.font = '600 12px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    const cx = (a.x + b.x) / 2 - uy * 20;
    const cy = (a.y + b.y) / 2 - ux * 20;
    const w = ctx.measureText(caption).width + 12;
    ctx.fillStyle = 'rgba(8, 12, 24, 0.85)';
    roundRect(ctx, cx - w / 2, cy - 9, w, 18, 5);
    ctx.fill();
    ctx.fillStyle = THEME.dipole;
    ctx.fillText(caption, cx, cy + 0.5);
    ctx.restore();
  }

  // --- подписи атомов -------------------------------------------------------

  private drawLabels(progress: number): void {
    const analysis = this.analysis!;
    const ctx = this.ctx;
    const mol = analysis.molecule;
    const showHydrogens = mol.atoms.length <= 14 || this.options.style !== 'spacefill';

    const items = mol.atoms
      .map((atom, index) => ({ atom, index, p: this.projected[index] }))
      .filter(({ p }) => p.visible)
      .sort((a, b) => b.p.z - a.p.z);

    ctx.save();
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    for (const { atom, index, p } of items) {
      if (atom.el === 'H' && !showHydrogens) continue;
      const appear = atomAppearance(progress, index, mol.atoms.length);
      if (appear < 0.6) continue;

      const r = this.atomRadius(atom.el) * p.scale;
      const fontSize = Math.max(9, Math.min(22, r * 0.95));
      if (fontSize < 9.5) continue;

      ctx.font = `700 ${fontSize.toFixed(1)}px "Inter", system-ui, -apple-system, sans-serif`;
      const luminance = luminanceOf(element(atom.el).color);
      ctx.fillStyle = luminance > 0.55 ? 'rgba(10, 14, 26, 0.92)' : 'rgba(255, 255, 255, 0.95)';
      ctx.globalAlpha = (appear - 0.6) / 0.4;
      ctx.fillText(atom.el, p.x, p.y);

      const analysisAtom = analysis.atoms[index];
      if (Math.abs(analysisAtom.charge) > 0.01) {
        const sign = analysisAtom.charge > 0 ? '+' : '−';
        const magnitude = Math.abs(analysisAtom.charge);
        const text = Math.abs(magnitude - 1) < 0.01 ? sign
          : Number.isInteger(magnitude) ? `${magnitude}${sign}` : `δ${sign}`;
        ctx.font = `700 ${(fontSize * 0.62).toFixed(1)}px "Inter", system-ui, sans-serif`;
        ctx.fillStyle = analysisAtom.charge > 0 ? '#fca5a5' : '#93c5fd';
        ctx.fillText(text, p.x + r * 0.78, p.y - r * 0.72);
      }
      ctx.globalAlpha = 1;
    }
    ctx.restore();
  }
}

// ---------------------------------------------------------------------------
// вспомогательные функции
// ---------------------------------------------------------------------------

/** Появление атома при загрузке молекулы: волной от центра, с лёгким «перелётом». */
function atomAppearance(progress: number, index: number, total: number): number {
  if (progress >= 1) return 1;
  const delay = (index / Math.max(1, total)) * 0.45;
  const t = Math.max(0, Math.min(1, (progress - delay) / 0.55));
  const eased = 1 - Math.pow(1 - t, 3);
  return eased * (1 + 0.12 * Math.sin(t * Math.PI));
}

/** Интерполяция направлений по дуге большого круга. */
function slerp(a: Vec3, b: Vec3, t: number): Vec3 {
  const cosine = Math.max(-1, Math.min(1, dot(a, b)));
  const omega = Math.acos(cosine);
  if (omega < 1e-4) return a;
  if (Math.PI - omega < 1e-4) {
    // векторы противоположны — идём через произвольную перпендикулярную ось
    const axis = norm(cross(a, Math.abs(a.x) < 0.9 ? v3(1, 0, 0) : v3(0, 1, 0)));
    return applyMat(rotationAxis(axis, omega * t), a);
  }
  const s = Math.sin(omega);
  return add(mul(a, Math.sin((1 - t) * omega) / s), mul(b, Math.sin(t * omega) / s));
}

function roundRect(ctx: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, r: number): void {
  ctx.beginPath();
  ctx.moveTo(x + r, y);
  ctx.arcTo(x + w, y, x + w, y + h, r);
  ctx.arcTo(x + w, y + h, x, y + h, r);
  ctx.arcTo(x, y + h, x, y, r);
  ctx.arcTo(x, y, x + w, y, r);
  ctx.closePath();
}

function luminanceOf(hex: string): number {
  const c = rgb(hex);
  return (0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b) / 255;
}
