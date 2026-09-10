/* Генератор презентации для защиты. Запуск: node docs/make-presentation.cjs */
const path = require('path');
const pptxgen = require('pptxgenjs');

const pres = new pptxgen();
pres.layout = 'LAYOUT_WIDE';           // 13.33 × 7.5 дюйма
pres.author = 'Индивидуальный проект';
pres.title = 'Компьютерное моделирование строения молекул';

const W = 13.33, H = 7.5;
const M = 0.7;                          // поле

// --- палитра: та же «лабораторная» тема, что и в самой программе ------------
const BG = '0B1220';
const CARD = '16203A';
const TEAL = '5EEAD4';
const VIOLET = 'A78BFA';
const AMBER = 'FBBF24';
const ROSE = 'FB7185';
const TEXT = 'E8EEFF';
const DIM = '9AA6C4';
const FAINT = '6B7796';

const EL = {
  H: 'F2F5FA', C: '6E6E7C', N: '4A6CFF', O: 'FF3B30', S: 'FFD429',
  F: '8FE04A', Cl: '33DD55', Xe: '42B8D0', B: 'FFB5B5', P: 'FF8C1A',
};

const HEAD = 'Cambria';
const BODY = 'Calibri';
const MONO = 'Courier New';

// ---------------------------------------------------------------------------
// помощники
// ---------------------------------------------------------------------------

function slide(dark = true) {
  const s = pres.addSlide();
  s.background = { color: dark ? BG : '101A30' };
  return s;
}

function title(s, text, sub) {
  s.addText(text, {
    x: M, y: 0.45, w: W - 2 * M, h: 0.8,
    fontFace: HEAD, fontSize: 34, bold: true, color: TEXT,
    align: 'left', valign: 'middle', isTextBox: true, margin: 0,
  });
  if (sub) {
    s.addText(sub, {
      x: M, y: 1.24, w: W - 2 * M, h: 0.42,
      fontFace: BODY, fontSize: 15, color: DIM,
      align: 'left', valign: 'top', isTextBox: true, margin: 0,
    });
  }
}

/** Карточка-подложка. */
function card(s, x, y, w, h, fill = CARD) {
  s.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.12,
    fill: { color: fill },
    line: { color: '2A3A5E', width: 1 },
  });
}

/** Кружок с номером — повторяющийся мотив всей презентации. */
function badge(s, x, y, text, color = TEAL, d = 0.46) {
  s.addText(String(text), {
    shape: pres.ShapeType.ellipse,
    x, y, w: d, h: d,
    fill: { color },
    fontFace: BODY, fontSize: 15, bold: true, color: '0B1220',
    align: 'center', valign: 'middle', isTextBox: true, margin: 0,
  });
}

function body(s, text, x, y, w, h, opts = {}) {
  s.addText(text, {
    x, y, w, h,
    fontFace: opts.face || BODY,
    fontSize: opts.size ?? 15,
    color: opts.color || DIM,
    bold: !!opts.bold, italic: !!opts.italic,
    align: opts.align || 'left', valign: opts.valign || 'top',
    lineSpacingMultiple: opts.lsm ?? 1.15,
    isTextBox: true, margin: 0,
  });
}

function bullets(s, items, x, y, w, h, opts = {}) {
  s.addText(items.map((it, k) => ({
    text: typeof it === 'string' ? it : it.text,
    options: {
      bullet: true, breakLine: k < items.length - 1,
      color: (typeof it === 'object' && it.color) || opts.color || DIM,
      bold: typeof it === 'object' && !!it.bold,
    },
  })), {
    x, y, w, h,
    fontFace: BODY, fontSize: opts.size ?? 15, color: opts.color || DIM,
    paraSpaceAfter: opts.gap ?? 9,
    lineSpacingMultiple: 1.1,
    isTextBox: true, margin: 0,
  });
}

/** Крупное число с подписью. */
function stat(s, x, y, w, value, label, color = TEAL) {
  s.addText(value, {
    x, y, w, h: 0.85,
    fontFace: HEAD, fontSize: 46, bold: true, color,
    align: 'center', valign: 'middle', isTextBox: true, margin: 0,
  });
  s.addText(label, {
    x, y: y + 0.85, w, h: 0.5,
    fontFace: BODY, fontSize: 12.5, color: DIM,
    align: 'center', valign: 'top', isTextBox: true, margin: 0,
  });
}

/**
 * Молекула из настоящих фигур PowerPoint: линии-связи и кружки-атомы.
 * atoms: [{ x, y, el, r? }] в условных единицах, начало — в центре рисунка.
 */
function molecule(s, cx, cy, scale, atoms, bonds, opts = {}) {
  const px = (a) => cx + a.x * scale;
  const py = (a) => cy + a.y * scale;
  const rad = (a) => (a.r ?? (a.el === 'H' ? 0.16 : 0.23)) * (opts.rScale ?? 1);

  // связи
  for (const [i, j, order = 1] of bonds) {
    const a = atoms[i], b = atoms[j];
    const x1 = px(a), y1 = py(a), x2 = px(b), y2 = py(b);
    const dx = x2 - x1, dy = y2 - y1;
    const len = Math.hypot(dx, dy) || 1;
    const nx = -dy / len, ny = dx / len;
    const count = order >= 2 ? Math.min(3, Math.round(order)) : 1;
    const sep = 0.055;

    for (let k = 0; k < count; k++) {
      const off = count === 1 ? 0 : (k - (count - 1) / 2) * sep * 2;
      const ax = x1 + nx * off, ay = y1 + ny * off;
      const bx = x2 + nx * off, by = y2 + ny * off;
      s.addShape(pres.ShapeType.line, {
        x: Math.min(ax, bx), y: Math.min(ay, by),
        w: Math.abs(bx - ax), h: Math.abs(by - ay),
        flipH: bx < ax, flipV: by < ay,
        line: { color: opts.bondColor || '8A93AD', width: opts.bondWidth ?? 3.2 },
      });
    }
  }

  // атомы
  for (const a of atoms) {
    const r = rad(a);
    s.addText(opts.hideLabels ? '' : a.el, {
      shape: pres.ShapeType.ellipse,
      x: px(a) - r, y: py(a) - r, w: r * 2, h: r * 2,
      fill: { color: EL[a.el] || '888888' },
      line: { color: '0B1220', width: 1.2 },
      fontFace: BODY, fontSize: r > 0.2 ? 13 : 10, bold: true,
      color: a.el === 'H' || a.el === 'S' || a.el === 'F' ? '10192E' : 'FFFFFF',
      align: 'center', valign: 'middle', isTextBox: true, margin: 0,
    });
  }
}

/** Облачко неподелённой пары. */
function lonePair(s, x, y, w = 0.34, h = 0.24) {
  s.addShape(pres.ShapeType.ellipse, {
    x: x - w / 2, y: y - h / 2, w, h,
    fill: { color: '7DD3FC', transparency: 55 },
    line: { color: '7DD3FC', width: 0.75 },
  });
}

/** Подпись валентного угла. */
function angleLabel(s, x, y, text, color = TEAL) {
  s.addText(text, {
    shape: pres.ShapeType.roundRect, rectRadius: 0.06,
    x: x - 0.42, y: y - 0.16, w: 0.84, h: 0.32,
    fill: { color: '0B1220' }, line: { color, width: 1 },
    fontFace: MONO, fontSize: 11, bold: true, color,
    align: 'center', valign: 'middle', isTextBox: true, margin: 0,
  });
}

function table(s, header, rows, x, y, colW, opts = {}) {
  const headRow = header.map((c) => ({
    text: c,
    options: {
      bold: true, color: TEXT, fill: { color: '1E2A48' },
      fontSize: opts.size ?? 12, align: 'center', valign: 'middle',
    },
  }));
  const dataRows = rows.map((r, ri) => r.map((c, ci) => ({
    text: String(c),
    options: {
      color: (opts.colorFor && opts.colorFor(ri, ci)) || (ci === 0 ? TEXT : DIM),
      fill: { color: ri % 2 ? '141E36' : '0F182C' },
      fontSize: opts.size ?? 12,
      align: ci === 0 ? 'left' : 'center', valign: 'middle',
      bold: !!(opts.boldFor && opts.boldFor(ri, ci)),
    },
  })));

  s.addTable([headRow, ...dataRows], {
    x, y, colW,
    rowH: opts.rowH ?? 0.31,
    border: { type: 'solid', color: '2A3A5E', pt: 0.75 },
    fontFace: BODY,
    margin: 0.06,
  });
}

// ===========================================================================
// СЛАЙД 1 — титульный
// ===========================================================================
{
  const s = slide();
  s.addShape(pres.ShapeType.ellipse, {
    x: -2.4, y: -2.6, w: 7.6, h: 7.6,
    fill: { color: '13B39B', transparency: 88 }, line: { color: '0B1220', width: 0 },
  });
  s.addShape(pres.ShapeType.ellipse, {
    x: 9.2, y: 3.6, w: 6.4, h: 6.4,
    fill: { color: '7C5CE0', transparency: 88 }, line: { color: '0B1220', width: 0 },
  });

  s.addText('ИНДИВИДУАЛЬНЫЙ ПРОЕКТ · 11 КЛАСС', {
    x: M, y: 1.15, w: 8.4, h: 0.35,
    fontFace: BODY, fontSize: 13, bold: true, color: TEAL, charSpacing: 2,
    isTextBox: true, margin: 0,
  });
  s.addText('Компьютерное моделирование\nстроения молекул', {
    x: M, y: 1.6, w: 8.6, h: 1.9,
    fontFace: HEAD, fontSize: 42, bold: true, color: TEXT,
    lineSpacingMultiple: 1.05, isTextBox: true, margin: 0,
  });
  s.addText('Валентные углы и гибридизация вычисляются по теории Гиллеспи — программа не хранит готовых координат', {
    x: M, y: 3.6, w: 8.2, h: 0.9,
    fontFace: BODY, fontSize: 16, color: DIM, italic: true,
    lineSpacingMultiple: 1.2, isTextBox: true, margin: 0,
  });

  s.addShape(pres.ShapeType.line, {
    x: M, y: 4.75, w: 3.2, h: 0, line: { color: '2A3A5E', width: 1.5 },
  });
  s.addText('Выполнил: ученик 11 класса\nРуководитель:', {
    x: M, y: 4.95, w: 5.2, h: 0.9,
    fontFace: BODY, fontSize: 14, color: DIM, lineSpacingMultiple: 1.3,
    isTextBox: true, margin: 0,
  });
  s.addText('2026', {
    x: M, y: 6.3, w: 2, h: 0.4,
    fontFace: BODY, fontSize: 14, color: FAINT, isTextBox: true, margin: 0,
  });

  // молекула воды как эмблема
  molecule(s, 10.9, 3.5, 1.25, [
    { x: 0, y: 0.35, el: 'O', r: 0.34 },
    { x: -0.78, y: -0.36, el: 'H', r: 0.23 },
    { x: 0.78, y: -0.36, el: 'H', r: 0.23 },
  ], [[0, 1], [0, 2]], { bondWidth: 4.5 });
  lonePair(s, 10.55, 5.05, 0.44, 0.3);
  lonePair(s, 11.25, 5.05, 0.44, 0.3);
  angleLabel(s, 10.9, 3.15, '104,5°');

  s.addNotes('Здравствуйте. Тема моего проекта — компьютерное моделирование строения молекул. '
    + 'Главная особенность работы в том, что программа не показывает заранее посчитанные картинки, '
    + 'а вычисляет геометрию молекулы сама, решая физическую задачу.');
}

// ===========================================================================
// СЛАЙД 2 — зачем
// ===========================================================================
{
  const s = slide();
  title(s, 'Всё решает валентный угол', 'Две трёхатомные молекулы. Разница в строении — и совершенно разные вещества');

  card(s, M, 2.0, 5.9, 4.5);
  card(s, W - M - 5.9, 2.0, 5.9, 4.5);

  molecule(s, 3.65, 3.4, 1.15, [
    { x: 0, y: 0.3, el: 'O', r: 0.3 },
    { x: -0.75, y: -0.35, el: 'H', r: 0.2 },
    { x: 0.75, y: -0.35, el: 'H', r: 0.2 },
  ], [[0, 1], [0, 2]], { bondWidth: 4 });
  lonePair(s, 3.35, 4.55, 0.4, 0.28);
  lonePair(s, 3.95, 4.55, 0.4, 0.28);
  angleLabel(s, 3.65, 3.05, '104,5°');

  s.addText('Вода H₂O — угловая', {
    x: M + 0.3, y: 5.0, w: 5.3, h: 0.4,
    fontFace: HEAD, fontSize: 19, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'Полярна. Кипит при +100 °C, растворяет соли, образует водородные связи. Основа жизни.',
    M + 0.3, 5.45, 5.3, 0.9, { size: 14 });

  molecule(s, 9.7, 3.6, 1.15, [
    { x: -1.05, y: 0, el: 'O', r: 0.28 },
    { x: 0, y: 0, el: 'C', r: 0.3 },
    { x: 1.05, y: 0, el: 'O', r: 0.28 },
  ], [[0, 1, 2], [1, 2, 2]], { bondWidth: 4 });
  angleLabel(s, 9.7, 2.85, '180,0°');

  s.addText('Углекислый газ CO₂ — линейная', {
    x: W - M - 5.6, y: 5.0, w: 5.3, h: 0.4,
    fontFace: HEAD, fontSize: 19, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'Неполярна: диполи связей компенсируются. Газ, возгоняется при −78 °C, в воде малорастворим.',
    W - M - 5.6, 5.45, 5.3, 0.9, { size: 14 });

  s.addNotes('Почему это вообще важно. Вот две молекулы из трёх атомов. У воды угол 104,5 градуса, '
    + 'у углекислого газа — 180. Из-за этого вода полярна, а углекислый газ нет. И отсюда — вся разница '
    + 'в свойствах: одна жидкость и основа жизни, другой газ.');
}

// ===========================================================================
// СЛАЙД 3 — цель и задачи
// ===========================================================================
{
  const s = slide();
  title(s, 'Цель и задачи');

  card(s, M, 1.85, W - 2 * M, 1.35, '17284A');
  s.addText('ЦЕЛЬ', {
    x: M + 0.35, y: 2.05, w: 1.4, h: 0.3,
    fontFace: BODY, fontSize: 11.5, bold: true, color: TEAL, charSpacing: 1.6,
    isTextBox: true, margin: 0,
  });
  body(s, 'Разработать программу, которая по формуле вещества строит трёхмерную модель молекулы, вычисляет валентные углы и определяет тип гибридизации центральных атомов на основании теории Гиллеспи.',
    M + 0.35, 2.4, W - 2 * M - 0.7, 0.8, { size: 15.5, color: TEXT });

  const tasks = [
    'Изучить теорию гибридизации и теорию ОЭПВО',
    'Разработать способ ввода молекулы',
    'Реализовать расчёт неподелённых пар и стерического числа',
    'Вычислять геометрию как задачу минимизации энергии',
    'Строить трёхмерные координаты, включая молекулы с циклами',
    'Сделать наглядное отображение результата',
    'Проверить точность и найти границы применимости',
  ];
  tasks.forEach((task, k) => {
    const y = 3.5 + k * 0.53;
    badge(s, M + 0.1, y, k + 1, k < 3 ? TEAL : VIOLET, 0.38);
    body(s, task, M + 0.72, y + 0.03, 11.4, 0.4, { size: 15, color: TEXT });
  });

  s.addNotes('Цель работы и семь задач. Первые три — теоретическая подготовка и ввод данных, '
    + 'остальные — собственно вычисления, отображение и проверка.');
}

// ===========================================================================
// СЛАЙД 4 — гибридизация
// ===========================================================================
{
  const s = slide();
  title(s, 'Теория 1. Гибридизация атомных орбиталей',
    'Атомные орбитали смешиваются — образуются одинаковые по форме и энергии гибридные орбитали');

  const kinds = [
    { name: 'sp', mix: 's + p', geo: 'линейная', ang: '180°', ex: 'CO₂, C₂H₂', color: TEAL },
    { name: 'sp²', mix: 's + 2p', geo: 'плоский треугольник', ang: '120°', ex: 'BF₃, C₂H₄', color: VIOLET },
    { name: 'sp³', mix: 's + 3p', geo: 'тетраэдр', ang: '109°28′', ex: 'CH₄, NH₃, H₂O', color: AMBER },
  ];

  kinds.forEach((k, idx) => {
    const x = M + idx * 4.05;
    card(s, x, 2.0, 3.75, 4.45);

    s.addText(k.name, {
      x: x + 0.25, y: 2.22, w: 1.5, h: 0.6,
      fontFace: HEAD, fontSize: 30, bold: true, color: k.color,
      isTextBox: true, margin: 0,
    });
    body(s, k.mix, x + 1.85, 2.42, 1.7, 0.35, { size: 14, color: FAINT, align: 'right' });

    // схема расположения орбиталей
    const cx = x + 1.87, cy = 3.85;
    const dirs = idx === 0
      ? [[-1, 0], [1, 0]]
      : idx === 1
        ? [[0, -1], [-0.87, 0.5], [0.87, 0.5]]
        : [[0, -1], [-0.9, 0.42], [0.62, 0.62], [0.28, -0.2]];
    dirs.forEach(([dx, dy]) => {
      const x2 = cx + dx * 0.72, y2 = cy + dy * 0.72;
      s.addShape(pres.ShapeType.line, {
        x: Math.min(cx, x2), y: Math.min(cy, y2),
        w: Math.abs(x2 - cx), h: Math.abs(y2 - cy),
        flipH: x2 < cx, flipV: y2 < cy,
        line: { color: k.color, width: 3 },
      });
      s.addShape(pres.ShapeType.ellipse, {
        x: x2 - 0.13, y: y2 - 0.13, w: 0.26, h: 0.26,
        fill: { color: k.color, transparency: 30 }, line: { color: k.color, width: 1 },
      });
    });
    s.addShape(pres.ShapeType.ellipse, {
      x: cx - 0.2, y: cy - 0.2, w: 0.4, h: 0.4,
      fill: { color: '2A3A5E' }, line: { color: k.color, width: 1.5 },
    });

    s.addText(k.ang, {
      x: x + 0.25, y: 4.95, w: 3.25, h: 0.45,
      fontFace: MONO, fontSize: 20, bold: true, color: TEXT,
      align: 'center', isTextBox: true, margin: 0,
    });
    body(s, k.geo, x + 0.25, 5.42, 3.25, 0.35, { size: 14, color: TEXT, align: 'center' });
    body(s, k.ex, x + 0.25, 5.8, 3.25, 0.35, { size: 13, color: FAINT, align: 'center' });
  });

  s.addText('В школьном курсе тип гибридизации обычно просто заучивают для каждой молекулы. Задача программы — вычислять его.', {
    x: M, y: 6.65, w: W - 2 * M, h: 0.4,
    fontFace: BODY, fontSize: 13.5, color: FAINT, italic: true, isTextBox: true, margin: 0,
  });

  s.addNotes('Первая теоретическая опора — гибридизация. Орбитали смешиваются, и получаются '
    + 'одинаковые гибридные орбитали с определённой ориентацией. Три основных типа дают '
    + 'три характерных угла: 180, 120 и 109 градусов 28 минут.');
}

// ===========================================================================
// СЛАЙД 5 — ОЭПВО
// ===========================================================================
{
  const s = slide();
  title(s, 'Теория 2. Отталкивание электронных пар (ОЭПВО)',
    'Электронные пары вокруг центрального атома расходятся как можно дальше друг от друга');

  card(s, M, 2.05, 6.0, 2.15, '17284A');
  s.addText('Стерическое число', {
    x: M + 0.35, y: 2.25, w: 5.3, h: 0.35,
    fontFace: BODY, fontSize: 13, bold: true, color: TEAL, isTextBox: true, margin: 0,
  });
  s.addText('СЧ = σ-связи + неподелённые пары', {
    x: M + 0.35, y: 2.65, w: 5.3, h: 0.5,
    fontFace: MONO, fontSize: 18, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'Определяет и тип гибридизации, и геометрию расположения электронных пар',
    M + 0.35, 3.25, 5.3, 0.75, { size: 13.5 });

  card(s, M, 4.35, 6.0, 2.15, '17284A');
  s.addText('Число неподелённых пар', {
    x: M + 0.35, y: 4.55, w: 5.3, h: 0.35,
    fontFace: BODY, fontSize: 13, bold: true, color: VIOLET, isTextBox: true, margin: 0,
  });
  s.addText('НЭП = (n − q − m) / 2', {
    x: M + 0.35, y: 4.95, w: 5.3, h: 0.5,
    fontFace: MONO, fontSize: 18, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'n — валентные электроны, q — формальный заряд, m — сумма кратностей связей',
    M + 0.35, 5.55, 5.3, 0.75, { size: 13.5 });

  const x2 = M + 6.4;
  s.addText('Правила отталкивания', {
    x: x2, y: 2.05, w: 5.5, h: 0.4,
    fontFace: HEAD, fontSize: 19, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  s.addText('НЭП — НЭП   >   НЭП — связь   >   связь — связь', {
    x: x2, y: 2.5, w: 5.5, h: 0.45,
    fontFace: MONO, fontSize: 14, bold: true, color: AMBER, isTextBox: true, margin: 0,
  });

  bullets(s, [
    'Кратная связь — одна электронная группа',
    'Неподелённая пара «занимает больше места», чем связывающая',
    'Кратная связь занимает больше места, чем одинарная',
    'Форму молекулы определяют только ядра — пары «не видны»',
  ], x2, 3.1, 5.5, 1.6, { size: 14.5 });

  card(s, x2, 4.85, 5.5, 1.65);
  body(s, 'Пример: кислород в воде', x2 + 0.3, 5.05, 5.0, 0.3, { size: 13, color: TEAL, bold: true });
  s.addText('n = 6,  q = 0,  m = 2   →   НЭП = 2', {
    x: x2 + 0.3, y: 5.4, w: 5.0, h: 0.35,
    fontFace: MONO, fontSize: 14, color: TEXT, isTextBox: true, margin: 0,
  });
  s.addText('СЧ = 2 + 2 = 4   →   sp³,  угловая', {
    x: x2 + 0.3, y: 5.8, w: 5.0, h: 0.35,
    fontFace: MONO, fontSize: 14, color: TEXT, isTextBox: true, margin: 0,
  });

  s.addNotes('Вторая опора — теория Гиллеспи. Она позволяет предсказать тип гибридизации без '
    + 'квантовой химии. Ключевая величина — стерическое число: сумма сигма-связей и неподелённых пар. '
    + 'Число пар считается по схеме Льюиса. Для кислорода в воде получается стерическое число 4, '
    + 'то есть sp3, а форма — угловая, потому что две пары не видны.');
}

// ===========================================================================
// СЛАЙД 6 — главная идея
// ===========================================================================
{
  const s = slide();
  title(s, 'Главная идея работы',
    'Фразу «пары расходятся как можно дальше» можно понимать буквально — как задачу на минимум энергии');

  card(s, M, 2.05, W - 2 * M, 1.55, '17284A');
  s.addText('E  =  Σ  (wᵢ · wⱼ) / rᵢⱼ⁶', {
    x: M, y: 2.25, w: W - 2 * M, h: 0.75,
    fontFace: MONO, fontSize: 30, bold: true, color: TEAL,
    align: 'center', isTextBox: true, margin: 0,
  });
  s.addText('N точек размещаются на сфере так, чтобы энергия их взаимного отталкивания была минимальной',
    {
      x: M, y: 3.02, w: W - 2 * M, h: 0.4,
      fontFace: BODY, fontSize: 14, color: DIM, align: 'center', isTextBox: true, margin: 0,
    });

  s.addText('w — «размер» электронной группы. Через него в модель входит вся теория Гиллеспи:', {
    x: M, y: 3.85, w: W - 2 * M, h: 0.35,
    fontFace: BODY, fontSize: 15, color: TEXT, isTextBox: true, margin: 0,
  });

  const weights = [
    { v: '1,00', l: 'связывающая пара\nодинарной связи', c: TEAL },
    { v: '+0,18', l: 'прибавка за каждую\nединицу кратности', c: VIOLET },
    { v: '1,45', l: 'неподелённая пара —\nона «крупнее»', c: AMBER },
  ];
  weights.forEach((w, k) => {
    const x = M + k * 4.05;
    card(s, x, 4.35, 3.75, 1.45);
    s.addText(w.v, {
      x: x + 0.2, y: 4.5, w: 1.5, h: 0.6,
      fontFace: HEAD, fontSize: 27, bold: true, color: w.c, isTextBox: true, margin: 0,
    });
    body(s, w.l, x + 1.75, 4.55, 1.85, 0.9, { size: 12.5 });
  });

  s.addText('Веса подобраны так, чтобы модель воспроизводила справочные углы: H₂O — 104,5°   ·   NH₃ — 106,8°   ·   H₂C=O — 117,8°', {
    x: M, y: 6.0, w: W - 2 * M, h: 0.4,
    fontFace: BODY, fontSize: 14, color: TEXT, align: 'center', isTextBox: true, margin: 0,
  });
  s.addText('Минимум ищется градиентным спуском, 24 запуска из разных начальных положений — чтобы не попасть в локальный минимум', {
    x: M, y: 6.45, w: W - 2 * M, h: 0.4,
    fontFace: BODY, fontSize: 13, color: FAINT, italic: true, align: 'center', isTextBox: true, margin: 0,
  });

  s.addNotes('А теперь главное. Фразу «пары расходятся как можно дальше» я понял буквально — '
    + 'как физическую задачу. Точки на сфере отталкиваются по закону обратной шестой степени, '
    + 'и надо найти минимум энергии. Вся теория Гиллеспи входит сюда через один параметр — '
    + 'размер группы w. Неподелённая пара крупнее, кратная связь крупнее. Три числа — и всё.');
}

// ===========================================================================
// СЛАЙД 7 — что получается само
// ===========================================================================
{
  const s = slide();
  title(s, 'Правила учебника получаются сами',
    'Ни одно из них не запрограммировано отдельно — всё это следствия минимума энергии');

  const items = [
    {
      y: 2.0, hTitle: 0.62,
      t: '2 группы — линия, 3 — треугольник,\n4 — тетраэдр, 6 — октаэдр',
      d: 'Базовые геометрии выходят прямо из минимизации',
    },
    {
      y: 3.25, hTitle: 0.32,
      t: 'Неподелённые пары сжимают углы',
      d: 'CH₄ 109,5°  →  NH₃ 106,8°  →  H₂O 104,5°',
    },
    {
      y: 4.1, hTitle: 0.32,
      t: 'Кратная связь расширяет свой угол',
      d: 'В H₂C=O угол H–C–H сжат до 117,8°',
    },
  ];
  items.forEach((it, k) => {
    badge(s, M, it.y + 0.06, k + 1, TEAL, 0.42);
    body(s, it.t, M + 0.62, it.y, 5.4, it.hTitle, { size: 15, color: TEXT, bold: true });
    body(s, it.d, M + 0.62, it.y + it.hTitle, 5.4, 0.35, { size: 13, color: FAINT });
  });

  card(s, M, 4.95, 6.1, 1.75, '1A2440');
  body(s, 'И даже «трудные» правила для гипервалентных частиц — те, которые в учебнике даются готовым списком.',
    M + 0.3, 5.15, 5.5, 1.3, { size: 14, color: AMBER });

  // XeF4 — плоский квадрат
  const cx = 10.4, cy = 3.35;
  molecule(s, cx, cy, 1.0, [
    { x: 0, y: 0, el: 'Xe', r: 0.3 },
    { x: -1.15, y: 0, el: 'F', r: 0.2 },
    { x: 1.15, y: 0, el: 'F', r: 0.2 },
    { x: 0, y: -1.15, el: 'F', r: 0.2 },
    { x: 0, y: 1.15, el: 'F', r: 0.2 },
  ], [[0, 1], [0, 2], [0, 3], [0, 4]], { bondWidth: 3.4 });
  lonePair(s, cx - 0.72, cy - 0.72, 0.42, 0.3);
  lonePair(s, cx + 0.72, cy + 0.72, 0.42, 0.3);
  angleLabel(s, cx + 0.72, cy - 0.5, '90,0°');

  s.addText('XeF₄ — плоский квадрат', {
    x: 8.4, y: 4.85, w: 4.2, h: 0.35,
    fontFace: HEAD, fontSize: 17, bold: true, color: TEXT,
    align: 'center', isTextBox: true, margin: 0,
  });
  body(s, 'Две неподелённые пары САМИ становятся в транс-положение — по разные стороны октаэдра. Программа выводит это, а не берёт из таблицы.',
    8.15, 5.3, 4.7, 1.3, { size: 13.5, align: 'center' });

  s.addNotes('И вот что оказалось самым интересным. Все правила, которые в учебнике перечислены списком, '
    + 'получаются из этой одной задачи автоматически. Я не программировал отдельно ни тетраэдр, '
    + 'ни правило о том, что пара уходит в экватор бипирамиды, ни то, что две пары в октаэдре '
    + 'становятся напротив друг друга. Вот тетрафторид ксенона: программа сама поставила пары в транс, '
    + 'и получился плоский квадрат.');
}

// ===========================================================================
// СЛАЙД 8 — конвейер
// ===========================================================================
{
  const s = slide();
  title(s, 'Как работает программа', 'Шесть этапов от строки с формулой до трёхмерной модели');

  const steps = [
    { n: '1', t: 'Разбор ввода', d: 'Формула, название или SMILES → граф молекулы', c: TEAL },
    { n: '2', t: 'Схема Льюиса', d: 'Неподелённые пары, стерическое число, делокализация', c: TEAL },
    { n: '3', t: 'Минимизация', d: 'Направления групп из минимума энергии на сфере', c: VIOLET },
    { n: '4', t: 'Сборка', d: 'Обход графа, шаблоны циклов, выбор конформации', c: VIOLET },
    { n: '5', t: 'Релаксация', d: 'Устранение сближений, уточнение длин связей', c: AMBER },
    { n: '6', t: 'Отображение', d: 'Собственный трёхмерный рендер на Canvas 2D', c: AMBER },
  ];

  steps.forEach((st, k) => {
    const col = k % 3, row = Math.floor(k / 3);
    const x = M + col * 4.05;
    const y = 2.05 + row * 2.35;
    card(s, x, y, 3.75, 2.0);
    badge(s, x + 0.28, y + 0.28, st.n, st.c, 0.44);
    s.addText(st.t, {
      x: x + 0.85, y: y + 0.3, w: 2.7, h: 0.4,
      fontFace: HEAD, fontSize: 17, bold: true, color: TEXT, isTextBox: true, margin: 0,
    });
    body(s, st.d, x + 0.28, y + 0.95, 3.2, 0.95, { size: 13.5 });
  });

  s.addText('Химическое ядро не содержит ни одного обращения к браузеру — при переносе на C++ или Python переписывается только отображение', {
    x: M, y: 6.75, w: W - 2 * M, h: 0.4,
    fontFace: BODY, fontSize: 13, color: FAINT, italic: true, align: 'center', isTextBox: true, margin: 0,
  });

  s.addNotes('Устройство программы. Шесть этапов. Отдельно отмечу, что вся химия отделена от графики: '
    + 'вычислительное ядро не знает ничего про браузер, поэтому его можно перенести на другой язык, '
    + 'переписав только рисование.');
}

// ===========================================================================
// СЛАЙД 9 — ввод и изомерия
// ===========================================================================
{
  const s = slide();
  title(s, 'Ввод: формула не определяет строение',
    'Одной брутто-формуле может отвечать несколько разных веществ — и это превращено в возможность');

  const ways = [
    { k: 'Брутто-формула', v: 'C2H6O', d: 'покажет ВСЕ изомеры', c: TEAL },
    { k: 'Название', v: 'бензол', d: 'поиск по базе и синонимам', c: VIOLET },
    { k: 'SMILES', v: 'CC(=O)O', d: 'молекула, которой нет в базе', c: AMBER },
  ];
  ways.forEach((w, k) => {
    const y = 2.05 + k * 0.82;
    card(s, M, y, 6.0, 0.7, '17284A');
    body(s, w.k, M + 0.28, y + 0.13, 1.85, 0.4, { size: 13.5, color: w.c, bold: true });
    s.addText(w.v, {
      x: M + 2.15, y: y + 0.11, w: 1.85, h: 0.42,
      fontFace: MONO, fontSize: 15, bold: true, color: TEXT, isTextBox: true, margin: 0,
    });
    body(s, w.d, M + 4.05, y + 0.15, 1.85, 0.4, { size: 12.5 });
  });

  s.addText('C₂H₆O  →  два разных вещества', {
    x: M, y: 4.75, w: 6.0, h: 0.45,
    fontFace: HEAD, fontSize: 20, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'Программа находит все отвечающие формуле вещества и предлагает выбрать. Структурная изомерия видна наглядно.',
    M, 5.3, 6.0, 1.0, { size: 14 });

  const bx = 7.35;
  card(s, bx, 2.05, 5.28, 2.15);
  molecule(s, bx + 1.5, 3.05, 0.78, [
    { x: -1.1, y: 0.25, el: 'C', r: 0.19 },
    { x: 0, y: -0.25, el: 'C', r: 0.19 },
    { x: 1.1, y: 0.25, el: 'O', r: 0.19 },
    { x: 2.0, y: -0.2, el: 'H', r: 0.13 },
  ], [[0, 1], [1, 2], [2, 3]], { bondWidth: 2.6 });
  s.addText('Этанол', {
    x: bx + 2.9, y: 2.35, w: 2.2, h: 0.35,
    fontFace: HEAD, fontSize: 16, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'CH₃–CH₂–OH\nжидкость, кипит при +78 °C\nводородные связи', bx + 2.9, 2.75, 2.2, 1.2, { size: 12.5 });

  card(s, bx, 4.4, 5.28, 2.15);
  molecule(s, bx + 1.5, 5.4, 0.78, [
    { x: -1.1, y: 0.25, el: 'C', r: 0.19 },
    { x: 0, y: -0.25, el: 'O', r: 0.19 },
    { x: 1.1, y: 0.25, el: 'C', r: 0.19 },
  ], [[0, 1], [1, 2]], { bondWidth: 2.6 });
  s.addText('Диметиловый эфир', {
    x: bx + 2.9, y: 4.7, w: 2.2, h: 0.35,
    fontFace: HEAD, fontSize: 16, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'CH₃–O–CH₃\nгаз, кипит при −24 °C\nводородных связей нет', bx + 2.9, 5.1, 2.2, 1.2, { size: 12.5 });

  s.addNotes('О вводе. Тут есть принципиальная тонкость: брутто-формула не определяет строение. '
    + 'C2H6O — это и этанол, и диметиловый эфир: одна жидкость, другой газ. Я сделал из этого ограничения '
    + 'возможность — программа показывает все изомеры. А для однозначного задания структуры используется '
    + 'нотация SMILES, и это позволяет построить молекулу, которой в базе вообще нет.');
}

// ===========================================================================
// СЛАЙД 10 — циклы
// ===========================================================================
{
  const s = slide();
  title(s, 'Молекулы с циклами', 'Обход графа не умеет замыкать кольцо — циклы укладываются по геометрическим шаблонам');

  // бензол
  card(s, M, 2.05, 3.9, 3.85);
  const hex = [];
  for (let k = 0; k < 6; k++) {
    const a = (Math.PI / 3) * k - Math.PI / 2;
    hex.push({ x: Math.cos(a), y: Math.sin(a), el: 'C', r: 0.17 });
  }
  molecule(s, M + 1.95, 3.55, 0.72, hex,
    [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]], { bondWidth: 3 });
  s.addShape(pres.ShapeType.ellipse, {
    x: M + 1.95 - 0.38, y: 3.55 - 0.38, w: 0.76, h: 0.76,
    fill: { color: '0B1220', transparency: 100 }, line: { color: TEAL, width: 2, dashType: 'dash' },
  });
  s.addText('Бензол — плоский', {
    x: M + 0.2, y: 4.7, w: 3.5, h: 0.35,
    fontFace: HEAD, fontSize: 16, bold: true, color: TEXT, align: 'center', isTextBox: true, margin: 0,
  });
  body(s, 'Правильный шестиугольник, углы 120°, все связи C–C равны 1,39 Å',
    M + 0.2, 5.08, 3.5, 0.7, { size: 12.5, align: 'center' });

  // циклогексан-кресло
  card(s, M + 4.2, 2.05, 3.9, 3.85);
  const chair = [
    { x: -1.0, y: -0.18, el: 'C', r: 0.17 }, { x: -0.5, y: 0.35, el: 'C', r: 0.17 },
    { x: 0.5, y: 0.35, el: 'C', r: 0.17 }, { x: 1.0, y: -0.18, el: 'C', r: 0.17 },
    { x: 0.5, y: -0.62, el: 'C', r: 0.17 }, { x: -0.5, y: -0.62, el: 'C', r: 0.17 },
  ];
  molecule(s, M + 6.15, 3.55, 0.95, chair,
    [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]], { bondWidth: 3 });
  s.addText('Циклогексан — «кресло»', {
    x: M + 4.4, y: 4.7, w: 3.5, h: 0.35,
    fontFace: HEAD, fontSize: 16, bold: true, color: TEXT, align: 'center', isTextBox: true, margin: 0,
  });
  body(s, 'Единственный ненапряжённый цикл: углы 111,4° почти равны тетраэдрическому',
    M + 4.4, 5.08, 3.5, 0.7, { size: 12.5, align: 'center' });

  // циклопропан
  card(s, M + 8.4, 2.05, 3.9, 3.85);
  molecule(s, M + 10.35, 3.6, 0.85, [
    { x: 0, y: -0.75, el: 'C', r: 0.17 },
    { x: -0.65, y: 0.38, el: 'C', r: 0.17 },
    { x: 0.65, y: 0.38, el: 'C', r: 0.17 },
  ], [[0, 1], [1, 2], [2, 0]], { bondWidth: 3, bondColor: ROSE });
  angleLabel(s, M + 10.35, 3.15, '60,0°', ROSE);
  s.addText('Циклопропан — напряжён', {
    x: M + 8.6, y: 4.7, w: 3.5, h: 0.35,
    fontFace: HEAD, fontSize: 16, bold: true, color: TEXT, align: 'center', isTextBox: true, margin: 0,
  });
  body(s, 'Угол 60° вместо 109,5°: угловое напряжение около 115 кДж/моль',
    M + 8.6, 5.08, 3.5, 0.7, { size: 12.5, align: 'center' });

  bullets(s, [
    'Наименьшие циклы ищутся так: связь мысленно удаляется, и находится кратчайший путь между её концами',
    'Конденсированные кольца (нафталин, антрацен) пристраиваются по общей связи',
    'Мостиковые системы (адамантан) строит обход графа, а расхождения устраняет релаксация',
  ], M, 6.15, W - 2 * M, 1.1, { size: 13.5, gap: 4 });

  s.addNotes('Отдельная сложность — циклы. Простой обход графа не умеет замкнуть кольцо: '
    + 'последняя связь окажется неправильной длины. Поэтому циклы сначала ищутся в графе, '
    + 'а потом укладываются по шаблонам: ароматические — плоским многоугольником, '
    + 'насыщенный шестичленный — креслом. Обратите внимание на циклопропан: там угол 60 градусов, '
    + 'и это не ошибка программы, а настоящее угловое напряжение.');
}

// ===========================================================================
// СЛАЙД 11 — делокализация
// ===========================================================================
{
  const s = slide();
  title(s, 'Делокализация: формула Льюиса обманывает',
    'Программа сама обнаруживает случаи, когда выбор «где рисовать двойную связь» произволен');

  card(s, M, 2.15, 5.6, 4.3);
  s.addText('Как рисуют', {
    x: M + 0.3, y: 2.35, w: 5.0, h: 0.35,
    fontFace: BODY, fontSize: 13, bold: true, color: FAINT, isTextBox: true, margin: 0,
  });
  molecule(s, M + 2.8, 4.0, 1.0, [
    { x: 0, y: 0, el: 'N', r: 0.24 },
    { x: 0, y: -1.15, el: 'O', r: 0.22 },
    { x: -1.0, y: 0.6, el: 'O', r: 0.22 },
    { x: 1.0, y: 0.6, el: 'O', r: 0.22 },
  ], [[0, 1, 2], [0, 2], [0, 3]], { bondWidth: 3.2 });
  body(s, 'одна двойная связь и две одинарные — три равноправных варианта',
    M + 0.3, 5.55, 5.0, 0.7, { size: 13.5, align: 'center' });

  card(s, W - M - 5.6, 2.15, 5.6, 4.3, '17284A');
  s.addText('Как на самом деле', {
    x: W - M - 5.3, y: 2.35, w: 5.0, h: 0.35,
    fontFace: BODY, fontSize: 13, bold: true, color: TEAL, isTextBox: true, margin: 0,
  });
  molecule(s, W - M - 2.8, 4.0, 1.0, [
    { x: 0, y: 0, el: 'N', r: 0.24 },
    { x: 0, y: -1.15, el: 'O', r: 0.22 },
    { x: -1.0, y: 0.6, el: 'O', r: 0.22 },
    { x: 1.0, y: 0.6, el: 'O', r: 0.22 },
  ], [[0, 1], [0, 2], [0, 3]], { bondWidth: 3.2, bondColor: TEAL });
  angleLabel(s, W - M - 2.8, 3.3, '120,0°');
  body(s, 'все три связи одинаковы: 1,24 Å, порядок 1⅓, правильный треугольник',
    W - M - 5.3, 5.55, 5.0, 0.7, { size: 13.5, align: 'center', color: TEXT });

  s.addText('Правило срабатывает для нитрат-, нитрит-, карбонат-, сульфат-, ацетат-ионов, нитрогруппы, озона — и НЕ срабатывает для уксусной кислоты, где второй кислород несёт водород и потому неравноценен первому', {
    x: M, y: 6.6, w: W - 2 * M, h: 0.55,
    fontFace: BODY, fontSize: 13, color: FAINT, align: 'center',
    lineSpacingMultiple: 1.15, isTextBox: true, margin: 0,
  });

  s.addNotes('Ещё одна тонкость. Формулу нитрат-иона рисуют с одной двойной связью и двумя одинарными. '
    + 'Но на самом деле все три связи одинаковы. Программа обнаруживает такие случаи сама: '
    + 'если центральный атом связан с несколькими одинаковыми концевыми атомами, но кратности разные, '
    + 'значит выбор произволен, и их надо усреднить. Важно, что для уксусной кислоты правило '
    + 'не срабатывает — там второй кислород с водородом, он неравноценен.');
}

// ===========================================================================
// СЛАЙД 12 — результаты
// ===========================================================================
{
  const s = slide();
  title(s, 'Результат: средняя ошибка 1,7°', 'Сравнение со справочными данными структурной химии, 66 углов');

  stat(s, M, 1.9, 3.8, '1,7°', 'средняя ошибка\nв области применимости');
  stat(s, M + 4.0, 1.9, 3.8, '89', 'молекул в базе\nданных', VIOLET);
  stat(s, M + 8.0, 1.9, 3.8, '25', 'молекул построено\nвне базы', AMBER);

  table(s,
    ['Молекула', 'Угол', 'ОЭПВО', 'Опыт', 'Δ'],
    [
      ['Метан', 'H–C–H', '109,5°', '109,5°', '0,0'],
      ['Аммиак', 'H–N–H', '106,8°', '107,0°', '−0,2'],
      ['Вода', 'H–O–H', '104,5°', '104,5°', '0,0'],
      ['Углекислый газ', 'O–C–O', '180,0°', '180,0°', '0,0'],
      ['Трифторид бора', 'F–B–F', '120,0°', '120,0°', '0,0'],
      ['Нитрат-ион', 'O–N–O', '120,0°', '120,0°', '0,0'],
      ['Этилен', 'H–C–H', '117,8°', '117,4°', '+0,4'],
      ['Формальдегид', 'H–C–H', '117,8°', '116,5°', '+1,3'],
      ['Гексафторид серы', 'F–S–F', '90,0°', '90,0°', '0,0'],
      ['Тетрафторид ксенона', 'F–Xe–F', '90,0°', '90,0°', '0,0'],
    ],
    M, 3.55, [2.35, 1.1, 1.0, 0.95, 0.8],
    { size: 11.5, rowH: 0.29, colorFor: (r, c) => (c === 4 ? TEAL : null) },
  );

  card(s, M + 6.75, 3.55, 5.15, 3.2, '17284A');
  body(s, 'Для учебной модели, не использующей никаких квантовохимических расчётов, это очень хороший результат.',
    M + 7.05, 3.85, 4.55, 0.9, { size: 15, color: TEXT });
  body(s, 'Все базовые геометрии — тетраэдр, плоский треугольник, октаэдр, тригональная бипирамида — воспроизводятся точно.',
    M + 7.05, 4.85, 4.55, 1.1, { size: 14 });
  body(s, 'Сравнение проводилось по данным структурной химии: газовая электронография и микроволновая спектроскопия.',
    M + 7.05, 5.95, 4.55, 0.7, { size: 12.5, color: FAINT });

  s.addNotes('Результаты. Средняя ошибка предсказания валентного угла — 1,7 градуса по 66 сравнениям '
    + 'со справочными данными. Для модели, в которой нет ни одного квантовохимического расчёта, '
    + 'это очень хорошо. В таблице видно, что базовые геометрии воспроизводятся точно.');
}

// ===========================================================================
// СЛАЙД 13 — границы модели
// ===========================================================================
{
  const s = slide();
  title(s, 'Где теория перестаёт работать',
    'Эти расхождения не скрываются — программа помечает такие атомы и объясняет причину');

  card(s, M, 2.0, 6.05, 2.15, '2A1F14');
  s.addText('H₂S', {
    x: M + 0.35, y: 2.2, w: 1.3, h: 0.55,
    fontFace: HEAD, fontSize: 26, bold: true, color: AMBER, isTextBox: true, margin: 0,
  });
  s.addText('предсказано 104,5°   —   опыт 92,1°', {
    x: M + 1.7, y: 2.3, w: 4.1, h: 0.4,
    fontFace: MONO, fontSize: 13, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'У серы 3s- и 3p-орбитали сильно различаются по энергии, смешивать их невыгодно. Гибридизация почти не идёт, связи образуют чистые p-орбитали — а угол между ними 90°.',
    M + 0.35, 2.85, 5.4, 1.2, { size: 13.5 });

  card(s, W - M - 6.05, 2.0, 6.05, 2.15, '2A1620');
  s.addText('C₃H₆', {
    x: W - M - 5.7, y: 2.2, w: 1.3, h: 0.55,
    fontFace: HEAD, fontSize: 26, bold: true, color: ROSE, isTextBox: true, margin: 0,
  });
  s.addText('предсказано 109,5°   —   опыт 60,0°', {
    x: W - M - 4.35, y: 2.3, w: 4.1, h: 0.4,
    fontFace: MONO, fontSize: 13, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'Циклопропан: геометрия кольца жёстко задаёт 60°. Разница почти в 50° — это и есть угловое напряжение, из-за которого малые циклы легко раскрываются.',
    W - M - 5.7, 2.85, 5.4, 1.2, { size: 13.5 });

  table(s,
    ['Молекула', 'Угол', 'ОЭПВО', 'Опыт', 'Δ', 'Причина'],
    [
      ['Сероводород', 'H–S–H', '104,5°', '92,1°', '+12,4', 'слабая гибридизация'],
      ['Фосфин', 'H–P–H', '106,8°', '93,5°', '+13,3', 'слабая гибридизация'],
      ['Диметилсульфоксид', 'C–S–C', '105,8°', '96,6°', '+9,2', 'слабая гибридизация'],
      ['Циклопропан', 'C–C–C', '109,5°', '60,0°', '+49,5', 'напряжение цикла'],
      ['Циклобутан', 'C–C–C', '109,5°', '88,0°', '+21,5', 'напряжение цикла'],
      ['Трифторид азота', 'F–N–F', '106,8°', '102,2°', '+4,6', 'правило Бента'],
      ['Диметиловый эфир', 'C–O–C', '104,5°', '111,7°', '−7,2', 'объём заместителей'],
    ],
    M, 4.45, [2.5, 1.1, 1.0, 0.95, 0.85, 3.0],
    { size: 12, rowH: 0.3, colorFor: (r, c) => (c === 4 ? AMBER : null) },
  );

  s.addNotes('И, на мой взгляд, самое ценное. Программа позволяет увидеть, где теория перестаёт работать. '
    + 'Самый яркий случай — сероводород: предсказание 104,5, а на опыте 92,1. Причина известна: '
    + 'у серы гибридизация почти не идёт, и связи образуют почти чистые p-орбитали, а между ними 90 градусов. '
    + 'Второй случай — циклы: там угол задан геометрией кольца. Все такие расхождения программа помечает '
    + 'и объясняет, а не прячет.');
}

// ===========================================================================
// СЛАЙД 14 — проверка
// ===========================================================================
{
  const s = slide();
  title(s, 'Как проверялась правильность', 'Две автоматические программы-теста, запускаются одной командой');

  card(s, M, 2.05, 5.9, 4.2);
  badge(s, M + 0.32, 2.3, '1', TEAL, 0.44);
  s.addText('Вся база — 89 молекул', {
    x: M + 0.9, y: 2.32, w: 4.7, h: 0.4,
    fontFace: HEAD, fontSize: 18, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  bullets(s, [
    'разбор и построение без ошибок',
    'длины связей соответствуют сумме ковалентных радиусов',
    'атомы нигде не налезают друг на друга',
    'для 27 эталонных молекул — тип гибридизации, форма, число циклов, плоскостность кольца',
    'сравнение углов со справочными данными',
  ], M + 0.32, 2.95, 5.25, 3.1, { size: 13.5, gap: 6 });

  card(s, W - M - 5.9, 2.05, 5.9, 4.2);
  badge(s, W - M - 5.58, 2.3, '2', VIOLET, 0.44);
  s.addText('25 молекул вне базы', {
    x: W - M - 5.0, y: 2.32, w: 4.7, h: 0.4,
    fontFace: HEAD, fontSize: 18, bold: true, color: TEXT, isTextBox: true, margin: 0,
  });
  body(s, 'Аспирин, парацетамол, кофеин, антрацен, индол, пиримидин, адамантан, норборнан, глицерин, аланин, гептафторид иода…',
    W - M - 5.58, 2.95, 5.25, 1.0, { size: 13.5, color: TEXT });
  bullets(s, [
    'строятся только по строке SMILES — готовых координат нет',
    'проверяется отсутствие наложений и правильность длин связей',
    'отдельно проверяются понятные сообщения об ошибках ввода',
  ], W - M - 5.58, 4.1, 5.25, 1.9, { size: 13.5, gap: 6 });

  card(s, M, 6.45, W - 2 * M, 0.62, '13301F');
  s.addText('Все проверки проходят успешно', {
    x: M, y: 6.5, w: W - 2 * M, h: 0.52,
    fontFace: BODY, fontSize: 16, bold: true, color: '86EFAC',
    align: 'center', valign: 'middle', isTextBox: true, margin: 0,
  });

  s.addNotes('Правильность проверяется автоматически. Первый тест прогоняет всю базу — 89 молекул, '
    + 'второй строит 25 молекул, которых в базе нет вообще, только по строке SMILES. '
    + 'Это и есть доказательство, что программа считает, а не показывает заготовки.');
}

// ===========================================================================
// СЛАЙД 15 — демонстрация
// ===========================================================================
{
  const s = slide();
  s.addShape(pres.ShapeType.ellipse, {
    x: 8.6, y: -2.2, w: 8.0, h: 8.0,
    fill: { color: '13B39B', transparency: 90 }, line: { color: BG, width: 0 },
  });

  title(s, 'Демонстрация', 'Программа — один файл размером 135 КБ. Открывается двойным щелчком, без установки и без интернета');

  const demos = [
    { q: 'вода', d: 'угол 104,5°, две неподелённые пары — видно, почему угол меньше 109,5°' },
    { q: 'C2H6O', d: 'формула даёт два изомера — этанол и диметиловый эфир' },
    { q: 'XeF4', d: 'пары сами становятся в транс — получается плоский квадрат' },
    { q: 'бензол', d: 'плоское кольцо, все связи одинаковы, порядок 1½' },
    { q: 'CC(=O)Oc1ccccc1C(=O)O', d: 'аспирин — молекулы нет в базе, строится по SMILES' },
  ];

  demos.forEach((d, k) => {
    const y = 2.25 + k * 0.86;
    card(s, M, y, 11.9, 0.72, k === 4 ? '1E2A48' : CARD);
    s.addText(d.q, {
      x: M + 0.3, y: y + 0.14, w: 4.0, h: 0.45,
      fontFace: MONO, fontSize: k === 4 ? 13 : 15, bold: true, color: k === 4 ? AMBER : TEAL,
      valign: 'middle', isTextBox: true, margin: 0,
    });
    body(s, d.d, M + 4.5, y + 0.18, 7.1, 0.45, { size: 13.5, color: TEXT });
  });

  s.addText('Управление: перетаскивание — поворот, колесо — масштаб, щелчок по атому — подробный разбор', {
    x: M, y: 6.75, w: W - 2 * M, h: 0.4,
    fontFace: BODY, fontSize: 13, color: FAINT, italic: true, align: 'center', isTextBox: true, margin: 0,
  });

  s.addNotes('Сейчас я покажу программу вживую. Введу воду — посмотрим на угол и неподелённые пары. '
    + 'Потом формулу C2H6O — она даст два изомера. Потом тетрафторид ксенона — плоский квадрат. '
    + 'И в конце введу строку SMILES аспирина: этой молекулы в базе нет, программа построит её с нуля.');
}

// ===========================================================================
// СЛАЙД 16 — выводы
// ===========================================================================
{
  const s = slide();
  s.addShape(pres.ShapeType.ellipse, {
    x: -3.0, y: 2.6, w: 8.0, h: 8.0,
    fill: { color: '7C5CE0', transparency: 90 }, line: { color: BG, width: 0 },
  });
  title(s, 'Выводы');

  const conclusions = [
    {
      t: 'Цель достигнута',
      d: 'Программа по формуле строит трёхмерную модель, вычисляет валентные углы и определяет гибридизацию.',
      c: TEAL,
    },
    {
      t: 'Теория Гиллеспи работает как задача на минимум энергии',
      d: 'Средняя ошибка 1,7°. Все правила учебника — тетраэдр, сжатие углов парами, экваториальное положение пары, транс-расположение двух пар — получаются автоматически, а не программируются по отдельности.',
      c: VIOLET,
    },
    {
      t: 'Границы модели объяснимы',
      d: 'Все систематические расхождения имеют известные химические причины: слабая гибридизация у элементов 3-го периода, угловое напряжение циклов, правило Бента, объём заместителей.',
      c: AMBER,
    },
    {
      t: 'Практическая значимость',
      d: 'Пригодна для уроков химии в 11 классе и подготовки к олимпиадам: не требует установки, работает без интернета, запускается с флешки.',
      c: ROSE,
    },
  ];

  conclusions.forEach((c, k) => {
    const y = 1.75 + k * 1.28;
    s.addShape(pres.ShapeType.ellipse, {
      x: M + 0.02, y: y + 0.12, w: 0.22, h: 0.22,
      fill: { color: c.c }, line: { color: c.c, width: 0 },
    });
    body(s, c.t, M + 0.5, y, 11.5, 0.38, { size: 17, color: TEXT, bold: true, face: HEAD });
    body(s, c.d, M + 0.5, y + 0.42, 11.5, 0.85, { size: 13.5 });
  });

  s.addText('Спасибо за внимание', {
    x: M, y: 6.9, w: W - 2 * M, h: 0.45,
    fontFace: HEAD, fontSize: 19, bold: true, color: TEAL,
    align: 'center', isTextBox: true, margin: 0,
  });

  s.addNotes('Выводы. Цель достигнута. Главный результат — не сама программа, а подтверждение того, '
    + 'что теория Гиллеспи, сформулированная как задача на минимум энергии, воспроизводит валентные углы '
    + 'со средней ошибкой 1,7 градуса, причём все правила учебника получаются из неё сами. '
    + 'А все расхождения объясняются известными химическими причинами. Спасибо за внимание.');
}

// ---------------------------------------------------------------------------
const out = path.join(__dirname, 'Презентация.pptx');
pres.writeFile({ fileName: out }).then(() => console.log('готово:', out));
