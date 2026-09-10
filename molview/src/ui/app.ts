import type { Analysis, AngleRecord, AtomAnalysis } from '../chem/types';
import { analyze } from '../chem/analyze';
import { element } from '../chem/periodic';
import { categories, database, isomersOf, moleculeFromEntry, resolveQuery, suggestions, type DbEntry } from '../chem/database';
import { parseSmiles } from '../chem/smiles';
import { formatBondOrder } from '../chem/resonance';
import { hybridExplanation, geometryFor } from '../chem/vsepr';
import { MoleculeRenderer, DEFAULT_OPTIONS, type ViewOptions, type AngleMode, type Style } from '../render/renderer';

/** Названия переключателей, значение которых логическое. */
type BooleanOption = 'showLonePairs' | 'showOrbitals' | 'showLabels' | 'showDipole' | 'autoRotate';

const $ = <T extends HTMLElement>(selector: string): T => {
  const el = document.querySelector<T>(selector);
  if (!el) throw new Error(`Не найден элемент ${selector}`);
  return el;
};

export class App {
  private renderer: MoleculeRenderer;
  private analysis: Analysis | null = null;
  private currentEntry: DbEntry | null = null;
  private options: ViewOptions = { ...DEFAULT_OPTIONS };
  private showHydrogenAngles = false;
  private suggestIndex = -1;
  private suggestItems: DbEntry[] = [];

  private input = $<HTMLInputElement>('#query');
  private suggestBox = $<HTMLDivElement>('#suggest');
  private panel = $<HTMLElement>('#panel');
  private heading = $<HTMLElement>('#stage-heading');
  private legend = $<HTMLElement>('#legend');

  constructor() {
    this.renderer = new MoleculeRenderer($<HTMLCanvasElement>('#scene'));
    this.renderer.onSelectChange = () => { this.renderPanel(); this.syncCards(); };
    this.renderer.onHoverChange = () => this.syncCards();

    this.buildCatalog();
    this.bindSearch();
    this.bindToolbar();
    this.bindMisc();

    window.addEventListener('resize', () => this.renderer.resize());
    this.renderer.start();

    this.load('water');
  }

  // =========================================================================
  // Загрузка молекулы
  // =========================================================================

  load(id: string): void {
    const entry = database().find((e) => e.id === id);
    if (!entry) return;
    this.showEntry(entry);
  }

  private showEntry(entry: DbEntry): void {
    this.currentEntry = entry;
    this.customSmiles = null;
    const analysis = analyze(moleculeFromEntry(entry));
    this.apply(analysis);
    this.input.value = entry.name;
    this.markCatalog(entry.id);
  }

  private customSmiles: string | null = null;

  private showSmiles(smiles: string): void {
    this.currentEntry = null;
    this.customSmiles = smiles;
    const mol = parseSmiles(smiles);
    this.apply(analyze(mol));
    this.markCatalog(null);
  }

  private apply(analysis: Analysis): void {
    this.analysis = analysis;
    const heavy = analysis.molecule.atoms.filter((a) => a.el !== 'H').length;
    this.showHydrogenAngles = heavy <= 4;
    this.renderer.setAnalysis(analysis);
    this.renderHeading();
    this.renderLegend();
    this.renderPanel();
  }

  // =========================================================================
  // Каталог
  // =========================================================================

  private buildCatalog(): void {
    const host = $<HTMLElement>('#catalog');
    host.innerHTML = '';
    for (const group of categories()) {
      const box = document.createElement('div');
      box.className = 'cat-group';
      const head = document.createElement('div');
      head.className = 'cat-head';
      head.textContent = group.name;
      box.appendChild(head);

      for (const entry of group.entries) {
        const item = document.createElement('div');
        item.className = 'cat-item';
        item.dataset.id = entry.id;
        item.innerHTML = `<span>${escapeHtml(entry.name)}</span><span class="f">${escapeHtml(entry.formulaPretty)}</span>`;
        item.addEventListener('click', () => this.showEntry(entry));
        box.appendChild(item);
      }
      host.appendChild(box);
    }
  }

  private markCatalog(id: string | null): void {
    document.querySelectorAll<HTMLElement>('.cat-item').forEach((el) => {
      const active = el.dataset.id === id;
      el.classList.toggle('active', active);
      if (active) el.scrollIntoView({ block: 'nearest' });
    });
  }

  // =========================================================================
  // Поиск
  // =========================================================================

  private bindSearch(): void {
    this.input.addEventListener('input', () => this.updateSuggestions());
    this.input.addEventListener('focus', () => this.updateSuggestions());

    this.input.addEventListener('blur', () => {
      window.setTimeout(() => { this.suggestBox.hidden = true; }, 140);
    });

    this.input.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
        e.preventDefault();
        const delta = e.key === 'ArrowDown' ? 1 : -1;
        this.suggestIndex = Math.max(-1, Math.min(this.suggestItems.length - 1, this.suggestIndex + delta));
        this.highlightSuggestion();
      } else if (e.key === 'Enter') {
        e.preventDefault();
        if (this.suggestIndex >= 0 && this.suggestItems[this.suggestIndex]) {
          this.showEntry(this.suggestItems[this.suggestIndex]);
          this.suggestBox.hidden = true;
          this.input.blur();
        } else {
          this.submit();
        }
      } else if (e.key === 'Escape') {
        this.suggestBox.hidden = true;
        this.input.blur();
      }
    });

    $('#search-go').addEventListener('click', () => this.submit());
  }

  private updateSuggestions(): void {
    const value = this.input.value.trim();
    this.suggestItems = value.length >= 1 ? suggestions(value) : [];
    this.suggestIndex = -1;

    if (this.suggestItems.length === 0) { this.suggestBox.hidden = true; return; }

    this.suggestBox.innerHTML = '';
    this.suggestItems.forEach((entry, index) => {
      const row = document.createElement('div');
      row.className = 'suggest-item';
      row.innerHTML =
        `<span class="suggest-name">${escapeHtml(entry.name)}</span>` +
        `<span class="suggest-formula">${escapeHtml(entry.formulaPretty)}</span>` +
        `<span class="suggest-cat">${escapeHtml(entry.cat)}</span>`;
      row.addEventListener('mousedown', (e) => {
        e.preventDefault();
        this.showEntry(entry);
        this.suggestBox.hidden = true;
      });
      row.addEventListener('mouseenter', () => { this.suggestIndex = index; this.highlightSuggestion(); });
      this.suggestBox.appendChild(row);
    });
    this.suggestBox.hidden = false;
  }

  private highlightSuggestion(): void {
    [...this.suggestBox.children].forEach((el, index) => {
      el.classList.toggle('active', index === this.suggestIndex);
    });
  }

  private submit(): void {
    const query = this.input.value.trim();
    if (!query) return;
    this.suggestBox.hidden = true;

    const result = resolveQuery(query);
    if (result.kind === 'entries') {
      this.showEntry(result.entries[0]);
      if (result.entries.length > 1 && result.reason === 'formula') {
        this.flashIsomerNotice(result.entries);
      }
    } else if (result.kind === 'smiles') {
      try {
        this.showSmiles(result.smiles);
      } catch (err) {
        this.showError('Не удалось построить молекулу', (err as Error).message);
      }
    } else {
      this.showError(result.message, result.hint);
    }
  }

  private flashIsomerNotice(entries: DbEntry[]): void {
    const box = document.createElement('div');
    box.className = 'note rise';
    box.style.margin = '10px 16px 0';
    box.innerHTML =
      `<b>Формуле ${escapeHtml(entries[0].formulaPretty)} отвечает ${entries.length} ${plural(entries.length, 'вещество', 'вещества', 'веществ')}.</b> ` +
      'Брутто-формула не определяет строение однозначно — это и есть структурная изомерия. ' +
      'Выберите изомер в разделе ниже.';
    this.panel.insertBefore(box, this.panel.firstChild);
  }

  private showError(message: string, hint?: string): void {
    this.panel.innerHTML =
      `<div class="error-box rise"><b>${escapeHtml(message)}</b>` +
      (hint ? `<div style="margin-top:7px">${escapeHtml(hint)}</div>` : '') +
      '</div>';
  }

  // =========================================================================
  // Панель управления сценой
  // =========================================================================

  private bindToolbar(): void {
    document.querySelectorAll<HTMLButtonElement>('[data-style]').forEach((btn) => {
      btn.addEventListener('click', () => {
        document.querySelectorAll('[data-style]').forEach((b) => b.classList.remove('active'));
        btn.classList.add('active');
        this.options.style = btn.dataset.style as Style;
        this.renderer.setOptions({ style: this.options.style });
      });
    });

    document.querySelectorAll<HTMLButtonElement>('[data-toggle]').forEach((btn) => {
      btn.addEventListener('click', () => {
        const key = btn.dataset.toggle as BooleanOption;
        const next = !this.options[key];
        this.options[key] = next;
        btn.classList.toggle('active', next);
        this.renderer.setOptions({ [key]: next });
        this.renderLegend();
      });
    });

    document.querySelectorAll<HTMLButtonElement>('[data-angles]').forEach((btn) => {
      btn.addEventListener('click', () => {
        document.querySelectorAll('[data-angles]').forEach((b) => b.classList.remove('active'));
        btn.classList.add('active');
        this.options.showAngles = btn.dataset.angles as AngleMode;
        this.renderer.setOptions({ showAngles: this.options.showAngles });
      });
    });

    $('#btn-reset').addEventListener('click', () => this.renderer.resetView());

    $('#btn-png').addEventListener('click', () => {
      const link = document.createElement('a');
      const name = this.currentEntry?.id ?? this.analysis?.formula ?? 'molecule';
      link.download = `${name}.png`;
      link.href = this.renderer.toDataUrl();
      link.click();
    });
  }

  private bindMisc(): void {
    const help = $<HTMLElement>('#help');
    $('#btn-help').addEventListener('click', () => { help.hidden = false; });
    $('#help-close').addEventListener('click', () => { help.hidden = true; });
    help.addEventListener('click', (e) => { if (e.target === help) help.hidden = true; });

    $('#btn-random').addEventListener('click', () => {
      const all = database();
      this.showEntry(all[Math.floor(Math.random() * all.length)]);
    });

    window.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') help.hidden = true;
      if (e.target === this.input) return;
      if (e.key === '/') { e.preventDefault(); this.input.focus(); this.input.select(); }
    });
  }

  // =========================================================================
  // Заголовок и легенда
  // =========================================================================

  private renderHeading(): void {
    const a = this.analysis;
    if (!a) return;
    const main = mainCenter(a);
    const name = this.currentEntry?.name ?? 'Построено по SMILES';
    const subtitle = this.currentEntry ? '' : `<div class="h-smiles">${escapeHtml(this.customSmiles ?? '')}</div>`;
    const geometry = main
      ? `${main.molecularGeom}, ${main.hybrid}-гибридизация центрального атома`
      : 'двухатомная частица';

    this.heading.innerHTML =
      `<div class="h-name">${escapeHtml(name)}</div>` +
      `<div class="h-formula">${escapeHtml(a.formulaHtml)}</div>` +
      subtitle +
      `<div class="h-geom">${escapeHtml(geometry)}</div>`;
  }

  private renderLegend(): void {
    const rows: string[] = [];
    rows.push(row('#5eead4', 'валентный угол'));
    if (this.analysis?.rings.length) rows.push(row('#fb7185', 'угол задан циклом'));
    if (this.options.showLonePairs) rows.push(row('#7dd3fc', 'неподелённая пара'));
    if (this.options.showOrbitals) rows.push(row('#c084fc', 'гибридная орбиталь'));
    if (this.options.showDipole) rows.push(row('#fbbf24', 'дипольный момент'));
    this.legend.innerHTML = rows.join('');

    function row(color: string, text: string): string {
      return `<div class="legend-row"><span class="legend-swatch" style="background:${color}"></span>${text}</div>`;
    }
  }

  // =========================================================================
  // Правая панель
  // =========================================================================

  private renderPanel(): void {
    const a = this.analysis;
    if (!a) return;
    const parts: string[] = [];

    // --- шапка ---
    const name = this.currentEntry?.name ?? 'Молекула из SMILES';
    parts.push('<div class="panel-head rise">');
    parts.push(`<div class="panel-title">${escapeHtml(name)}</div>`);
    if (!this.currentEntry && this.customSmiles) {
      parts.push(`<div style="font-family:var(--mono);font-size:12px;color:var(--text-faint);margin-top:2px">${escapeHtml(this.customSmiles)}</div>`);
    }
    parts.push(`<div class="panel-formula">${escapeHtml(a.formulaHtml)}</div>`);
    parts.push('<div class="panel-meta">');
    if (this.currentEntry) parts.push(`<span class="chip">${escapeHtml(this.currentEntry.cat)}</span>`);
    parts.push(`<span class="chip">M = ${a.mass.toFixed(2)} г/моль</span>`);
    for (const [hybrid, count] of Object.entries(a.hybridSummary).sort()) {
      parts.push(`<span class="chip violet">${hybrid} × ${count}</span>`);
    }
    if (a.rings.length > 0) {
      const aromatic = a.rings.filter((r) => r.aromatic).length;
      parts.push(`<span class="chip rose">${a.rings.length} ${plural(a.rings.length, 'цикл', 'цикла', 'циклов')}${aromatic ? `, ${aromatic} ароматических` : ''}</span>`);
    }
    parts.push(`<span class="chip ${a.polar ? 'amber' : 'accent'}">${a.polar ? 'полярная' : 'неполярная'}</span>`);
    parts.push('</div></div>');

    // --- пояснение из справочника ---
    if (this.currentEntry?.note) {
      parts.push(`<div class="section"><div class="note">${escapeHtml(this.currentEntry.note)}</div></div>`);
    }

    // --- атомы ---
    const focus = this.renderer.selectedAtom;
    const shown = a.atoms.filter((atom) =>
      focus !== null ? atom.id === focus : atom.isCentral && atom.el !== 'H');

    if (shown.length > 0) {
      parts.push('<div class="section"><div class="section-title">');
      parts.push(focus !== null ? 'Выбранный атом' : 'Центральные атомы');
      parts.push('</div>');
      const groups = focus !== null ? shown : dedupeByShape(shown, a);
      groups.forEach((atom) => parts.push(this.atomCard(atom, a)));
      if (focus !== null) {
        parts.push('<button class="isomer-btn" data-clear-selection>← ко всем центрам</button>');
      }
      parts.push('</div>');
    }

    // --- углы ---
    parts.push(this.anglesSection(a, focus));

    // --- диполь ---
    parts.push(this.dipoleSection(a));

    // --- изомеры ---
    if (this.currentEntry) {
      const isomers = isomersOf(this.currentEntry);
      if (isomers.length > 0) {
        parts.push('<div class="section"><div class="section-title">Изомеры — та же формула, другое строение</div><div class="isomers">');
        for (const iso of isomers) {
          parts.push(`<button class="isomer-btn" data-isomer="${iso.id}">${escapeHtml(iso.name)}</button>`);
        }
        parts.push('</div></div>');
      }
    }

    this.panel.innerHTML = parts.join('');
    this.bindPanelEvents();
  }

  private atomCard(atom: AtomAnalysis, a: Analysis): string {
    const info = element(atom.el);
    const active = this.renderer.selectedAtom === atom.id || this.renderer.hoveredAtom === atom.id;
    const geometry = geometryFor(atom.sigma, atom.lonePairs);

    const rows: [string, string, boolean][] = [
      ['σ-связей', String(atom.sigma), true],
      ['Неподелённых пар', String(atom.lonePairs), true],
      ['Стерическое число', String(atom.steric), true],
      ['Тип по Гиллеспи', atom.axe, false],
      ['Геометрия электронных пар', geometry.electronGeom, false],
      ['Форма молекулы', geometry.molecularGeom, false],
    ];
    if (atom.idealAngle !== null) rows.push(['Идеальный угол', `${atom.idealAngle}°`, true]);
    if (atom.predictedAngle !== null) rows.push(['Предсказанный угол', `${atom.predictedAngle.toFixed(1)}°`, true]);

    const bonds = a.molecule.bonds
      .filter((b) => b.a === atom.id || b.b === atom.id)
      .filter((b) => b.order !== 1);
    const orderNote = bonds.length > 0
      ? `<div class="formula-line">Порядок связей: ${bonds
          .map((b) => {
            const other = b.a === atom.id ? b.b : b.a;
            return `${atom.el}–${a.molecule.atoms[other].el} = <b>${formatBondOrder(b.order)}</b>`;
          })
          .filter((s, i, arr) => arr.indexOf(s) === i)
          .join(', ')}</div>`
      : '';

    return `
      <div class="atom-card rise ${active ? 'active' : ''}" data-atom="${atom.id}">
        <div class="atom-head">
          <div class="atom-badge" style="background:${info.color}">${escapeHtml(atom.el)}</div>
          <div>
            <div class="atom-title">${escapeHtml(capitalize(info.name))}${chargeSuffix(atom.charge)}</div>
            <div class="atom-sub">атом №${atom.id + 1} · ${escapeHtml(geometry.hint.split(':')[0])}</div>
          </div>
          <div class="atom-hyb">${escapeHtml(atom.hybrid)}</div>
        </div>
        <dl class="atom-grid">
          ${rows.map(([k, v, num]) => `<dt>${escapeHtml(k)}</dt><dd class="${num ? 'num' : ''}">${escapeHtml(v)}</dd>`).join('')}
        </dl>
        <div class="formula-line">
          СЧ = ${atom.sigma} σ-связ${ending(atom.sigma)} + ${atom.lonePairs} пар${pairEnding(atom.lonePairs)} = <b>${atom.steric}</b>
          → <b>${escapeHtml(atom.hybrid)}</b>${atom.steric >= 2 && atom.steric <= 7 ? `: ${escapeHtml(hybridExplanation(atom.steric))}` : ''}
        </div>
        ${orderNote}
        ${atom.warning ? `<div class="warn">${escapeHtml(atom.warning)}</div>` : ''}
      </div>`;
  }

  private anglesSection(a: Analysis, focus: number | null): string {
    let records = a.angles;
    if (focus !== null) records = records.filter((r) => r.center === focus);
    if (!this.showHydrogenAngles) {
      const heavyOnly = records.filter(
        (r) => a.molecule.atoms[r.i].el !== 'H' && a.molecule.atoms[r.j].el !== 'H',
      );
      if (heavyOnly.length > 0) records = heavyOnly;
    }
    if (records.length === 0) {
      return '<div class="section"><div class="section-title">Валентные углы</div>' +
        '<div class="empty" style="padding:14px 0">В этой частице нет ни одного атома с двумя связями — валентного угла не существует.</div></div>';
    }

    // одинаковые по смыслу углы (например, шесть углов H-C-H в этане)
    // сворачиваем в одну строку с указанием кратности
    const grouped = new Map<string, { record: AngleRecord; count: number }>();
    for (const r of [...records].sort((x, y) => x.center - y.center || x.actual - y.actual)) {
      const key = [
        a.molecule.atoms[r.i].el, r.center, a.molecule.atoms[r.j].el,
        r.actual.toFixed(1), r.predicted?.toFixed(1), r.experimental,
      ].join('|');
      const existing = grouped.get(key);
      if (existing) existing.count++;
      else grouped.set(key, { record: r, count: 1 });
    }
    const rows = [...grouped.values()].map(({ record, count }) => this.angleRow(record, a, count)).join('');

    const hasExperiment = records.some((r) => r.experimental !== null);
    return `
      <div class="section">
        <div class="section-title">Валентные углы${focus !== null ? ' выбранного атома' : ''}</div>
        <table class="angles">
          <thead><tr>
            <th>Угол</th><th>ОЭПВО</th><th>модель</th>${hasExperiment ? '<th>опыт</th><th>Δ</th>' : ''}
          </tr></thead>
          <tbody>${rows}</tbody>
        </table>
        <div style="margin-top:8px;font-size:11px;color:var(--text-faint);line-height:1.5">
          «ОЭПВО» — предсказание теории для изолированного атома,
          «модель» — угол, измеренный по построенной трёхмерной структуре${hasExperiment ? ', «опыт» — справочные данные' : ''}.
        </div>
      </div>`;
  }

  private angleRow(r: AngleRecord, a: Analysis, count: number): string {
    const multiplier = count > 1
      ? ` <span style="color:var(--text-faint)">×${count}</span>`
      : '';
    const label =
      `${a.molecule.atoms[r.i].el}–${a.molecule.atoms[r.center].el}<sub>${r.center + 1}</sub>–${a.molecule.atoms[r.j].el}${multiplier}`;
    const cells: string[] = [
      `<td>${label}</td>`,
      `<td class="muted">${r.predicted !== null ? r.predicted.toFixed(1) + '°' : '—'}</td>`,
      `<td>${r.actual.toFixed(1)}°</td>`,
    ];
    if (r.experimental !== null) {
      const delta = r.actual - r.experimental;
      const cls = Math.abs(delta) < 3 ? 'good' : Math.abs(delta) < 8 ? 'warnv' : 'bad';
      cells.push(`<td>${r.experimental.toFixed(1)}°</td>`);
      cells.push(`<td class="${cls}">${delta >= 0 ? '+' : ''}${delta.toFixed(1)}</td>`);
    } else if (a.angles.some((x) => x.experimental !== null)) {
      cells.push('<td class="muted">—</td><td class="muted">—</td>');
    }
    return `<tr class="${r.inRing ? 'ring' : ''}" data-angle-center="${r.center}">${cells.join('')}</tr>`;
  }

  private dipoleSection(a: Analysis): string {
    const measured = this.currentEntry?.dipole;
    // Расчёт надёжно отвечает на вопрос «полярна или нет» (это следствие
    // симметрии), но не даёт значения в дебаях. Поэтому шкала — безразмерная,
    // а число приводится по справочнику.
    const level = a.dipoleValue < 0.06 ? 0 : Math.min(1, a.dipoleValue / 3.2);
    const wording =
      a.dipoleValue < 0.06 ? 'диполя нет'
        : a.dipoleValue < 0.9 ? 'слабо полярная'
          : a.dipoleValue < 2.0 ? 'умеренно полярная'
            : 'сильно полярная';

    return `
      <div class="section">
        <div class="section-title">Полярность</div>
        <div style="font-size:12.5px;line-height:1.6;color:var(--text-dim)">
          ${a.polar
            ? 'Диполи связей и неподелённых пар <b style="color:var(--text)">не компенсируют</b> друг друга — молекула полярна.'
            : 'Диполи связей расположены симметрично и <b style="color:var(--text)">полностью компенсируются</b> — молекула неполярна.'}
        </div>
        <div style="margin-top:10px;display:grid;grid-template-columns:auto 1fr auto;gap:6px 10px;align-items:center;font-size:12px">
          <span style="color:var(--text-faint)">расчёт</span>
          <span class="bar"><i style="width:${(level * 100).toFixed(1)}%"></i></span>
          <span>${wording}</span>
          ${measured !== undefined ? `
            <span style="color:var(--text-faint)">эксперимент</span>
            <span class="bar"><i style="width:${Math.min(100, (measured / 4) * 100).toFixed(1)}%;background:var(--amber)"></i></span>
            <span style="font-family:var(--mono)">μ = ${measured.toFixed(2)} Д</span>` : ''}
        </div>
        <div style="margin-top:8px;font-size:11px;color:var(--text-faint);line-height:1.5">
          Программа складывает векторы диполей связей (по разности электроотрицательностей)
          и вкладов неподелённых пар. Такой расчёт надёжно отвечает на вопрос
          «полярна молекула или нет» — это следствие симметрии, — но не даёт значения
          в дебаях: точный расчёт требует квантовохимических методов.
          ${measured !== undefined ? 'Измеренное значение приведено по справочнику.' : ''}
        </div>
      </div>`;
  }

  private bindPanelEvents(): void {
    this.panel.querySelectorAll<HTMLElement>('[data-atom]').forEach((card) => {
      const id = Number(card.dataset.atom);
      card.addEventListener('click', () => this.renderer.setSelected(id));
      card.addEventListener('mouseenter', () => this.renderer.setHovered(id));
      card.addEventListener('mouseleave', () => this.renderer.setHovered(null));
    });

    this.panel.querySelectorAll<HTMLElement>('[data-angle-center]').forEach((row) => {
      const id = Number(row.dataset.angleCenter);
      row.addEventListener('click', () => this.renderer.setSelected(id));
      row.addEventListener('mouseenter', () => this.renderer.setHovered(id));
      row.addEventListener('mouseleave', () => this.renderer.setHovered(null));
    });

    this.panel.querySelectorAll<HTMLElement>('[data-isomer]').forEach((btn) => {
      btn.addEventListener('click', () => this.load(btn.dataset.isomer!));
    });

    this.panel.querySelector('[data-clear-selection]')?.addEventListener('click', () => {
      this.renderer.setSelected(null);
    });
  }

  /** Подсветка карточек в такт наведению на сцене — без полной перерисовки панели. */
  private syncCards(): void {
    const active = this.renderer.selectedAtom ?? this.renderer.hoveredAtom;
    this.panel.querySelectorAll<HTMLElement>('[data-atom]').forEach((card) => {
      card.classList.toggle('active', Number(card.dataset.atom) === active);
    });
  }
}

// ===========================================================================
// вспомогательные функции
// ===========================================================================

/** Главный центральный атом — по нему подписывается заголовок сцены. */
function mainCenter(a: Analysis): AtomAnalysis | undefined {
  return a.atoms
    .filter((x) => x.isCentral && x.el !== 'H')
    .sort((x, y) => y.sigma + y.lonePairs - (x.sigma + x.lonePairs))[0];
}

/** Одинаковые по окружению атомы показываем один раз — иначе бензол даёт шесть копий. */
function dedupeByShape(atoms: AtomAnalysis[], a: Analysis): AtomAnalysis[] {
  const seen = new Map<string, AtomAnalysis>();
  for (const atom of atoms) {
    const neighbors = a.molecule.bonds
      .filter((b) => b.a === atom.id || b.b === atom.id)
      .map((b) => {
        const other = b.a === atom.id ? b.b : b.a;
        return `${a.molecule.atoms[other].el}${b.order.toFixed(2)}`;
      })
      .sort()
      .join(',');
    const key = `${atom.el}|${atom.charge}|${atom.axe}|${neighbors}`;
    if (!seen.has(key)) seen.set(key, atom);
  }
  return [...seen.values()];
}

const escapeHtml = (s: string): string =>
  s.replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]!));

const capitalize = (s: string): string => s.charAt(0).toUpperCase() + s.slice(1);

function plural(n: number, one: string, few: string, many: string): string {
  const mod100 = n % 100;
  const mod10 = n % 10;
  if (mod100 >= 11 && mod100 <= 14) return many;
  if (mod10 === 1) return one;
  if (mod10 >= 2 && mod10 <= 4) return few;
  return many;
}

const ending = (n: number): string => plural(n, 'ь', 'и', 'ей');
const pairEnding = (n: number): string => plural(n, 'а', 'ы', '');

function chargeSuffix(charge: number): string {
  if (Math.abs(charge) < 0.01) return '';
  if (!Number.isInteger(charge)) return ` <span style="color:var(--text-faint)">(q = ${charge.toFixed(2)})</span>`;
  const sign = charge > 0 ? '+' : '−';
  const magnitude = Math.abs(charge) > 1 ? Math.abs(charge) : '';
  return ` <span style="color:${charge > 0 ? '#fca5a5' : '#93c5fd'}">${magnitude}${sign}</span>`;
}
