/** Оформление сцены. Все цвета собраны здесь, чтобы тему можно было менять целиком. */
export const THEME = {
  background: '#080b14',
  backgroundGlow: '#141c33',
  fog: '#0a0e1a',
  accent: '#5eead4',
  accentWarm: '#fbbf24',
  accentViolet: '#a78bfa',
  lonePair: '#7dd3fc',
  orbital: '#c084fc',
  dipole: '#fbbf24',
  angleArc: '#5eead4',
  angleArcRing: '#fb7185',
  text: '#e6ecff',
  textDim: '#8b97b8',
};

// --- работа с цветом --------------------------------------------------------

export interface Rgb { r: number; g: number; b: number }

export function hexToRgb(hex: string): Rgb {
  const h = hex.replace('#', '');
  const full = h.length === 3 ? h.split('').map((c) => c + c).join('') : h;
  return {
    r: parseInt(full.slice(0, 2), 16),
    g: parseInt(full.slice(2, 4), 16),
    b: parseInt(full.slice(4, 6), 16),
  };
}

export const rgbToCss = (c: Rgb, alpha = 1): string =>
  alpha >= 1
    ? `rgb(${Math.round(c.r)}, ${Math.round(c.g)}, ${Math.round(c.b)})`
    : `rgba(${Math.round(c.r)}, ${Math.round(c.g)}, ${Math.round(c.b)}, ${alpha})`;

export function mix(a: Rgb, b: Rgb, t: number): Rgb {
  return { r: a.r + (b.r - a.r) * t, g: a.g + (b.g - a.g) * t, b: a.b + (b.b - a.b) * t };
}

export const lighten = (c: Rgb, t: number): Rgb => mix(c, { r: 255, g: 255, b: 255 }, t);
export const darken = (c: Rgb, t: number): Rgb => mix(c, { r: 8, g: 10, b: 20 }, t);

const cache = new Map<string, Rgb>();
export function rgb(hex: string): Rgb {
  let c = cache.get(hex);
  if (!c) { c = hexToRgb(hex); cache.set(hex, c); }
  return c;
}
