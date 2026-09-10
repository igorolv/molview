import { copyFileSync, mkdirSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { defineConfig, type Plugin } from 'vite';
import { viteSingleFile } from 'vite-plugin-singlefile';

const projectDir = dirname(fileURLToPath(import.meta.url));

// Каталог, который раздаёт GitHub Pages (Settings → Pages → Source: main, /docs).
// Лежит в корне репозитория и, в отличие от molview/dist, версионируется.
const publishDir = resolve(projectDir, '..', 'docs');

// Убираем type="module"/crossorigin из итогового HTML.
// Модульные скрипты подчиняются правилам CORS, а у file:// источник «null»,
// поэтому при открытии файла с диска (особенно на iOS/Safari) они не
// выполняются и страница остаётся пустой. Классический <script> таких
// ограничений не имеет.
function classicScript(): Plugin {
  return {
    name: 'classic-script',
    enforce: 'post',
    generateBundle(_options, bundle) {
      for (const file of Object.values(bundle)) {
        if (file.type !== 'asset' || !file.fileName.endsWith('.html')) continue;
        file.source = String(file.source)
          .replace(/<script type="module"(\s+crossorigin)?>/g, '<script>')
          .replace(/\s*<link rel="modulepreload"[^>]*>/g, '');
      }
    },
  };
}

// После сборки кладём готовый однофайловый index.html в ../docs.
function publishToDocs(): Plugin {
  return {
    name: 'publish-to-docs',
    apply: 'build',
    closeBundle() {
      mkdirSync(publishDir, { recursive: true });
      copyFileSync(resolve(projectDir, 'dist', 'index.html'), resolve(publishDir, 'index.html'));
      // Отключаем обработку Jekyll: он лишний и может испортить статику.
      writeFileSync(resolve(publishDir, '.nojekyll'), '');
      this.info(`опубликовано: ${resolve(publishDir, 'index.html')}`);
    },
  };
}

// Сборка в один автономный HTML-файл: открывается двойным кликом,
// без сервера и без доступа в интернет.
export default defineConfig({
  base: './',
  plugins: [viteSingleFile(), classicScript(), publishToDocs()],
  build: {
    target: 'es2019',
    assetsInlineLimit: 100_000_000,
    cssCodeSplit: false,
    reportCompressedSize: false,
    modulePreload: false,
    rollupOptions: {
      output: {
        format: 'iife',
        inlineDynamicImports: true,
      },
    },
  },
});
