import './styles.css';
import { App } from './ui/app';

/**
 * Точка входа. Программа полностью автономна: справочник элементов и молекул
 * зашит в сборку, обращений в сеть нет.
 */
declare global {
  interface Window { molview?: App }
}

function boot(): void {
  try {
    // Экземпляр доступен из консоли браузера — удобно при отладке и на защите:
    // molview.load('benzene') покажет нужную молекулу без мыши.
    window.molview = new App();
  } catch (err) {
    document.body.innerHTML =
      `<div style="padding:40px;font-family:system-ui;color:#fecdd3">
         <h2>Не удалось запустить программу</h2>
         <pre style="white-space:pre-wrap">${String((err as Error).message)}</pre>
       </div>`;
    throw err;
  }
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', boot);
} else {
  boot();
}
