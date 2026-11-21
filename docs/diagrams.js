// Compact line icons share a 24-unit drawing grid with the application's mark.
const iconPaths = {
  mark: 'M12 1.2 22 6.8V17.2L12 22.8 2 17.2V6.8ZM2 6.8 12 12.4 22 6.8M12 12.4v10.4M12 5.5l6 3.4v6.3l-6 3.4-6-3.4V8.9ZM6 8.9l6 3.4 6-3.4M2 17.2l4-2M22 17.2l-4-2M12 1.2v4.3',
  help: 'M21 12a9 9 0 1 1-18 0 9 9 0 1 1 18 0ZM9.4 9.1a2.7 2.7 0 1 1 4.4 2.2c-1.3.9-1.8 1.3-1.8 2.9M12 17v.1',
  external: 'M7 17 18 6M8 6h10v10',
  headset: 'M5.5 5.7h13c1.4 0 2.3 1 2.6 2.5l1 6.9c.3 2-.9 3.4-2.8 3.4h-3c-1.3 0-2-1-2.6-1.8-.6-.9-1-1.3-1.7-1.3s-1.1.4-1.7 1.3c-.6.8-1.3 1.8-2.6 1.8h-3c-1.9 0-3.1-1.4-2.8-3.4l1-6.9c.3-1.5 1.2-2.5 2.6-2.5ZM6 11h3M15 11h3',
  plus: 'M12 5v14M5 12h14', minus: 'M5 12h14', close: 'M6 6l12 12M18 6 6 18',
  focus: 'M8 3H4v5m12-5h4v5M4 16v4h4m12-4v4h-4M15 12a3 3 0 1 1-6 0 3 3 0 1 1 6 0Z',
  orbit: 'M19 12a7 7 0 1 1-14 0 7 7 0 1 1 14 0ZM4.5 11.5C-4 18 5 19 14 14.5c5.5-2.8 9.2-6.4 7.4-7.5-1.1-.7-3-.2-4.2.3',
  mirror: 'M2 20h20M3 6l9 12 9-13M6 18l-2 2m7-2-2 2m7-2-2 2',
  glass: 'M3 7h18v10H3ZM1 2l7 5 4 10 8 5',
  room: 'M3 5h18v14H3ZM3 16l6-6 4 7 5-10 3 3',
  relay: 'M3 5h8m2 14h8M2 13l5-7 9 12 6-8',
  cell24: 'M12 2 21 7v10l-9 5-9-5V7ZM12 2 7 12l5 10 5-10ZM3 7l14 5-14 5m18-10L7 12l14 5',
  tir: 'M2 6h20M2 15l6-8 7 11 7-9M5.5 7l2.5-1 1.5 2',
  light: 'M15 18H9M10 21h4M8.5 14a6 6 0 1 1 7 0l-1 2h-5Z'
};
class UiIcon extends HTMLElement {
  static get observedAttributes() { return ['name']; }
  connectedCallback() { this.paint(); }
  attributeChangedCallback() { if (this.isConnected) this.paint(); }
  paint() {
    this.setAttribute('aria-hidden', 'true');
    const canvas = this.firstElementChild || this.appendChild(document.createElement('canvas'));
    canvas.width = canvas.height = 72;
    const ctx = canvas.getContext('2d');
    const name = this.getAttribute('name');
    ctx.scale(3, 3);
    ctx.strokeStyle = ctx.fillStyle = getComputedStyle(this).color;
    ctx.lineCap = ctx.lineJoin = 'round';
    ctx.lineWidth = name === 'mark' ? 1.2 : 1.65;
    const shape = new Path2D(iconPaths[name] || '');
    ctx.stroke(shape);
  }
}
customElements.define('ui-icon', UiIcon);
