import { createRenderer } from './renderer.js';
import { initialLanguage, localizePage, saveLanguage, translate } from './i18n.js';
import { sceneCopy } from './scene-copy.js';
import { createViewer } from './viewer.js';
import './diagrams.js';

const byId = id => document.getElementById(id);
const presets = [...document.querySelectorAll('[data-scene]')];
const planes = [...document.querySelectorAll('[data-plane]')];
const state = { scene: 'mirror', w: 0, angle: 0, plane: 'xw', bounces: 5, lightOn: true, language: initialLanguage() };
const canvas = byId('scene-canvas');
const wInput = byId('w'), angleInput = byId('angle'), bounceInput = byId('bounces');
const dialog = byId('guide-dialog');
let xrStatus = 'checking', noticeKey = '', viewer;
const languageButtons = [...document.querySelectorAll('[data-language]')];
const t = (key, values) => translate(state.language, key, values);

function updateRange(input, value) {
  input.value = value;
  const fraction = (input.valueAsNumber - Number(input.min)) / (Number(input.max) - Number(input.min));
  input.style.setProperty('--fill', (fraction * 100) + '%');
}

function updateControls() {
  if (xrStatus === 'active') return;
  updateRange(wInput, state.w.toFixed(2));
  const w = wInput.valueAsNumber;
  const label = (w > 0 ? '+' : '') + w.toFixed(2).replace('-', '−');
  byId('w-value').value = label;
  byId('viewport-w').textContent = 'ΔW = ' + label;
  updateRange(angleInput, Math.round(state.angle * 180 / Math.PI));
  byId('angle-value').value = angleInput.value.replace('-', '−') + '°';
  updateRange(bounceInput, state.bounces);
  byId('bounce-value').value = state.bounces;
  byId('depth-readout').textContent = t('depthReadout', {count: state.bounces});
  byId('light').checked = state.lightOn;
  byId('compare-light').textContent = t(state.lightOn ? 'lightOffAction' : 'lightOnAction');
  const scene = sceneCopy(state.language, state.scene);
  byId('scene-title').textContent = scene.title;
  byId('experiment-title').textContent = scene.heading;
  byId('experiment-detail').textContent = scene.detail;
  for (const button of presets) button.setAttribute('aria-pressed', button.dataset.scene === state.scene);
  for (const button of planes) button.setAttribute('aria-pressed', button.dataset.plane === state.plane);
  viewer?.requestRender();
}

function resetCamera() {
  state.angle = state.w = 0;
  state.plane = 'xw';
  if (xrStatus !== 'active') viewer?.resetView();
  updateControls();
}

function changeScene(scene) {
  state.scene = scene;
  resetCamera();
}

function setW(value) {
  state.w = Math.min(1.5, Math.max(-1.5, value));
  updateControls();
}

function toggleLight() {
  state.lightOn = !state.lightOn;
  updateControls();
}

function showGuide(key = '') {
  noticeKey = key;
  byId('xr-notice').textContent = key ? t(key) : '';
  byId('xr-notice').hidden = !key;
  if (!dialog.open) dialog.showModal();
}

function updateVRButton() {
  const keys = { checking: 'vrChecking', starting: 'vrStarting', active: 'vrExit' };
  byId('vr').disabled = xrStatus === 'checking' || xrStatus === 'starting';
  byId('vr-label').textContent = t(keys[xrStatus] || 'vrEnter');
}

function updateXR(status) {
  xrStatus = status;
  updateVRButton();
  updateControls();
}

function updateLanguage(language) {
  state.language = language;
  localizePage(language);
  for (const button of languageButtons) button.setAttribute('aria-pressed', button.dataset.language === language);
  for (const button of presets) {
    const scene = sceneCopy(language, button.dataset.scene);
    button.querySelector('span').textContent = scene.label;
    button.title = scene.title;
  }
  if (noticeKey) byId('xr-notice').textContent = t(noticeKey);
  if (!byId('render-error').hidden) byId('render-error').textContent = t('webglUnavailable');
  updateVRButton();
  updateControls();
}

for (const button of languageButtons) button.addEventListener('click', () => {
  const language = button.dataset.language;
  saveLanguage(language);
  updateLanguage(language);
});
byId('help').addEventListener('click', () => showGuide());
byId('dialog-close').addEventListener('click', () => dialog.close());
dialog.addEventListener('click', event => {
  if (event.target !== dialog) return;
  const r = dialog.getBoundingClientRect();
  if (event.clientX < r.left || event.clientX > r.right || event.clientY < r.top || event.clientY > r.bottom) dialog.close();
});
byId('vr').addEventListener('click', () => {
  if (xrStatus === 'available' || xrStatus === 'active') viewer.enterVR();
  else if (xrStatus === 'insecure') showGuide('vrInsecure');
  else if (!viewer) showGuide('vrNoWebGL');
  else showGuide('vrUnsupported');
});
for (const button of presets) button.addEventListener('click', () => changeScene(button.dataset.scene));
for (const button of planes) button.addEventListener('click', () => {
  state.plane = button.dataset.plane;
  updateControls();
});
wInput.addEventListener('input', () => setW(wInput.valueAsNumber));
angleInput.addEventListener('input', () => {
  state.angle = angleInput.valueAsNumber * Math.PI / 180;
  updateControls();
});
bounceInput.addEventListener('input', () => {
  state.bounces = bounceInput.valueAsNumber;
  updateControls();
});
byId('compare-light').addEventListener('click', toggleLight);
byId('light').addEventListener('change', event => {
  state.lightOn = event.currentTarget.checked;
  updateControls();
});
byId('zoom-in').addEventListener('click', () => viewer?.zoomBy(.85));
byId('zoom-out').addEventListener('click', () => viewer?.zoomBy(1/.85));
byId('camera-reset').addEventListener('click', resetCamera);
canvas.addEventListener('keydown', event => {
  if (event.key === 'ArrowRight' || event.key === 'ArrowLeft') {
    event.preventDefault();
    setW(state.w + (event.key === 'ArrowRight' ? .05 : -.05));
  } else if (event.key === ' ' && viewer) {
    event.preventDefault();
    if (!event.repeat) toggleLight();
  }
});

updateLanguage(state.language);
try {
  const renderer = createRenderer(canvas);
  viewer = createViewer(canvas, renderer, {
    getState: () => state,
    onXRStatus: updateXR,
    onXRSelect: toggleLight,
    onXRAxis: (axis, dt) => setW(state.w + axis * dt * .8),
    onXRRotate: (axis, dt) => {
      const turn = 2 * Math.PI;
      state.angle = (state.angle + axis * dt + Math.PI + turn) % turn - Math.PI;
    },
    onXRSqueeze: handedness => {
      if (handedness === 'left') {
        const next = (planes.findIndex(button => button.dataset.plane === state.plane) + 1) % planes.length;
        state.plane = planes[next].dataset.plane;
      } else {
        const next = (presets.findIndex(button => button.dataset.scene === state.scene) + 1) % presets.length;
        changeScene(presets[next].dataset.scene);
      }
    },
    onError: showGuide
  });
} catch (error) {
  console.error(error);
  byId('render-error').textContent = t('webglUnavailable');
  byId('render-error').hidden = false;
  for (const button of document.querySelectorAll('.controls button, .controls input, .camera-tools button, #compare-light')) button.disabled = true;
  updateXR('unsupported');
}
