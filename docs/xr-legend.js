import { translate } from './i18n.js';
import { sceneCopy } from './scene-copy.js';

// The sampler must be complete for every WebGL draw. The desktop view keeps
// the one-pixel texture; the canvas and full-size storage are created in VR.
export function createLegend(gl) {
  const texture = gl.createTexture();
  gl.activeTexture(gl.TEXTURE0);
  gl.bindTexture(gl.TEXTURE_2D, texture);
  for (const filter of [gl.TEXTURE_MIN_FILTER, gl.TEXTURE_MAG_FILTER]) {
    gl.texParameteri(gl.TEXTURE_2D, filter, gl.LINEAR);
  }
  for (const wrap of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T]) {
    gl.texParameteri(gl.TEXTURE_2D, wrap, gl.CLAMP_TO_EDGE);
  }
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);

  let canvas, ctx, previous, lastUpdate = -Infinity;

  function update(state, frameTime) {
    if (frameTime >= lastUpdate && frameTime - lastUpdate < 100) return;
    lastUpdate = frameTime;
    const w = (Math.round(state.w * 100) / 100).toFixed(2);
    const angle = Math.round(state.angle * 180 / Math.PI);
    const language = state.language || 'ja';
    const key = [language, state.scene, w, state.plane, angle, state.lightOn].join(' ');
    if (key === previous) return;
    if (!canvas) {
      canvas = document.createElement('canvas');
      canvas.width = 768;
      canvas.height = 224;
      ctx = canvas.getContext('2d', {alpha: false});
    }
    ctx.fillStyle = '#171a29';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = '#b9b4dc';
    ctx.font = '500 30px system-ui, sans-serif';
    const t = (key, values) => translate(language, key, values);
    ctx.fillText(sceneCopy(language, state.scene).title, 28, 42, 402);
    ctx.fillStyle = state.lightOn ? '#ace4d8' : '#a5a8b8';
    ctx.font = '500 26px system-ui, sans-serif';
    ctx.fillText(t(state.lightOn ? 'hudLightOn' : 'hudLightOff'), 454, 42, 287);

    ctx.font = '500 43px system-ui, sans-serif';
    ctx.fillStyle = '#eeeafa';
    const signedW = (Number(w) > 0 ? '+' : '') + w.replace('-', '−');
    ctx.fillText(t('hudCamera', {w: signedW}), 26, 103, 404);
    ctx.fillStyle = '#9de0d8';
    ctx.fillText(t('hudRotation', {plane: state.plane.toUpperCase(), angle}), 454, 103, 287);
    ctx.fillStyle = '#383c51';
    ctx.fillRect(26, 123, 716, 2);

    ctx.fillStyle = '#d9dbe8';
    ctx.font = '500 23px system-ui, sans-serif';
    ctx.fillText(t('hudStick'), 28, 161, 412);
    ctx.fillText(t('hudTrigger'), 454, 161, 287);
    ctx.fillText(t('hudLeftGrip'), 28, 200, 412);
    ctx.fillText(t('hudRightGrip'), 454, 200, 287);

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, texture);
    if (previous === undefined) {
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, canvas);
    } else {
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, gl.RGBA, gl.UNSIGNED_BYTE, canvas);
    }
    previous = key;
  }
  return {update, dispose: () => gl.deleteTexture(texture)};
}
