const initialDistance = 4.9;
const initialPitch = 0;
const halfFov = 39 * Math.PI / 360;
const unitsPerMeter = 2;

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

function desktopProjection(aspect) {
  // Preserve room for the whole scene in a narrow phone viewport.
  const fy = Math.min(1, aspect / 1.22) / Math.tan(halfFov);
  const near = 0.05;
  const far = 100;
  return new Float32Array([
    fy / aspect, 0, 0, 0,
    0, fy, 0, 0,
    0, 0, (far + near) / (near - far), -1,
    0, 0, 2 * far * near / (near - far), 0
  ]);
}

function orbitTransform(yaw, pitch, distance) {
  const x = Math.sin(yaw), z = Math.cos(yaw);
  const up = Math.sin(pitch), forward = Math.cos(pitch);
  return new Float32Array([
    z, 0, -x, 0,
    -up * x, forward, -up * z, 0,
    forward * x, up, forward * z, 0,
    distance * forward * x, distance * up, distance * forward * z, 1
  ]);
}

// Anchor the exhibit to the user's initial horizontal heading while keeping
// gravity vertical. Subsequent head motion and stereo offsets are preserved.
function anchorFromPose(m) {
  const heading = Math.hypot(m[8], m[10]) > 0.001
    ? Math.atan2(m[8], m[10]) : Math.atan2(-m[2], m[0]);
  return { x: m[12], y: m[13], z: m[14], cos: Math.cos(heading), sin: Math.sin(heading) };
}

function eyeTransform(m, anchor, out) {
  const { cos, sin } = anchor;
  for (let i = 0; i < 12; i += 4) {
    out[i] = cos * m[i] - sin * m[i + 2];
    out[i + 1] = m[i + 1];
    out[i + 2] = sin * m[i] + cos * m[i + 2];
  }
  const dx = m[12] - anchor.x, dz = m[14] - anchor.z;
  out[12] = unitsPerMeter * (cos * dx - sin * dz);
  out[13] = unitsPerMeter * (m[13] - anchor.y) + 0.15;
  out[14] = unitsPerMeter * (sin * dx + cos * dz) + 3.8;
  out[15] = 1;
  return out;
}

function controllerAxis(sources, component) {
  let strongest = 0;
  for (const source of sources) {
    const axes = source.gamepad?.axes;
    if (!axes?.length) continue;
    // xr-standard reserves slots 0/1 for touchpads and 2/3 for sticks.
    const axis = axes[(axes.length >= 4 ? 2 : 0) + component];
    if (Math.abs(axis) > Math.abs(strongest)) strongest = axis;
  }
  const magnitude = Math.abs(strongest);
  return magnitude <= 0.15 ? 0 : Math.sign(strongest) * (magnitude - 0.15) / 0.85;
}

export function createViewer(canvas, renderer, {
  getState, onXRStatus, onXRSelect, onXRAxis, onXRRotate, onXRSqueeze, onError
}) {
  const { gl } = renderer;
  let yaw = 0, pitch = initialPitch, distance = initialDistance;
  let desktopFrame = 0, lastTime = null;
  let session = null, starting = false, referenceSpace = null, layer = null, anchor = null;
  let contextLost = false, supportQuery = 0;
  const eye = new Float32Array(16);

  function elapsed(time) {
    const dt = lastTime === null ? 0 : clamp((time - lastTime) / 1000, 0, 0.1);
    lastTime = time;
    return dt;
  }

  function requestRender() {
    if (!session && !contextLost && !desktopFrame) {
      desktopFrame = window.requestAnimationFrame(renderDesktop);
    }
  }

  function renderDesktop() {
    desktopFrame = 0;
    if (session || contextLost) return;
    const rect = canvas.getBoundingClientRect();
    if (!rect.width || !rect.height) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 1.75);
    const width = Math.max(1, Math.round(rect.width * dpr));
    const height = Math.max(1, Math.round(rect.height * dpr));
    if (canvas.width !== width || canvas.height !== height) {
      canvas.width = width;
      canvas.height = height;
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    const refining = renderer.draw({
      viewport: { x: 0, y: 0, width, height },
      projection: desktopProjection(width / height),
      transform: orbitTransform(yaw, pitch, distance),
      state: getState(), progressive: true
    });
    if (refining) requestRender();
  }

  function resetView() {
    yaw = 0;
    pitch = initialPitch;
    distance = initialDistance;
    anchor = null;
    requestRender();
  }

  function zoomBy(factor) {
    distance = clamp(distance * factor, 2.8, 10);
    requestRender();
  }

  function orbit(dx, dy) {
    yaw -= dx;
    pitch = clamp(pitch + dy, -1.4, 1.4);
    requestRender();
  }

  const pointers = new Map();
  function span() {
    const [a, b] = pointers.values();
    return Math.hypot(a.x - b.x, a.y - b.y);
  }
  canvas.style.touchAction = 'none';
  canvas.addEventListener('pointerdown', event => {
    if (event.button !== 0) return;
    pointers.set(event.pointerId, { x: event.clientX, y: event.clientY });
    canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener('pointermove', event => {
    const point = pointers.get(event.pointerId);
    if (!point) return;
    const dx = event.clientX - point.x, dy = event.clientY - point.y;
    const oldSpan = pointers.size === 2 ? span() : 0;
    point.x = event.clientX;
    point.y = event.clientY;
    if (pointers.size === 1) orbit(dx * 0.006, dy * 0.006);
    if (oldSpan > 1 && span() > 1) zoomBy(oldSpan / span());
  });
  for (const type of ['pointerup', 'pointercancel', 'lostpointercapture']) {
    canvas.addEventListener(type, event => pointers.delete(event.pointerId));
  }
  canvas.addEventListener('wheel', event => {
    event.preventDefault();
    const unit = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? canvas.clientHeight : 1;
    zoomBy(Math.exp(clamp(event.deltaY * unit * 0.0012, -1, 1)));
  }, { passive: false });
  canvas.addEventListener('keydown', event => {
    switch (event.key) {
      case 'a': case 'A': orbit(-0.09, 0); break;
      case 'd': case 'D': orbit(0.09, 0); break;
      case 'ArrowUp': orbit(0, 0.09); break;
      case 'ArrowDown': orbit(0, -0.09); break;
      case '+': case '=': zoomBy(0.9); break;
      case '-': zoomBy(1 / 0.9); break;
      case 'Home': resetView(); break;
      default: return;
    }
    event.preventDefault();
  });
  new ResizeObserver(requestRender).observe(canvas);
  window.addEventListener('resize', requestRender);

  async function detectXR() {
    const id = ++supportQuery;
    if (starting || session) return;
    if (!window.isSecureContext) return onXRStatus('insecure');
    if (!navigator.xr) return onXRStatus('unsupported');
    try {
      const ok = await navigator.xr.isSessionSupported('immersive-vr');
      if (id === supportQuery && !starting && !session) onXRStatus(ok ? 'available' : 'unsupported');
    } catch {
      if (id === supportQuery && !starting && !session) onXRStatus('unsupported');
    }
  }

  function renderXR(time, frame) {
    if (frame.session !== session) return;
    session.requestAnimationFrame(renderXR);
    if (contextLost || session.visibilityState === 'hidden') {
      lastTime = null;
      return;
    }
    const pose = frame.getViewerPose(referenceSpace);
    if (!pose) return;
    anchor ??= anchorFromPose(pose.transform.matrix);
    const dt = elapsed(time);
    if (session.visibilityState === 'visible') {
      const horizontal = controllerAxis(session.inputSources, 0);
      const up = -controllerAxis(session.inputSources, 1);
      if (horizontal) onXRAxis(horizontal, dt);
      if (up) onXRRotate?.(up, dt);
    }
    const state = getState();
    gl.bindFramebuffer(gl.FRAMEBUFFER, layer.framebuffer);
    for (const view of pose.views) {
      const viewport = layer.getViewport(view);
      if (viewport) renderer.draw({
        viewport, projection: view.projectionMatrix,
        transform: eyeTransform(view.transform.matrix, anchor, eye), state, immersive: true, frameTime: time
      });
    }
  }

  function release(endedSession) {
    if (!endedSession || session !== endedSession) return;
    session = referenceSpace = layer = anchor = null;
    lastTime = null;
    requestRender();
    if (!starting) detectXR();
  }

  async function enterVR() {
    if (starting) return;
    if (session) {
      try { await session.end(); }
      catch { onError('vrExitFailed'); }
      return;
    }
    if (contextLost || !navigator.xr || !window.isSecureContext) return;
    starting = true;
    onXRStatus('starting');
    let pending;
    try {
      // Request within the button's activation, before any asynchronous setup.
      pending = await navigator.xr.requestSession('immersive-vr');
      session = pending;
      pending.addEventListener('end', () => release(pending), { once: true });
      window.cancelAnimationFrame(desktopFrame);
      desktopFrame = 0;
      await gl.makeXRCompatible();
      if (session !== pending) return;
      const space = await pending.requestReferenceSpace('local');
      if (session !== pending) return;
      referenceSpace = space;
      referenceSpace.addEventListener('reset', () => { anchor = null; });
      layer = new XRWebGLLayer(pending, gl, {
        alpha: false, depth: false, antialias: false, framebufferScaleFactor: 0.8
      });
      pending.updateRenderState({ baseLayer: layer });
      pending.addEventListener('select', () => onXRSelect());
      pending.addEventListener('squeeze', event => onXRSqueeze(event.inputSource?.handedness));
      anchor = lastTime = null;
      onXRStatus('active');
      pending.requestAnimationFrame(renderXR);
    } catch (error) {
      if (pending && session === pending) {
        try { await pending.end(); } catch { /* Already ended. */ }
        release(pending);
      }
      onError(error.name === 'NotAllowedError'
        ? 'vrPermissionDenied'
        : 'vrStartFailed');
    } finally {
      starting = false;
      if (!session) detectXR();
    }
  }

  canvas.addEventListener('webglcontextlost', event => {
    // Switching graphics adapters during makeXRCompatible can trigger this pair.
    event.preventDefault();
    contextLost = true;
    window.cancelAnimationFrame(desktopFrame);
    desktopFrame = 0;
    lastTime = null;
    if (session && !starting) {
      // The immersive layer's framebuffer was invalidated by the context loss.
      session.end().catch(() => {});
      onError('graphicsInterrupted');
    }
  });
  canvas.addEventListener('webglcontextrestored', () => {
    try {
      renderer.restore();
      contextLost = false;
      requestRender();
    } catch {
      onError('graphicsRestoreFailed');
    }
  });
  onXRStatus('checking');
  navigator.xr?.addEventListener('devicechange', detectXR);
  detectXR();
  requestRender();
  return { enterVR, resetView, zoomBy, requestRender };
}
