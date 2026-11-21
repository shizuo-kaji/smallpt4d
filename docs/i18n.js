const messages = {
  ja: {
    pageTitle: '4D レイトレーサー · smallpt4D',
    description: '4次元空間を進む光の反射・屈折・散乱を両眼で観察する WebXR パストレーサー。',
    navigation: 'ナビゲーション',
    language: '表示言語',
    help: '操作ガイド',
    heading: '4D レイトレーサー',
    introduction: '鏡とガラスを通して、4次元空間を進む光を観察する。',
    canvasLabel: '4次元パストレーシングの映像。ドラッグで視点を回転。左右の矢印キーでカメラを W 方向に移動。',
    pathTracing: '4次元パストレーシング',
    viewControls: '視点の操作',
    zoomIn: '拡大',
    zoomOut: '縮小',
    resetView: '基準視点に戻す',
    dragHint: 'ドラッグで視点を回転',
    swipeHint: 'スワイプで視点を回転',
    secondaryRays: '二次光線',
    depthReadout: '最大 {count} 回',
    controls: 'シーンとカメラの操作',
    scenes: 'シーン',
    chooseScene: 'シーンを選択',
    bounces: '二次光線の追跡回数',
    orangeLight: '橙色の光源',
    cameraW: 'カメラの W 移動量',
    rotation: 'カメラの4次元回転',
    rotationPlane: '4次元の回転平面',
    rotationAngle: '4次元の回転角度',
    close: '閉じる',
    opticsHeading: '鏡に映る4次元目',
    opticsExplanation: '各眼から出た一次光線は、視点と2次元の撮像面が張る3次元空間を進みます。超球面で反射・屈折・散乱した光線は、この空間の外にある光源にも届きます。',
    opticalExperiment: '「鏡面」の基準視点で、橙色の光源を点灯して「二次光線の追跡回数」を 0 にします。次に 1 に上げると、1回の反射で光源に届き、鏡に橙色の像が現れます。',
    browserControls: 'ブラウザーでの操作',
    orbitZoom: '視点を回す・拡大する',
    orbitZoomHow: 'ドラッグ・ホイール／ピンチ',
    homeHow: '映像内の照準ボタン',
    keyboard: '描画領域を選択中',
    keyboardW: '← → でカメラを W 方向へ移動',
    keyboardLight: 'スペースで光源を点灯・消灯',
    xrControls: 'VR コントローラー',
    sticks: 'スティック左右／上下',
    stickHow: 'カメラの W 移動／4次元回転',
    trigger: 'トリガー',
    triggerHow: '光源の点灯・消灯',
    grips: '左／右グリップ',
    gripsHow: '回転平面／シーンの切り替え',
    lightOffAction: '光源を消す',
    lightOnAction: '光源を点ける',
    vrEnter: 'VR で見る',
    vrChecking: 'VR を確認中',
    vrStarting: 'VR を起動中',
    vrExit: 'VR を終了',
    vrInsecure: 'VR を利用するには、ヘッドセットのブラウザーで HTTPS のページを開いてください。',
    vrNoWebGL: 'VR 表示には WebGL と WebXR に対応したブラウザーが必要です。',
    vrUnsupported: 'このブラウザーでは VR ヘッドセットを利用できません。WebXR 対応のヘッドセットで同じページを開くと、シーンを立体視できます。',
    webglUnavailable: 'WebGL を開始できませんでした。ブラウザーのグラフィック設定をご確認ください。',
    vrExitFailed: 'VR を終了できませんでした。ヘッドセットのシステムメニューから終了してください。',
    vrPermissionDenied: 'VR の利用が許可されませんでした。もう一度「VR で見る」を選び、利用を許可してください。',
    vrStartFailed: 'VR を開始できませんでした。ヘッドセットの接続を確認して、もう一度お試しください。',
    graphicsInterrupted: '描画が中断されました。シーンが再表示されたら、もう一度「VR で見る」を選んでください。',
    graphicsRestoreFailed: '描画を再開できませんでした。ページを再読み込みしてください。',
    hudLightOn: '橙色の光源　点灯',
    hudLightOff: '橙色の光源　消灯',
    hudCamera: 'カメラ ΔW　{w}',
    hudRotation: '{plane} 回転　{angle}°',
    hudStick: 'スティック　左右で W 移動・上下で 4D 回転',
    hudTrigger: 'トリガー　光源の点灯／消灯',
    hudLeftGrip: '左グリップ　回転平面',
    hudRightGrip: '右グリップ　シーン切替'
  },
  en: {
    pageTitle: '4D Ray Tracer · smallpt4D',
    description: 'A WebXR path tracer for exploring how light reflects, refracts and scatters in four dimensions.',
    navigation: 'Navigation',
    language: 'Display language',
    help: 'Guide',
    heading: '4D Ray Tracer',
    introduction: 'Explore how light travels in four dimensions through mirrors and glass.',
    canvasLabel: 'Four-dimensional path-traced view. Drag to orbit. Use the left and right arrow keys to move the camera along W.',
    pathTracing: '4D PATH TRACING',
    viewControls: 'View controls',
    zoomIn: 'Zoom in',
    zoomOut: 'Zoom out',
    resetView: 'Restore home view',
    dragHint: 'Drag to orbit',
    swipeHint: 'Swipe to orbit',
    secondaryRays: 'Secondary rays',
    depthReadout: 'Up to {count}',
    controls: 'Scene and camera controls',
    scenes: 'Scenes',
    chooseScene: 'Choose a scene',
    bounces: 'Bounce limit per path',
    orangeLight: 'Orange emitter',
    cameraW: 'Camera offset along W',
    rotation: '4D camera rotation',
    rotationPlane: '4D rotation plane',
    rotationAngle: '4D rotation angle',
    close: 'Close',
    opticsHeading: 'The fourth dimension in a mirror',
    opticsExplanation: 'Primary rays travel through the three-dimensional space spanned by the eye and its two-dimensional image plane. When a ray reflects, refracts or scatters at a hypersphere, it can reach an emitter outside this space.',
    opticalExperiment: 'In the Mirror scene, restore the home view and turn on the orange emitter. Set the bounce limit to 0, then raise it to 1. A single reflection reaches the emitter and its orange image appears in the mirror.',
    browserControls: 'Browser controls',
    orbitZoom: 'Orbit and zoom',
    orbitZoomHow: 'Drag to orbit. Scroll or pinch to zoom.',
    homeHow: 'Select the crosshair in the viewport.',
    keyboard: 'With the viewport focused',
    keyboardW: '← → to move along W',
    keyboardLight: 'Space to toggle the emitter',
    xrControls: 'VR controllers',
    sticks: 'Thumbstick',
    stickHow: 'Left/right moves along W. Up/down rotates in 4D.',
    trigger: 'Trigger',
    triggerHow: 'Toggle the emitter',
    grips: 'Left / right grip',
    gripsHow: 'Cycle the rotation plane / scene',
    lightOffAction: 'Turn emitter off',
    lightOnAction: 'Turn emitter on',
    vrEnter: 'Enter VR',
    vrChecking: 'Checking VR',
    vrStarting: 'Starting VR',
    vrExit: 'Exit VR',
    vrInsecure: 'Open the HTTPS version of this page in your headset browser to enter VR.',
    vrNoWebGL: 'VR requires a browser with WebGL and WebXR support.',
    vrUnsupported: 'VR is unavailable in this browser. Open this page in a WebXR-compatible headset to view the scene in stereo.',
    webglUnavailable: 'WebGL could not start. Check your browser’s graphics settings.',
    vrExitFailed: 'Use the headset’s system menu to end this VR session.',
    vrPermissionDenied: 'VR permission was declined. Select “Enter VR” again and allow headset access.',
    vrStartFailed: 'VR could not start. Check your headset connection and try again.',
    graphicsInterrupted: 'Rendering was interrupted. When the scene reappears, select “Enter VR” again.',
    graphicsRestoreFailed: 'Rendering could not resume. Reload the page to try again.',
    hudLightOn: 'Orange emitter on',
    hudLightOff: 'Orange emitter off',
    hudCamera: 'Camera ΔW  {w}',
    hudRotation: '{plane} angle  {angle}°',
    hudStick: 'Stick   ← → move W   ↑ ↓ rotate in 4D',
    hudTrigger: 'Trigger   toggle emitter',
    hudLeftGrip: 'Left grip   rotation plane',
    hudRightGrip: 'Right grip   next scene'
  }
};

const preferenceKey = 'smallpt4d.language';

export function initialLanguage() {
  try {
    const saved = localStorage.getItem(preferenceKey);
    if (saved === 'ja' || saved === 'en') return saved;
  } catch { /* Storage can be unavailable in a private browsing session. */ }
  for (const tag of navigator.languages || [navigator.language]) {
    const language = tag.toLowerCase().split('-')[0];
    if (language === 'ja' || language === 'en') return language;
  }
  return 'ja';
}

export function saveLanguage(language) {
  try {
    localStorage.setItem(preferenceKey, language);
  } catch { /* The chosen language still applies to this session. */ }
}

export function translate(language, key, values = {}) {
  const message = messages[language][key];
  return message.replace(/\{(\w+)\}/g, (_, name) => values[name]);
}

export function localizePage(language) {
  document.documentElement.lang = language;
  for (const node of document.querySelectorAll('[data-i18n]')) {
    node.textContent = translate(language, node.dataset.i18n);
  }
  for (const node of document.querySelectorAll('[data-i18n-label]')) {
    node.setAttribute('aria-label', translate(language, node.dataset.i18nLabel));
  }
  for (const node of document.querySelectorAll('[data-i18n-tooltip]')) {
    const label = translate(language, node.dataset.i18nTooltip);
    node.setAttribute('aria-label', label);
    node.title = label;
  }
  document.querySelector('meta[name="description"]').content = translate(language, 'description');
}
