import { SCENE_TITLES } from './optics.js';

// The button names stay short enough for the three-column selector on phones.
const copy = {
  mirror: {
    ja: {
      label: '鏡面', title: '鏡面の超球',
      heading: '鏡だけに映る光源',
      detail: '基準視点で橙色の光源を消すと、鏡の中の像が消えます。光源を点け、追跡回数を 0 から 1 にすると、鏡で1回反射した光が届きます。'
    },
    en: {
      label: 'Mirror', title: SCENE_TITLES.mirror,
      heading: 'A lamp visible in the mirror',
      detail: 'At the home view, switch off the orange lamp and its reflection disappears. Turn it back on, then raise the bounce limit from 0 to 1 to see the light after one reflection.'
    }
  },
  glass: {
    ja: {
      label: 'ガラス', title: 'ガラスの超球レンズ',
      heading: 'ガラスを通る光',
      detail: '橙色の光源を点け、追跡回数を 1 から 2 に上げてください。中心の光線がガラスへの入射と出射を経て、直接見えない光源に届きます。'
    },
    en: {
      label: 'Glass', title: SCENE_TITLES.glass,
      heading: 'Through two glass interfaces',
      detail: 'With the orange lamp on, raise the bounce limit from 1 to 2. The central ray refracts as it enters the glass and again as it leaves, reaching a lamp hidden from direct view.'
    }
  },
  room: {
    ja: {
      label: '光の部屋', title: '4次元の部屋',
      heading: 'W 方向にも壁がある部屋',
      detail: '左右と奥、床と天井に加え、W の両側にも鏡の壁があります。カメラを W 方向に移動し、追跡回数を増やして壁や超球に映る像を比べてください。'
    },
    en: {
      label: 'Light room', title: SCENE_TITLES.room,
      heading: 'Walls on both sides of W',
      detail: 'The room has mirrored walls in both W directions as well as a floor, ceiling, back and side walls. Move the camera along W and increase the bounce limit to compare reflections in the walls and spheres.'
    }
  },
  relay: {
    ja: {
      label: '二段鏡', title: '二枚鏡の潜望鏡',
      heading: '二度反射して届く光',
      detail: '基準視点で橙色の光源を点け、追跡回数を 1 から 2 に上げると、画像の中心に小さな橙色の像が現れます。最初の鏡で一次光線の空間を外れ、二枚目で光源へ向かう光路です。'
    },
    en: {
      label: 'Mirror relay', title: SCENE_TITLES.relay,
      heading: 'A reflection inside a reflection',
      detail: 'With the orange emitter on, raise the bounce limit from 1 to 2. A small orange image appears at the center. The first mirror sends this ray out of the primary-ray space, and the second directs it to the emitter.'
    }
  },
  cell24: {
    ja: {
      label: '正24胞体', title: '正24胞体の星座',
      heading: '24個の頂点を光でたどる',
      detail: '24個の頂点の一つが橙色の光源です。基準視点では、その像が別の頂点の鏡面に映ります。追跡回数を 0 にして ZW 回転を 0° から 45° に動かすと、光源を直接見られます。'
    },
    en: {
      label: '24-cell', title: SCENE_TITLES.cell24,
      heading: 'Follow the vertices of the 24-cell',
      detail: 'One of the 24 vertex spheres is the orange emitter. At the home view it appears in a mirrored vertex. To see the emitter directly, set the bounce limit to 0 and turn the ZW angle from 0° to 45°.'
    }
  },
  tir: {
    ja: {
      label: '全反射', title: 'ガラスの臨界角',
      heading: 'ガラスの内側で折り返す光',
      detail: '基準視点で中心の光線はガラスの内側から約45°で境界に当たり、全反射します。W 移動量を −0.30 にすると入射角が臨界角（約42°）を下回り、中央の橙色の反射像が弱くなります。'
    },
    en: {
      label: 'Critical angle', title: SCENE_TITLES.tir,
      heading: 'Total internal reflection',
      detail: 'At the home view, the center ray meets the inside of the glass at about 45° and reflects completely. Set the W offset to −0.30. The incidence angle falls below the critical angle (about 42°), dimming the orange reflection at the center.'
    }
  }
};

export function sceneCopy(language, scene) {
  return copy[scene][language];
}
