// B embeds the tracked 3D world isometrically into R4. Apply it to the
// initial exhibit as well, so the witness experiment has an exact solution.
const tilt = Math.PI / 10;
const cz = Math.cos(tilt), sz = Math.sin(tilt);
function baseline([x, y, z, w]) {
  return [x, y, cz*z - sz*w, sz*z + cz*w];
}
function ball(name, point, radius, material, color, emission = [0, 0, 0], witness = false) {
  return {name, center: baseline(point), radius, material, color, emission, witness};
}
const silver = [.89, .92, .99];
const transparent = [.96, .984, .997];
const black = [0, 0, 0];
function lamp(point, radius) {
  return ball('Off-axis emitter', point, radius, 'diffuse', black, [5.0, 1.16, .16], true);
}
const floor = ball('Reflecting floor', [0, -41.15, 0, 0], 40, 'mirror', [.32, .34, .42]);
const s = Math.SQRT1_2;

// The D4 roots are the 24 vertices of a regular 24-cell. Tilt and translate
// the entire crystal rigidly to arrange an optical path between two vertices.
function build24Cell() {
  const a = .95, radius = .46;
  const sin = s * radius / (2 * a);
  const cos = Math.sqrt(1 - sin * sin);
  const shiftW = a * (sin + cos) - s * radius;
  const roll = Math.PI / 8;
  const objects = [];
  for (let i = 0; i < 4; i++) {
    for (let j = i + 1; j < 4; j++) {
      for (const signI of [-1, 1]) {
        for (const signJ of [-1, 1]) {
          const v = [0, 0, 0, 0];
          v[i] = a * signI;
          v[j] = a * signJ;
          const [x, y, z, w] = v;
          const point = [x*Math.cos(roll) - y*Math.sin(roll),
            x*Math.sin(roll) + y*Math.cos(roll), cos*z + sin*w, cos*w - sin*z + shiftW];
          if (z > 0 && w > 0) {
            objects.push(lamp(point, radius));
          } else if (z > 0 && w < 0) {
            objects.push(ball('Axial mirror vertex', point, radius, 'mirror', silver));
          } else if (y > 0 && (w !== 0 || z < 0)) {
            objects.push(ball('Glass vertex', point, radius, 'glass', transparent));
          } else if (x > 0 && y > 0) {
            objects.push(ball('Reference vertex', point, radius, 'diffuse', black, [.18, .55, .64]));
          } else {
            const color = x < 0 ? [.32, .85, .84] : y < 0 ? [.79, .52, .93] : silver;
            objects.push(ball('Mirror vertex', point, radius, 'mirror', color));
          }
        }
      }
    }
  }
  return objects;
}

export const SCENES = {
  mirror: [
    ball('Witness mirror', [0, 0, -s, -s], 1, 'mirror', silver),
    lamp([0, 0, 0, 2], .58),
    ball('Turquoise mirror', [-1.11, -.48, -.03, .15], .51, 'mirror', [.21, .85, .81]),
    ball('Pearl mirror', [.90, -.54, .26, -.10], .47, 'mirror', [.82, .54, .96]),
    ball('Distant chrome', [1.41, .27, -.96, .09], .28, 'mirror', [.61, .76, .95]),
    ball('Distant glass', [-1.15, .53, -1.38, .12], .31, 'glass', transparent),
    floor,
  ],
  glass: [
    ball('Glass lens', [0, 0, -.25, -.62], 1.12, 'glass', transparent),
    // The center ray exits the second interface in direction
    // (0, 0, -.91418, -.40532), hitting this lamp two units later.
    lamp([0, 0, -3.18234, -1.24190], .51),
    ball('Jade mirror', [-1.18, -.57, -.96, -.14], .52, 'mirror', [.20, .86, .74]),
    ball('Violet mirror', [1.14, -.58, -.60, .11], .54, 'mirror', [.72, .46, .91]),
    ball('Blue marker', [.54, .35, -2.3, -.38], .32, 'diffuse', black, [.18, .49, .82]),
    floor,
  ],
  room: [
    ball('Left wall', [-42.10, 0, 0, 0], 40, 'mirror', [.70, .28, .36]),
    ball('Right wall', [42.10, 0, 0, 0], 40, 'mirror', [.22, .53, .70]),
    ball('Back wall', [0, 0, -43.15, 0], 40, 'diffuse', [.45, .47, .54], [.018, .023, .043]),
    ball('Ceiling', [0, 41.70, 0, 0], 40, 'diffuse', [.12, .13, .15], [.48, .56, .78]),
    floor,
    ball('W− wall', [0, 0, 0, -42.0], 40, 'mirror', [.48, .38, .60]),
    ball('W+ wall', [0, 0, 0, 42.4], 40, 'mirror', [.40, .52, .54]),
    ball('Chrome', [-.74, -.40, -.71, -.30], .81, 'mirror', silver),
    ball('Glass', [.78, -.53, -.09, .31], .72, 'glass', transparent),
    ball('Matte', [.16, -.71, -1.87, .11], .43, 'diffuse', [.55, .72, .67]),
    lamp([-.64, .59, -.43, 1.30], .66),
  ],
  relay: [
    // Hit the first port at the origin and the second at W=1.30. Their
    // 45-degree normals turn −Z into +W, then +W into +X.
    ball('Entrance mirror', [0, 0, -1.22*s, -1.22*s], 1.22, 'mirror', silver),
    ball('Turning mirror', [-s, 0, 0, 1.30+s], 1, 'mirror', [.75, .91, .94]),
    lamp([1.15, 0, 0, 1.30], .42),
    ...Array.from({length: 6}, (_, i) => {
      const a = (i + .5) * Math.PI / 3;
      const p = [1.15*Math.cos(a), 1.05*Math.sin(a), -.54, -.04];
      return ball('Entrance fiducial', p, .13, 'mirror', [.48, .74, .83], [.02, .07, .09]);
    }),
    ball('Upper exit fiducial', [1.15, .57, 0, 1.30], .12, 'mirror', silver, [.08, .18, .20]),
    ball('Lower exit fiducial', [1.15, -.57, 0, 1.30], .12, 'mirror', silver, [.08, .18, .20]),
    floor,
  ],
  cell24: build24Cell(),
  tir: [
    // At home the eye is inside the glass. The parallel center ray's impact
    // parameter is 2.4, so sin(theta)=2.4/3.4 > 1/1.5 at the exit boundary.
    ball('Interior glass boundary', [0, 0, 3.10, -2.40], 3.4, 'glass', transparent,
      [.002, .003, .007]),
    lamp([0, 0, .69652537, -1.39999162], .43),
    ball('Turquoise beacon', [-.89, .31, .72, -1.42], .32, 'mirror',
      [.75, .91, .94], [.085, .41, .47]),
    ball('Violet beacon', [.84, -.28, .61, -1.35], .36, 'mirror',
      [.85, .71, .94], [.29, .14, .48]),
    ball('Blue beacon', [.58, .70, .68, -1.55], .21, 'mirror',
      [.72, .79, .94], [.12, .21, .40]),
    ball('Lower reference', [-.39, -.73, .82, -1.21], .17, 'mirror',
      silver, [.26, .32, .37]),
    ball('Exterior mirror', [-.80, -.17, -1.75, .15], .51, 'mirror', silver),
    ball('Exterior glass', [.78, .26, -1.82, .35], .37, 'glass', transparent),
  ],
};
export const SCENE_TITLES = {
  mirror: 'Off-Axis Mirror', glass: 'Glass Lens', room: 'Four-Dimensional Room',
  relay: 'Two-Mirror Relay', cell24: '24-Cell Constellation', tir: 'Critical Angle',
};

// Column-major isometry: camera-space origin and film axes are transformed
// together under the user's 4D orbit. W travel is in fixed world coordinates.
export function embeddingMatrix(state) {
  const axis = 'xyz'.indexOf(state.plane[0]);
  const c = Math.cos(state.angle), s = Math.sin(state.angle);
  const result = new Float32Array(16);
  for (let i = 0; i < 4; i++) {
    const p = baseline([Number(i===0), Number(i===1), Number(i===2), Number(i===3)]);
    const a = p[axis], w = p[3];
    p[axis] = c*a - s*w;
    p[3] = s*a + c*w;
    result.set(p, 4*i);
  }
  return result;
}

export function cameraFrame(transform, state) {
  const basis = embeddingMatrix(state);
  function column(j) {
    return Array.from({length: 4}, (_, k) => basis[k]*transform[j] +
      basis[4+k]*transform[j+1] + basis[8+k]*transform[j+2]);
  }
  const origin = column(12);
  origin[3] += state.w;
  return {origin, right: column(0), up: column(4), back: column(8)};
}
