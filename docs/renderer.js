import { createLegend } from './xr-legend.js';
import { SCENES, embeddingMatrix } from './optics.js';
export { SCENE_TITLES } from './optics.js';

const vertexSource = `
attribute vec2 aPosition;
varying vec2 vNdc;
void main() {
  vNdc = aPosition;
  gl_Position = vec4(aPosition, 0., 1.);
}`;

function sceneSource(objects) {
  const literal = n => n.toFixed(9);
  const vector = v => `vec${v.length}(${v.map(literal).join(',')})`;
  return objects.map(s => {
    const material = {diffuse: 0, mirror: 1, glass: 2}[s.material];
    const emission = (s.witness ? 'float(uLightOn)*' : '') + vector(s.emission);
    return `testBall(h, ro, rd, ${vector(s.center)}, ${literal(s.radius)}, ${material}, ${vector(s.color)}, ${emission});`;
  }).join('\n');
}

function fragmentSource(objects) { return `
precision highp float;
varying vec2 vNdc;
uniform mat4 uCamera, uBasis;
uniform vec4 uLens;
uniform vec2 uPixelNdc, uJitter;
uniform float uTravel, uSample, uPixel;
uniform int uBounces, uMode;
uniform bool uLightOn, uImmersive;
uniform sampler2D uLegend, uHistory;
const float FAR = 10000.;
const float EPS = .00015;
const float TAU = 6.28318530718;

struct Hit { float t; vec4 center; vec3 albedo; vec3 emission; int material; };
struct Path { vec4 ro; vec4 rd; vec3 weight; int depth; };

void testBall(inout Hit h, vec4 ro, vec4 rd, vec4 center, float radius,
              int material, vec3 albedo, vec3 emission) {
  vec4 delta = center - ro;
  float b = dot(delta, rd);
  vec4 perp = delta - b*rd;
  float disc = radius*radius - dot(perp,perp);
  if (disc <= 0.) return;
  float root = sqrt(disc);
  float t = b - root;
  if (t <= EPS) t = b + root;
  if (t > EPS && t < h.t) h = Hit(t,center,albedo,emission,material);
}

Hit nearest(vec4 ro, vec4 rd) {
  Hit h = Hit(FAR, vec4(0.), vec3(0.), vec3(0.), 0);
  ${sceneSource(objects)}
  return h;
}

float hash(vec2 p) {
  vec3 v = fract(vec3(p.xyx) * vec3(.1031,.11369,.13787));
  v += dot(v, v.yzx + 19.19);
  return fract((v.x+v.y)*v.z);
}
float random(inout vec2 rng) {
  rng += vec2(13.147, 7.193);
  return hash(rng);
}

// Analytic radiance field on the direction three-sphere. Large angular
// windows form studio softboxes, visible whenever a ray escapes toward one.
vec3 environment(vec4 worldDir) {
  const float c = .9510565163, s = .3090169944;
  vec4 d = vec4(worldDir.xy, c*worldDir.z+s*worldDir.w, -s*worldDir.z+c*worldDir.w);
  float canopy = smoothstep(-.69,.28,d.z+.91*abs(d.w));
  vec3 radiance = vec3(.0028,.0037,.0082)
    + vec3(.0014,.002,.0032) * smoothstep(-.45,.56,d.y);
  radiance += canopy * (vec3(.044,.061,.118) + vec3(.08,.11,.15)*(.5+.5*d.y));
  float softbox = smoothstep(.15,.48,-d.x+.60*d.y+.21*d.w)
    * smoothstep(-.20,.20,d.z+.82*d.w+.30*d.y)
    * (1.-smoothstep(.75,.95,d.y));
  radiance += vec3(1.80,2.13,2.78) * softbox;
  float strip = 1.-smoothstep(.035,.21,abs(d.y+.22*d.x-.40));
  radiance += canopy*strip*vec3(.81,.91,1.30);
  float violet = smoothstep(.16,.67,d.x-.32*d.y) * smoothstep(-.15,.55,d.w+.64*d.z);
  radiance += violet * vec3(.67,.28,.89);
  float teal = smoothstep(.31,.78,-d.y-.24*d.x) * smoothstep(-.28,.44,d.z-.74*d.w);
  radiance += teal * vec3(.07,.56,.69);
  if (d.z < -.62) {
    vec2 p = (d.xy+.35*d.w)/(-d.z)*43.;
    vec2 cell = floor(p);
    vec2 delta = fract(p)-.5-(vec2(hash(cell+3.),hash(cell+6.))-.5)*.7;
    float star = (1.-smoothstep(.012,.050,length(delta))) * smoothstep(.63,1.,hash(cell));
    radiance += star * vec3(.15,.18,.26);
  }
  return radiance;
}

// Uniform on S³ via two independent angles and a uniform squared radius.
// Folding across the tangent hyperplane reproduces smallpt4D's diffuse law.
vec4 hemisphere(vec4 normal, inout vec2 rng) {
  float u = random(rng), a = TAU*random(rng), b = TAU*random(rng);
  vec4 v = vec4(sqrt(1.-u)*vec2(cos(a),sin(a)), sqrt(u)*vec2(cos(b),sin(b)));
  return dot(v,normal) < 0. ? -v : v;
}

vec3 radiance(vec4 ro, vec4 rd, inout vec2 rng, out float primaryDistance) {
  Path p = Path(ro,rd,vec3(1.),0);
  Path saved0 = p, saved1 = p;
  int pending = 0;
  vec3 sum = vec3(0.);
  primaryDistance = FAR;
  // Split dielectrics at the first two path vertices, as in the native
  // renderers. At most two branches await evaluation and four leaves exist.
  for (int branch=0; branch<4; ++branch) {
    if (branch>0) {
      if (pending==0) break;
      p = saved0;
      if (pending==2) p = saved1;
      pending--;
    }
    for (int step=0; step<7; ++step) {
      Hit h = nearest(p.ro,p.rd);
      if (p.depth==0) primaryDistance = h.t;
      if (h.t==FAR) { sum += p.weight*environment(p.rd); break; }
      sum += p.weight*h.emission;
      if (p.depth>=uBounces || max(h.albedo.r,max(h.albedo.g,h.albedo.b))==0.) break;
      p.ro += h.t*p.rd;
      vec4 n = normalize(p.ro-h.center);
      bool entering = dot(n,p.rd)<0.;
      vec4 nl = entering ? n : -n;
      p.weight *= h.albedo;
      p.depth++;
      if (h.material==0) {
        p.rd = hemisphere(nl,rng);
      } else if (h.material==1) {
        p.rd = normalize(reflect(p.rd,n));
      } else {
        float eta = entering ? 1./1.5 : 1.5;
        float cosIn = clamp(-dot(p.rd,nl),0.,1.);
        float cosOut2 = 1.-eta*eta*(1.-cosIn*cosIn);
        vec4 reflected = normalize(reflect(p.rd,n));
        if (cosOut2<=0.) {
          p.rd = reflected;
        } else {
          vec4 transmitted = normalize(eta*p.rd+(eta*cosIn-sqrt(cosOut2))*nl);
          float cosF = entering ? cosIn : clamp(dot(transmitted,n),0.,1.);
          float fresnel = .04+.96*pow(1.-cosF,5.);
          if (p.depth<=2) {
            Path next = Path(p.ro,reflected,p.weight*fresnel,p.depth);
            if (pending==0) saved0=next; else saved1=next;
            pending++;
            p.weight *= 1.-fresnel;
            p.rd = transmitted;
          } else {
            // Select by the Fresnel weights, so the Monte Carlo probability
            // cancels its corresponding reflection/transmission coefficient.
            p.rd = random(rng)<fresnel ? reflected : transmitted;
          }
        }
      }
    }
  }
  return sum;
}

vec3 displayColor(vec3 light) {
  return pow(light/(1.+light), vec3(1./2.2));
}
// RGBE preserves HDR radiance on WebGL 1 render targets without extensions.
vec4 encodeHDR(vec3 light) {
  float e = ceil(log2(max(.0000001,max(light.r,max(light.g,light.b)))));
  return vec4(light*exp2(-e),(e+64.)/255.);
}
vec3 decodeHDR(vec4 rgba) { return rgba.rgb*exp2(floor(rgba.a*255.+.5)-64.); }

void legend(inout vec3 col, vec3 ro, vec3 rd, float surface) {
  const vec3 center = vec3(0.,-.80,1.50);
  const vec3 normal = vec3(0.,.38461538,.92307692);
  const vec3 up = vec3(0.,.92307692,-.38461538);
  const vec2 halfSize = vec2(.72,.21);
  float facing = dot(rd,normal);
  if (abs(facing)<.001) return;
  float t = dot(center-ro,normal)/facing;
  if (t<=.001 || t>=surface) return;
  vec3 hit = ro+t*rd-center;
  vec2 xy = vec2(hit.x,dot(hit,up));
  vec2 q = abs(xy)-halfSize+.046;
  float edge = min(max(q.x,q.y),0.)+length(max(q,0.))-.046;
  if (edge>0.) return;
  vec2 uv = .5+xy/(2.*halfSize)*vec2(-sign(facing),-1.);
  vec3 label = texture2D(uLegend,uv).rgb;
  float px = uPixel*t;
  label = mix(label,vec3(.42,.41,.55),smoothstep(-.012-px,-.006,edge));
  col = mix(col,label,1.-smoothstep(-px,0.,edge));
}

void main() {
  vec2 uv = vNdc*.5+.5;
  if (uMode==2) {
    gl_FragColor = vec4(displayColor(decodeHDR(texture2D(uHistory,uv))),1.);
    return;
  }
  vec3 eyeDir = normalize(vec3((vNdc+uJitter*uPixelNdc+uLens.zw)/uLens.xy,-1.));
  vec3 ray3 = normalize(mat3(uCamera)*eyeDir);
  vec3 eye3 = uCamera[3].xyz;
  vec4 ro = uBasis*vec4(eye3,0.) + vec4(0.,0.,0.,uTravel);
  vec4 rd = normalize(uBasis*vec4(ray3,0.));
  vec2 rng = gl_FragCoord.xy + (uSample+1.)*vec2(31.17,67.37);
  float surface;
  vec3 light = radiance(ro,rd,rng,surface);
  if (uMode==1) {
    if (uSample>0.) light = mix(decodeHDR(texture2D(uHistory,uv)),light,1./(uSample+1.));
    gl_FragColor = encodeHDR(light);
  } else {
    vec3 col = displayColor(light);
    if (uImmersive) legend(col,eye3,ray3,surface);
    gl_FragColor = vec4(col,1.);
  }
}`; }

function compile(gl, type, source) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader,source);
  gl.compileShader(shader);
  if (gl.getShaderParameter(shader,gl.COMPILE_STATUS)) return shader;
  const message = gl.getShaderInfoLog(shader);
  gl.deleteShader(shader);
  throw Error(message || 'Shader compilation failed');
}

function makeProgram(gl, objects) {
  const vert = compile(gl,gl.VERTEX_SHADER,vertexSource);
  const frag = compile(gl,gl.FRAGMENT_SHADER,fragmentSource(objects));
  const program = gl.createProgram();
  gl.attachShader(program,vert);
  gl.attachShader(program,frag);
  gl.linkProgram(program);
  gl.deleteShader(vert);
  gl.deleteShader(frag);
  if (!gl.getProgramParameter(program,gl.LINK_STATUS)) {
    const message=gl.getProgramInfoLog(program);
    gl.deleteProgram(program);
    throw Error(message || 'WebGL program could not be linked');
  }
  const u = {};
  for (const name of ['Camera','Basis','Lens','PixelNdc','Jitter','Travel','Sample',
    'Pixel','Bounces','Mode','LightOn','Immersive','Legend','History']) {
    u[name] = gl.getUniformLocation(program,'u'+name);
  }
  gl.useProgram(program);
  gl.uniform1i(u.Legend,0);
  gl.uniform1i(u.History,1);
  return {program,u,position:gl.getAttribLocation(program,'aPosition')};
}

function halton(i, base) {
  let result=0, scale=1;
  while (i>0) {
    scale/=base;
    result+=scale*(i%base);
    i=Math.floor(i/base);
  }
  return result;
}

export function createRenderer(canvas) {
  const gl = canvas.getContext('webgl',{alpha:false,antialias:false,depth:false,stencil:false});
  if (!gl) throw Error('WebGL is unavailable');
  let programs, triangle, label, history, framebuffer, historySize, previous, count;

  function restore() {
    programs = new Map();
    programs.set('mirror',makeProgram(gl,SCENES.mirror));
    triangle=gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER,triangle);
    gl.bufferData(gl.ARRAY_BUFFER,new Float32Array([-1,-1,3,-1,-1,3]),gl.STATIC_DRAW);
    label=createLegend(gl);
    framebuffer=gl.createFramebuffer();
    gl.activeTexture(gl.TEXTURE1);
    history=[gl.createTexture(),gl.createTexture()];
    for (const tex of history) {
      gl.bindTexture(gl.TEXTURE_2D,tex);
      gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MIN_FILTER,gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_MAG_FILTER,gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_S,gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D,gl.TEXTURE_WRAP_T,gl.CLAMP_TO_EDGE);
      gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,1,1,0,gl.RGBA,gl.UNSIGNED_BYTE,null);
    }
    historySize=[0,0];
    previous=null;
    count=0;
    // The alpha byte stores an exact radiance exponent in the history buffer.
    gl.disable(gl.DITHER);
  }

  function draw({viewport,projection,transform,state,immersive=false,frameTime,progressive=false}) {
    const {x,y,width,height}=viewport;
    if (!programs.has(state.scene)) programs.set(state.scene,makeProgram(gl,SCENES[state.scene]));
    const {program,u,position}=programs.get(state.scene);
    gl.useProgram(program);
    gl.bindBuffer(gl.ARRAY_BUFFER,triangle);
    gl.enableVertexAttribArray(position);
    gl.vertexAttribPointer(position,2,gl.FLOAT,false,0,0);
    gl.uniformMatrix4fv(u.Camera,false,transform);
    gl.uniformMatrix4fv(u.Basis,false,embeddingMatrix(state));
    gl.uniform4f(u.Lens,projection[0],projection[5],projection[8],projection[9]);
    gl.uniform2f(u.PixelNdc,2/width,2/height);
    gl.uniform1f(u.Pixel,2/(Math.abs(projection[5])*height));
    gl.uniform1f(u.Travel,state.w);
    gl.uniform1i(u.Bounces,state.bounces);
    gl.uniform1i(u.LightOn,state.lightOn ? 1 : 0);
    gl.uniform1i(u.Immersive,immersive ? 1 : 0);
    if (immersive) label.update(state,frameTime ?? performance.now());

    if (!progressive || immersive) {
      gl.viewport(x,y,width,height);
      gl.uniform2f(u.Jitter,0,0);
      gl.uniform1f(u.Sample,0);
      gl.uniform1i(u.Mode,0);
      gl.drawArrays(gl.TRIANGLES,0,3);
      return false;
    }

    gl.activeTexture(gl.TEXTURE1);
    if (historySize[0]!==width || historySize[1]!==height) {
      for (const tex of history) {
        gl.bindTexture(gl.TEXTURE_2D,tex);
        gl.texImage2D(gl.TEXTURE_2D,0,gl.RGBA,width,height,0,gl.RGBA,gl.UNSIGNED_BYTE,null);
      }
      historySize=[width,height];
      previous=null;
    }
    const key=[state.scene,state.w,state.plane,state.angle,state.lightOn,state.bounces,
      ...projection,...transform].join(',');
    if (key!==previous) { previous=key; count=0; }

    if (count<48) {
      const target=gl.getParameter(gl.FRAMEBUFFER_BINDING);
      gl.bindFramebuffer(gl.FRAMEBUFFER,framebuffer);
      gl.framebufferTexture2D(gl.FRAMEBUFFER,gl.COLOR_ATTACHMENT0,gl.TEXTURE_2D,history[(count+1)%2],0);
      gl.bindTexture(gl.TEXTURE_2D,history[count%2]);
      gl.viewport(0,0,width,height);
      gl.uniform1i(u.Mode,1);
      gl.uniform1f(u.Sample,count);
      gl.uniform2f(u.Jitter,count ? halton(count,2)-.5 : 0,count ? halton(count,3)-.5 : 0);
      gl.drawArrays(gl.TRIANGLES,0,3);
      gl.bindFramebuffer(gl.FRAMEBUFFER,target);
      count++;
    }
    gl.bindTexture(gl.TEXTURE_2D,history[count%2]);
    gl.viewport(x,y,width,height);
    gl.uniform1i(u.Mode,2);
    gl.drawArrays(gl.TRIANGLES,0,3);
    return count<48;
  }

  function dispose() {
    label.dispose();
    for (const tex of history) gl.deleteTexture(tex);
    gl.deleteFramebuffer(framebuffer);
    gl.deleteBuffer(triangle);
    for (const {program} of programs.values()) gl.deleteProgram(program);
  }
  restore();
  return {gl,draw,restore,dispose};
}
