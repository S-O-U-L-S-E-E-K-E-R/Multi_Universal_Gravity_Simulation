// Two stacked sheets of matter on a periodic grid. Gravity acts within each sheet and across the gap between them.
// Particle-mesh solver: matter is spread onto a grid, the pull is found with FFTs, then read back at each particle.
// Space expands: positions are measured on a grid that stretches with space (scale factor a), so as space grows the
// pull between lumps weakens (by 1/a³) and the expansion brakes motions (2H). Lumps stay a spread-out web.
const TwoSheet = (() => {
  const N = 128, LOG = 7; // grid cells per side (box is N cells wide)

  function fft1(re, im, inv) { // in-place radix-2 FFT of length N
    for (let i = 1, j = 0; i < N; i++) {
      let bit = N >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= N; len <<= 1) {
      const ang = (inv ? 2 : -2) * Math.PI / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < N; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const a = i + k, b = a + len / 2, xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
          re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
          const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
  }
  const rowR = new Float64Array(N), rowI = new Float64Array(N);
  function fft2(re, im, inv) {
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) { rowR[x] = re[y * N + x]; rowI[x] = im[y * N + x]; }
      fft1(rowR, rowI, inv);
      for (let x = 0; x < N; x++) { re[y * N + x] = rowR[x]; im[y * N + x] = rowI[x]; }
    }
    for (let x = 0; x < N; x++) {
      for (let y = 0; y < N; y++) { rowR[y] = re[y * N + x]; rowI[y] = im[y * N + x]; }
      fft1(rowR, rowI, inv);
      for (let y = 0; y < N; y++) { re[y * N + x] = rowR[y]; im[y * N + x] = rowI[y]; }
    }
    if (inv) { const s = 1 / (N * N); for (let i = 0; i < N * N; i++) { re[i] *= s; im[i] *= s; } }
  }

  function rng(seed) { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
  const wrap = (v) => ((v % N) + N) % N;

  function create(opts) {
    const o = Object.assign({ perSide: 56, hiddenMass: 5, gap: 3, cross: true, coolVisible: 2.0, coolHidden: 0, G: 2000, tau: 3.8, kick: 0, kickShare: 0.05, soft: 1.0, seed: 1 }, opts);
    const n = o.perSide * o.perSide, sp = N / o.perSide;
    // independent random starting lumps on each sheet: alignment can only come from gravity across the gap
    function sheet(seed, mass) {
      const r = rng(seed), x = new Float64Array(n), y = new Float64Array(n), vx = new Float64Array(n), vy = new Float64Array(n);
      for (let j = 0, k = 0; j < o.perSide; j++) for (let i = 0; i < o.perSide; i++, k++) {
        x[k] = wrap((i + 0.5) * sp + (r() - 0.5) * sp * 0.9); y[k] = wrap((j + 0.5) * sp + (r() - 0.5) * sp * 0.9);
      }
      return { x, y, vx, vy, m: mass / n };
    }
    const s = { o, t: 0, vis: sheet(o.seed * 7919 + 1, 1), hid: sheet(o.seed * 104729 + 2, o.hiddenMass),
      rhoV: new Float64Array(N * N), rhoH: new Float64Array(N * N), kS: null, kC: null, history: [] };
    s.f = {}; for (const k of ["vSelfX", "vSelfY", "vTotX", "vTotY", "hTotX", "hTotY", "pxV", "pyV"]) s.f[k] = new Float64Array(N * N);
    buildKernels(s);
    return s;
  }
  function kernelFT(depth2, soft) {
    const re = new Float64Array(N * N), im = new Float64Array(N * N);
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const dx = x <= N / 2 ? x : x - N, dy = y <= N / 2 ? y : y - N;
      re[y * N + x] = -1 / Math.sqrt(dx * dx + dy * dy + depth2 + soft * soft);
    }
    fft2(re, im, false);
    return { re, im };
  }
  function buildKernels(s) { s.kS = kernelFT(0, s.o.soft); s.kC = kernelFT(s.o.gap * s.o.gap, s.o.soft); }

  function deposit(p, rho, fx, fy) { // cloud-in-cell; optionally also deposits momentum
    rho.fill(0); if (fx) { fx.fill(0); fy.fill(0); }
    for (let k = 0; k < p.x.length; k++) {
      const gx = p.x[k] - 0.5, gy = p.y[k] - 0.5, i0 = Math.floor(gx), j0 = Math.floor(gy), tx = gx - i0, ty = gy - j0;
      const i1 = wrap(i0 + 1), j1 = wrap(j0 + 1), ii = wrap(i0), jj = wrap(j0);
      const w = [[ii, jj, (1 - tx) * (1 - ty)], [i1, jj, tx * (1 - ty)], [ii, j1, (1 - tx) * ty], [i1, j1, tx * ty]];
      for (const [a, b, ww] of w) { const c = b * N + a; rho[c] += ww * p.m; if (fx) { fx[c] += ww * p.m * p.vx[k]; fy[c] += ww * p.m * p.vy[k]; } }
    }
  }
  function interp(field, x, y) {
    const gx = x - 0.5, gy = y - 0.5, i0 = Math.floor(gx), j0 = Math.floor(gy), tx = gx - i0, ty = gy - j0;
    const ii = wrap(i0), jj = wrap(j0), i1 = wrap(i0 + 1), j1 = wrap(j0 + 1);
    return field[jj * N + ii] * (1 - tx) * (1 - ty) + field[jj * N + i1] * tx * (1 - ty) + field[j1 * N + ii] * (1 - tx) * ty + field[j1 * N + i1] * tx * ty;
  }
  function gradient(phi, gx, gy, G) { // pull = -G * slope of the well
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      gx[y * N + x] = -G * (phi[y * N + wrap(x + 1)] - phi[y * N + wrap(x - 1)]) / 2;
      gy[y * N + x] = -G * (phi[wrap(y + 1) * N + x] - phi[wrap(y - 1) * N + x]) / 2;
    }
  }
  const A = { re: new Float64Array(N * N), im: new Float64Array(N * N) }, B = { re: new Float64Array(N * N), im: new Float64Array(N * N) };
  const T = { re: new Float64Array(N * N), im: new Float64Array(N * N) };
  function conv(out, a, k, b, kb) { // out = a*k (+ b*kb), back in real space
    for (let i = 0; i < N * N; i++) {
      let r = a.re[i] * k.re[i] - a.im[i] * k.im[i], m = a.re[i] * k.im[i] + a.im[i] * k.re[i];
      if (b) { r += b.re[i] * kb.re[i] - b.im[i] * kb.im[i]; m += b.re[i] * kb.im[i] + b.im[i] * kb.re[i]; }
      T.re[i] = r; T.im[i] = m;
    }
    T.re[0] = 0; T.im[0] = 0; // only differences from the average density pull
    fft2(T.re, T.im, true);
    out.set(T.re);
  }
  const phi = new Float64Array(N * N);
  function forces(s) {
    const { o, f } = s;
    deposit(s.vis, s.rhoV, f.pxV, f.pyV); deposit(s.hid, s.rhoH);
    A.re.set(s.rhoV); A.im.fill(0); fft2(A.re, A.im, false);
    B.re.set(s.rhoH); B.im.fill(0); fft2(B.re, B.im, false);
    // keep this simulation's own copy: A and B are shared work buffers that wellField() also uses
    if (!s.ftV) { s.ftV = { re: new Float64Array(N * N), im: new Float64Array(N * N) }; s.ftH = { re: new Float64Array(N * N), im: new Float64Array(N * N) }; }
    s.ftV.re.set(A.re); s.ftV.im.set(A.im); s.ftH.re.set(B.re); s.ftH.im.set(B.im);
    conv(phi, A, s.kS); gradient(phi, f.vSelfX, f.vSelfY, o.G);                          // ordinary matter's own pull
    conv(phi, A, s.kS, o.cross ? B : null, s.kC); gradient(phi, f.vTotX, f.vTotY, o.G);   // + the hidden sheet across the gap
    conv(phi, B, s.kS, o.cross ? A : null, s.kC); gradient(phi, f.hTotX, f.hTotY, o.G);
  }
  function move(s, p, gx, gy, dt, cool, px, py, rho) {
    const { o } = s, mean = 1 / (N * N) * (p.m * p.x.length);
    const a = scaleFactor(s), H = (2 / 3) / (o.tau + s.t), weak = 1 / (a * a * a), drag = 2 * H;
    for (let k = 0; k < p.x.length; k++) {
      let ax = weak * interp(gx, p.x[k], p.y[k]) - drag * p.vx[k], ay = weak * interp(gy, p.x[k], p.y[k]) - drag * p.vy[k];
      if (cool > 0 && px) { // gas cooling: random motions die away where matter is dense, so it sinks to the centre
        const r = interp(rho, p.x[k], p.y[k]);
        if (r > 1e-12) {
          const lvx = interp(px, p.x[k], p.y[k]) / r, lvy = interp(py, p.x[k], p.y[k]) / r, g = cool * Math.min(4, r / mean);
          ax -= g * (p.vx[k] - lvx); ay -= g * (p.vy[k] - lvy);
        }
      }
      p.vx[k] += ax * dt; p.vy[k] += ay * dt;
      p.x[k] = wrap(p.x[k] + p.vx[k] * dt); p.y[k] = wrap(p.y[k] + p.vy[k] * dt);
    }
  }
  const scaleFactor = (s) => Math.pow(1 + s.t / s.o.tau, 2 / 3); // how much space has stretched since the start
  const hpx = new Float64Array(N * N), hpy = new Float64Array(N * N);
  // Star explosions: every so often a fixed share of each galaxy's gas gets the same outward kick (a fixed energy per
  // unit mass). Gravity alone then decides whether it escapes: shallow wells lose it, deep wells pull it back.
  function blowOut(s) {
    const o = s.o, p = s.vis;
    const gals = galaxies(s, 30);
    for (const g of gals) {
      const [cx, cy] = centre(p, g.x, g.y, 3);
      for (let k = 0; k < p.x.length; k++) {
        const dx = md(p.x[k], cx), dy = md(p.y[k], cy), r = Math.hypot(dx, dy);
        if (r > 2.5 || Math.random() > o.kickShare) continue;
        const ux = r > 1e-3 ? dx / r : Math.cos(k), uy = r > 1e-3 ? dy / r : Math.sin(k);
        p.vx[k] += o.kick * ux; p.vy[k] += o.kick * uy;
      }
    }
  }
  function step(s, dt) {
    forces(s);
    s.nStep = (s.nStep || 0) + 1;
    if (s.o.kick > 0 && s.t > 8 && s.nStep % 10 === 0) blowOut(s);
    move(s, s.vis, s.f.vTotX, s.f.vTotY, dt, s.o.coolVisible, s.f.pxV, s.f.pyV, s.rhoV);
    if (s.o.coolHidden > 0) deposit(s.hid, s.rhoH, hpx, hpy);
    move(s, s.hid, s.f.hTotX, s.f.hTotY, dt, s.o.coolHidden, hpx, hpy, s.rhoH);
    s.t += dt;
  }

  // ---------- measurements ----------
  const smV = new Float64Array(N * N), smH = new Float64Array(N * N);
  function smooth(out, ft, sigma) {
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const kx = 2 * Math.PI * (x <= N / 2 ? x : x - N) / N, ky = 2 * Math.PI * (y <= N / 2 ? y : y - N) / N, w = Math.exp(-(kx * kx + ky * ky) * sigma * sigma / 2), i = y * N + x;
      T.re[i] = ft.re[i] * w; T.im[i] = ft.im[i] * w;
    }
    fft2(T.re, T.im, true); out.set(T.re);
  }
  function alignment(s) { // correlation between the two sheets' lumps (0 = unrelated, 1 = identical pattern)
    smooth(smV, s.ftV, 2); smooth(smH, s.ftH, 2);
    let ma = 0, mb = 0; for (let i = 0; i < N * N; i++) { ma += smV[i]; mb += smH[i]; } ma /= N * N; mb /= N * N;
    let sab = 0, saa = 0, sbb = 0;
    for (let i = 0; i < N * N; i++) { const a = smV[i] - ma, b = smH[i] - mb; sab += a * b; saa += a * a; sbb += b * b; }
    return sab / Math.sqrt(saa * sbb + 1e-30);
  }
  function galaxies(s, count) { // brightest clumps of ordinary matter
    smooth(smV, s.ftV, 1.5);
    const peaks = [];
    for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
      const v = smV[y * N + x]; let top = true;
      for (let dy = -3; dy <= 3 && top; dy++) for (let dx = -3; dx <= 3; dx++) if ((dx || dy) && smV[wrap(y + dy) * N + wrap(x + dx)] > v) { top = false; break; }
      if (top) peaks.push({ x: x + 0.5, y: y + 0.5, v });
    }
    return peaks.sort((a, b) => b.v - a.v).slice(0, count);
  }
  const md = (a, b) => { let d = a - b; d -= N * Math.round(d / N); return d; };
  function centre(p, gx, gy, R) { // centre of mass of p near (gx, gy), refined twice
    let cx = gx, cy = gy;
    for (let it = 0; it < 3; it++) {
      let sx = 0, sy = 0, c = 0;
      for (let k = 0; k < p.x.length; k++) { const dx = md(p.x[k], cx), dy = md(p.y[k], cy); if (dx * dx + dy * dy < R * R) { sx += dx; sy += dy; c++; } }
      if (c) { cx = wrap(cx + sx / c); cy = wrap(cy + sy / c); }
    }
    return [cx, cy];
  }
  function halfRadius(p, cx, cy, R) {
    const d = [];
    for (let k = 0; k < p.x.length; k++) { const dx = md(p.x[k], cx), dy = md(p.y[k], cy), r = Math.hypot(dx, dy); if (r < R) d.push(r); }
    d.sort((a, b) => a - b);
    return { half: d.length ? d[Math.floor(d.length / 2)] : 0, mass: d.length * p.m };
  }
  function rotation(s, cx, cy, rMax) { // orbital speed around a galaxy, from the actual pull on the sheet
    const out = { r: [], tot: [], self: [] };
    for (let r = 0.5; r <= rMax; r += 0.5) {
      let at = 0, as = 0;
      for (let k = 0; k < 48; k++) {
        const a = (k / 48) * 2 * Math.PI, ux = Math.cos(a), uy = Math.sin(a), x = wrap(cx + r * ux), y = wrap(cy + r * uy);
        at += -(interp(s.f.vTotX, x, y) * ux + interp(s.f.vTotY, x, y) * uy);
        as += -(interp(s.f.vSelfX, x, y) * ux + interp(s.f.vSelfY, x, y) * uy);
      }
      out.r.push(r); out.tot.push(Math.sqrt(Math.max(0, r * at / 48))); out.self.push(Math.sqrt(Math.max(0, r * as / 48)));
    }
    return out;
  }
  function measure(s, R = 12) {
    if (!s.ftV) forces(s);
    const al = alignment(s), gal = galaxies(s, 12);
    const stats = gal.map((g) => {
      const [cx, cy] = centre(s.vis, g.x, g.y, 4);
      const v = halfRadius(s.vis, cx, cy, R), h = halfRadius(s.hid, cx, cy, R);
      // offset between the galaxy and the centre of the hidden matter around it
      const [hx, hy] = centre(s.hid, cx, cy, 6);
      return { cx, cy, visHalf: v.half, hidHalf: h.half, visMass: v.mass, hidMass: h.mass, offset: Math.hypot(md(hx, cx), md(hy, cy)) };
    });
    return { t: s.t, alignment: al, galaxies: stats };
  }
  // Generic well map for other simulations: positions in grid units, returns the well (phi) and its pull (gx, gy), G = 1.
  const pmRho = new Float64Array(N * N), pmPhi = new Float64Array(N * N), pmGx = new Float64Array(N * N), pmGy = new Float64Array(N * N);
  const pmKernels = {};
  function wellField(xs, ys, ms, count, soft) {
    pmKernels[soft] = pmKernels[soft] || kernelFT(0, soft);
    pmRho.fill(0);
    for (let k = 0; k < count; k++) { // mass-weighted cloud-in-cell
      const gx = xs[k] - 0.5, gy = ys[k] - 0.5, i0 = Math.floor(gx), j0 = Math.floor(gy), tx = gx - i0, ty = gy - j0;
      const ii = wrap(i0), jj = wrap(j0), i1 = wrap(i0 + 1), j1 = wrap(j0 + 1), m = ms[k];
      pmRho[jj * N + ii] += m * (1 - tx) * (1 - ty); pmRho[jj * N + i1] += m * tx * (1 - ty); pmRho[j1 * N + ii] += m * (1 - tx) * ty; pmRho[j1 * N + i1] += m * tx * ty;
    }
    A.re.set(pmRho); A.im.fill(0); fft2(A.re, A.im, false);
    conv(pmPhi, A, pmKernels[soft]); gradient(pmPhi, pmGx, pmGy, 1);
    return { phi: pmPhi, gx: pmGx, gy: pmGy, rho: pmRho };
  }
  return { N, create, step, measure, rotation, buildKernels, galaxies, centre, wellField, interp, scaleFactor };
})();
if (typeof module !== "undefined") module.exports = TwoSheet;
