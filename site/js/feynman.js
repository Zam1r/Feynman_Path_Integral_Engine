/**
 * Browser port of the Euclidean Brownian-bridge reweighting used in
 * cpp/feynman/direct.hpp (bridge_potential_weight / sample_bridges /
 * a light direct_ground_density for interactive |φ(1)|² estimates).
 */
(function (global) {
  "use strict";

  function mulberry32(seed) {
    let t = seed >>> 0;
    return function () {
      t = (t + 0x6d2b79f5) >>> 0;
      let r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
  }

  function makeNormal(rand) {
    let spare = null;
    return function () {
      if (spare !== null) {
        const v = spare;
        spare = null;
        return v;
      }
      let u, v, s;
      do {
        u = rand() * 2 - 1;
        v = rand() * 2 - 1;
        s = u * u + v * v;
      } while (s >= 1 || s === 0);
      const m = Math.sqrt((-2 * Math.log(s)) / s);
      spare = v * m;
      return u * m;
    };
  }

  function potentialHO(x) {
    return 0.5 * x * x;
  }

  function potentialAnharm(x) {
    const x2 = x * x;
    return 0.5 * x2 + 0.1 * x2 * x2;
  }

  function potentialDoubleWell(x) {
    const u = x * x - 2.25;
    return 0.8 * u * u;
  }

  const POTENTIALS = {
    ho: potentialHO,
    anharm: potentialAnharm,
    doublewell: potentialDoubleWell,
  };

  /**
   * Faithful port of bridge_potential_weight from direct.hpp.
   * Samples a free Brownian bridge pinned at xb and returns the
   * potential reweight exp(-Σ ε V(mid)/ħ).
   */
  function bridgePotentialWeight(xb, N, eps, mass, hbar, V, normal, pathOut) {
    const sigma2 = (eps * hbar) / mass;
    let xPrev = xb;
    let spot = 0;
    if (pathOut) {
      pathOut.length = 0;
      pathOut.push(xb);
    }
    for (let j = 1; j < N; j++) {
      const steps = N - j + 1;
      const mean = xPrev + (xb - xPrev) / steps;
      const variance = sigma2 * ((steps - 1) / steps);
      const x = mean + Math.sqrt(Math.max(variance, 0)) * normal();
      const mid = 0.5 * (x + xPrev);
      spot += (eps / hbar) * V(mid);
      xPrev = x;
      if (pathOut) pathOut.push(x);
    }
    const mid = 0.5 * (xb + xPrev);
    spot += (eps / hbar) * V(mid);
    if (pathOut) pathOut.push(xb);
    if (spot > 700) return 0;
    return Math.exp(-spot);
  }

  function sampleBridges(xb, N, T, mass, hbar, potName, npaths, seed) {
    const V = POTENTIALS[potName] || potentialHO;
    const eps = T / N;
    const rand = mulberry32(seed >>> 0);
    const normal = makeNormal(rand);
    const out = [];
    for (let i = 0; i < npaths; i++) {
      const path = [];
      const weight = bridgePotentialWeight(xb, N, eps, mass, hbar, V, normal, path);
      out.push({ x: path.slice(), weight });
    }
    return out;
  }

  function interpLinear(xs, ys, xq) {
    const n = xs.length;
    if (!n) return 0;
    if (xq <= xs[0]) return ys[0];
    if (xq >= xs[n - 1]) return ys[n - 1];
    let lo = 0;
    let hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (xs[mid] <= xq) lo = mid;
      else hi = mid;
    }
    const t = (xq - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] * (1 - t) + ys[hi] * t;
  }

  /**
   * Interactive direct MC density estimate (same algorithm as
   * direct_ground_density, fewer samples for browser latency).
   */
  function directGroundDensity(cfg) {
    const {
      N,
      T,
      mass = 1,
      hbar = 1,
      potName = "ho",
      xmin = -4,
      xmax = 4,
      dx = 0.2,
      samples = 400,
      batches = 4,
      seed = 42,
      targetX = 1,
    } = cfg;
    const V = POTENTIALS[potName] || potentialHO;
    const eps = T / N;
    const nx = Math.floor((xmax - xmin) / dx + 0.5) + 1;
    const xs = new Array(nx);
    for (let i = 0; i < nx; i++) xs[i] = xmin + i * dx;
    const per = Math.max(1, Math.floor(samples / batches));
    const batchG = Array.from({ length: batches }, () => new Array(nx).fill(0));

    for (let i = 0; i < nx; i++) {
      const xb = xs[i];
      const rand = mulberry32((seed + Math.imul(i + 1, 0x9e3779b9)) >>> 0);
      const normal = makeNormal(rand);
      for (let b = 0; b < batches; b++) {
        let sum = 0;
        for (let s = 0; s < per; s++) {
          sum += bridgePotentialWeight(xb, N, eps, mass, hbar, V, normal, null);
        }
        batchG[b][i] = sum / per;
      }
    }

    const psi2 = new Array(nx).fill(0);
    const err = new Array(nx).fill(0);
    for (let b = 0; b < batches; b++) {
      let Z = 0.5 * batchG[b][0] + 0.5 * batchG[b][nx - 1];
      for (let i = 1; i < nx - 1; i++) Z += batchG[b][i];
      Z *= dx;
      if (!(Z > 0)) continue;
      for (let i = 0; i < nx; i++) psi2[i] += batchG[b][i] / Z;
    }
    for (let i = 0; i < nx; i++) psi2[i] /= batches;

    for (let b = 0; b < batches; b++) {
      let Z = 0.5 * batchG[b][0] + 0.5 * batchG[b][nx - 1];
      for (let i = 1; i < nx - 1; i++) Z += batchG[b][i];
      Z *= dx;
      if (!(Z > 0)) continue;
      for (let i = 0; i < nx; i++) {
        const rho = batchG[b][i] / Z;
        const d = rho - psi2[i];
        err[i] += d * d;
      }
    }
    const scale = 1 / Math.sqrt(batches * (batches - 1));
    for (let i = 0; i < nx; i++) err[i] = Math.sqrt(err[i]) * scale;

    const analytic = (x) => Math.exp(-((x * x) / hbar)) / Math.sqrt(Math.PI * hbar);

    return {
      x: xs,
      psi2,
      err,
      psi2_at_target: interpLinear(xs, psi2, targetX),
      err_at_target: interpLinear(xs, err, targetX),
      analytic_at_target: analytic(targetX),
      analytic,
    };
  }

  function analyticPsi2(x, hbar) {
    return Math.exp(-((x * x) / hbar)) / Math.sqrt(Math.PI * hbar);
  }

  global.Feynman = {
    POTENTIALS,
    bridgePotentialWeight,
    sampleBridges,
    directGroundDensity,
    analyticPsi2,
    interpLinear,
    mulberry32,
    makeNormal,
  };
})(window);
