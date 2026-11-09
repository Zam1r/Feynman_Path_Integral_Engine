/**
 * Interactive lab: bridge sampling, live density estimate, N–T lattice lookup.
 */
(function () {
  "use strict";

  const F = window.Feynman;
  if (!F) return;

  const state = {
    N: 7,
    T: 4,
    hbar: 1,
    mass: 1,
    pot: "ho",
    xb: 1,
    nBridges: 24,
    seed: 99,
    tab: "bridges",
    density: null,
    bridges: null,
    exactAlpha: [],
    suiteSummary: {},
  };

  const els = {
    N: document.getElementById("ctrl-N"),
    T: document.getElementById("ctrl-T"),
    hbar: document.getElementById("ctrl-hbar"),
    pot: document.getElementById("ctrl-pot"),
    outN: document.getElementById("out-N"),
    outT: document.getElementById("out-T"),
    outH: document.getElementById("out-hbar"),
    outEps: document.getElementById("out-eps"),
    outMC: document.getElementById("out-mc"),
    outExact: document.getElementById("out-exact"),
    outAnalytic: document.getElementById("out-analytic"),
    status: document.getElementById("lab-status"),
    bridgeCanvas: document.getElementById("canvas-bridges"),
    densityCanvas: document.getElementById("canvas-density"),
    ntCanvas: document.getElementById("nt-canvas"),
    ntReadout: document.getElementById("nt-readout"),
    caption: document.getElementById("viz-caption"),
  };

  function fmt(x, d = 4) {
    if (x == null || Number.isNaN(x)) return "—";
    return Number(x).toFixed(d);
  }

  function setStatus(msg) {
    if (els.status) els.status.textContent = msg;
  }

  function syncOutputs() {
    const eps = state.T / state.N;
    if (els.outN) els.outN.textContent = String(state.N);
    if (els.outT) els.outT.textContent = fmt(state.T, 2);
    if (els.outH) els.outH.textContent = fmt(state.hbar, 2);
    if (els.outEps) els.outEps.textContent = fmt(eps, 3);
  }

  function sampleBridges() {
    state.bridges = F.sampleBridges(
      state.xb,
      state.N,
      state.T,
      state.mass,
      state.hbar,
      state.pot,
      state.nBridges,
      state.seed
    );
  }

  function estimateDensity() {
    // Browser-friendly sample count; suite uses far more.
    const samples = state.N <= 12 ? 480 : 320;
    state.density = F.directGroundDensity({
      N: state.N,
      T: state.T,
      mass: state.mass,
      hbar: state.hbar,
      potName: state.pot,
      xmin: -4,
      xmax: 4,
      dx: 0.25,
      samples,
      batches: 4,
      seed: 42 + state.N * 17 + Math.round(state.T * 10),
      targetX: 1,
    });
    if (els.outMC) {
      els.outMC.textContent =
        fmt(state.density.psi2_at_target, 4) + " ± " + fmt(state.density.err_at_target, 4);
    }
    if (els.outAnalytic) {
      els.outAnalytic.textContent = fmt(state.density.analytic_at_target, 5);
    }
    // Exact lattice for HO only, from suite table when N,T match; else blank.
    updateExactFromSuite();
  }

  function updateExactFromSuite() {
    if (!els.outExact) return;
    const hit = state.exactAlpha.find(
      (r) => Number(r.N) === state.N && Math.abs(Number(r.T) - state.T) < 1e-9
    );
    if (hit && state.pot === "ho" && Math.abs(state.hbar - 1) < 1e-9) {
      els.outExact.textContent = fmt(hit.psi2_1, 5);
    } else if (state.suiteSummary.lattice_psi2_1 && state.N === 7 && state.T === 4 && state.pot === "ho") {
      els.outExact.textContent = fmt(state.suiteSummary.lattice_psi2_1, 5);
    } else {
      els.outExact.textContent = state.pot === "ho" ? "n/a (off-grid)" : "n/a";
    }
  }

  function drawBridges() {
    const canvas = els.bridgeCanvas;
    if (!canvas || !state.bridges) return;
    const ctx = canvas.getContext("2d");
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const cssW = canvas.clientWidth || 640;
    const cssH = canvas.clientHeight || 360;
    canvas.width = Math.floor(cssW * dpr);
    canvas.height = Math.floor(cssH * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const pad = { l: 48, r: 18, t: 22, b: 40 };
    const W = cssW - pad.l - pad.r;
    const H = cssH - pad.t - pad.b;

    ctx.clearRect(0, 0, cssW, cssH);
    ctx.fillStyle = "#16120e";
    ctx.fillRect(0, 0, cssW, cssH);

    // axes range
    let xmin = Infinity;
    let xmax = -Infinity;
    state.bridges.forEach((p) => {
      p.x.forEach((v) => {
        if (v < xmin) xmin = v;
        if (v > xmax) xmax = v;
      });
    });
    const padX = Math.max(0.4, (xmax - xmin) * 0.08);
    xmin -= padX;
    xmax += padX;
    const tmax = state.T;

    const xMap = (x) => pad.l + ((x - xmin) / (xmax - xmin)) * W;
    const tMap = (tau) => pad.t + (tau / tmax) * H;

    // grid
    ctx.strokeStyle = "rgba(196,163,90,0.12)";
    ctx.lineWidth = 1;
    for (let i = 0; i <= state.N; i++) {
      const tau = (i / state.N) * tmax;
      const y = tMap(tau);
      ctx.beginPath();
      ctx.moveTo(pad.l, y);
      ctx.lineTo(pad.l + W, y);
      ctx.stroke();
    }

    // classical straight line xb
    ctx.strokeStyle = "rgba(196,163,90,0.35)";
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(xMap(state.xb), tMap(0));
    ctx.lineTo(xMap(state.xb), tMap(tmax));
    ctx.stroke();
    ctx.setLineDash([]);

    const maxW = Math.max(...state.bridges.map((p) => p.weight), 1e-12);
    state.bridges.forEach((p) => {
      const alpha = 0.12 + 0.55 * (p.weight / maxW);
      ctx.strokeStyle = `rgba(168, 192, 176, ${alpha})`;
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      p.x.forEach((xv, i) => {
        const tau = (i / (p.x.length - 1)) * tmax;
        const X = xMap(xv);
        const Y = tMap(tau);
        if (i === 0) ctx.moveTo(X, Y);
        else ctx.lineTo(X, Y);
      });
      ctx.stroke();
      // beads
      p.x.forEach((xv, i) => {
        const tau = (i / (p.x.length - 1)) * tmax;
        ctx.beginPath();
        ctx.arc(xMap(xv), tMap(tau), 2.1, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(239, 231, 216, ${0.15 + 0.45 * (p.weight / maxW)})`;
        ctx.fill();
      });
    });

    // labels
    ctx.fillStyle = "#8f8678";
    ctx.font = "12px IBM Plex Mono, monospace";
    ctx.fillText("x →", pad.l + W - 28, cssH - 14);
    ctx.save();
    ctx.translate(16, pad.t + H / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("imaginary time τ", 0, 0);
    ctx.restore();
    ctx.fillStyle = "#c9bfae";
    ctx.fillText(`pinned at x = ${state.xb}, N = ${state.N}, T = ${state.T}`, pad.l, 14);
  }

  function drawDensity() {
    const canvas = els.densityCanvas;
    if (!canvas || !state.density) return;
    const ctx = canvas.getContext("2d");
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const cssW = canvas.clientWidth || 640;
    const cssH = canvas.clientHeight || 360;
    canvas.width = Math.floor(cssW * dpr);
    canvas.height = Math.floor(cssH * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const pad = { l: 48, r: 18, t: 22, b: 40 };
    const W = cssW - pad.l - pad.r;
    const H = cssH - pad.t - pad.b;
    const xs = state.density.x;
    const ys = state.density.psi2;
    const analytic = xs.map((x) => state.density.analytic(x));
    const ymax = Math.max(...ys, ...analytic, 0.01) * 1.12;
    const xmin = xs[0];
    const xmax = xs[xs.length - 1];

    const xMap = (x) => pad.l + ((x - xmin) / (xmax - xmin)) * W;
    const yMap = (y) => pad.t + H - (y / ymax) * H;

    ctx.clearRect(0, 0, cssW, cssH);
    ctx.fillStyle = "#16120e";
    ctx.fillRect(0, 0, cssW, cssH);

    // axes
    ctx.strokeStyle = "rgba(196,163,90,0.25)";
    ctx.beginPath();
    ctx.moveTo(pad.l, pad.t);
    ctx.lineTo(pad.l, pad.t + H);
    ctx.lineTo(pad.l + W, pad.t + H);
    ctx.stroke();

    // analytic
    ctx.strokeStyle = "rgba(196,163,90,0.85)";
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    analytic.forEach((y, i) => {
      const X = xMap(xs[i]);
      const Y = yMap(y);
      if (i === 0) ctx.moveTo(X, Y);
      else ctx.lineTo(X, Y);
    });
    ctx.stroke();

    // MC
    ctx.strokeStyle = "rgba(168,192,176,0.95)";
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    ys.forEach((y, i) => {
      const X = xMap(xs[i]);
      const Y = yMap(y);
      if (i === 0) ctx.moveTo(X, Y);
      else ctx.lineTo(X, Y);
    });
    ctx.stroke();

    // marker at x=1
    const y1 = F.interpLinear(xs, ys, 1);
    ctx.fillStyle = "#efe7d8";
    ctx.beginPath();
    ctx.arc(xMap(1), yMap(y1), 3.5, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#8f8678";
    ctx.font = "12px IBM Plex Mono, monospace";
    ctx.fillText("x", pad.l + W - 12, cssH - 14);
    ctx.fillText("|φ|²", 8, pad.t + 8);
    ctx.fillStyle = "#c4a35a";
    ctx.fillText("analytic", pad.l + 8, pad.t + 18);
    ctx.fillStyle = "#a8c0b0";
    ctx.fillText("direct MC (browser)", pad.l + 8, pad.t + 34);
  }

  function drawNT() {
    const canvas = els.ntCanvas;
    if (!canvas || !state.exactAlpha.length) return;
    const ctx = canvas.getContext("2d");
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const cssW = canvas.clientWidth || 480;
    const cssH = 280;
    canvas.width = Math.floor(cssW * dpr);
    canvas.height = Math.floor(cssH * dpr);
    canvas.style.height = cssH + "px";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const pad = { l: 44, r: 16, t: 20, b: 36 };
    const W = cssW - pad.l - pad.r;
    const H = cssH - pad.t - pad.b;
    const rows = state.exactAlpha;
    const Ns = [...new Set(rows.map((r) => Number(r.N)))].sort((a, b) => a - b);
    const Ts = [...new Set(rows.map((r) => Number(r.T)))].sort((a, b) => a - b);
    const vals = rows.map((r) => Number(r.psi2_1));
    const vmin = Math.min(...vals);
    const vmax = Math.max(...vals);

    ctx.fillStyle = "#16120e";
    ctx.fillRect(0, 0, cssW, cssH);

    const cellW = W / Ts.length;
    const cellH = H / Ns.length;

    Ns.forEach((N, ni) => {
      Ts.forEach((T, ti) => {
        const row = rows.find((r) => Number(r.N) === N && Number(r.T) === T);
        if (!row) return;
        const v = Number(row.psi2_1);
        const t = (v - vmin) / (vmax - vmin + 1e-12);
        const selected = N === state.N && Math.abs(T - state.T) < 1e-9;
        ctx.fillStyle = `rgba(196, 163, 90, ${0.12 + t * 0.55})`;
        ctx.fillRect(pad.l + ti * cellW + 2, pad.t + ni * cellH + 2, cellW - 4, cellH - 4);
        if (selected) {
          ctx.strokeStyle = "#efe7d8";
          ctx.lineWidth = 1.5;
          ctx.strokeRect(pad.l + ti * cellW + 2, pad.t + ni * cellH + 2, cellW - 4, cellH - 4);
        }
        ctx.fillStyle = "#efe7d8";
        ctx.font = "11px IBM Plex Mono, monospace";
        ctx.fillText(v.toFixed(3), pad.l + ti * cellW + 8, pad.t + ni * cellH + cellH * 0.58);
      });
    });

    ctx.fillStyle = "#8f8678";
    ctx.font = "11px IBM Plex Mono, monospace";
    Ts.forEach((T, ti) => {
      ctx.fillText("T=" + T, pad.l + ti * cellW + 8, cssH - 14);
    });
    Ns.forEach((N, ni) => {
      ctx.fillText("N=" + N, 6, pad.t + ni * cellH + cellH * 0.58);
    });

    if (els.ntReadout) {
      const hit = rows.find(
        (r) => Number(r.N) === state.N && Math.abs(Number(r.T) - state.T) < 1e-9
      );
      if (hit) {
        els.ntReadout.textContent =
          `Exact HO lattice |φ(1)|² = ${Number(hit.psi2_1).toFixed(5)} at N=${state.N}, T=${state.T} (suite exact_alpha). Analytic continuum = 0.20755.`;
      } else {
        els.ntReadout.textContent =
          `No exact_alpha entry for N=${state.N}, T=${state.T}. Move sliders onto a tabulated lattice cell.`;
      }
    }
  }

  function updateCaption() {
    if (!els.caption) return;
    if (state.tab === "bridges") {
      els.caption.textContent =
        "Brownian bridges pinned at x = 1, reweighted by the potential. Opacity follows the potential weight — the same bridge_potential_weight used in the C++ engine.";
    } else {
      els.caption.textContent =
        "Live direct Monte Carlo |φ(x)|² from the browser port of direct_ground_density (fewer samples than the suite). Gold: analytic HO ground state. Marker: estimate at x = 1.";
    }
  }

  function refresh(heavy) {
    syncOutputs();
    sampleBridges();
    if (heavy) estimateDensity();
    if (state.tab === "bridges") drawBridges();
    else drawDensity();
    drawNT();
    updateCaption();
  }

  let debounceTimer = null;
  function onControlChange(heavy) {
    state.N = Number(els.N.value);
    state.T = Number(els.T.value);
    state.hbar = Number(els.hbar.value);
    if (els.pot) state.pot = els.pot.value;
    syncOutputs();
    sampleBridges();
    drawBridges();
    drawNT();
    setStatus("Sampling bridges…");
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => {
      if (heavy !== false) {
        setStatus("Estimating density…");
        estimateDensity();
        if (state.tab === "density") drawDensity();
        setStatus("Ready — browser MC (not suite sample size).");
      }
    }, 180);
  }

  function bindTabs() {
    document.querySelectorAll("[data-tab]").forEach((btn) => {
      btn.addEventListener("click", () => activateTab(btn.getAttribute("data-tab")));
    });
  }

  function applyPreset(p) {
    const btn = document.querySelector(`[data-preset="${p}"]`);
    document.querySelectorAll("[data-preset]").forEach((b) => b.classList.remove("active"));
    if (btn) btn.classList.add("active");
    if (p === "milestone") {
      els.N.value = 7;
      els.T.value = 4;
      els.hbar.value = 1;
      if (els.pot) els.pot.value = "ho";
    } else if (p === "fine") {
      els.N.value = 16;
      els.T.value = 8;
      els.hbar.value = 1;
      if (els.pot) els.pot.value = "ho";
    } else if (p === "smallh") {
      els.N.value = 12;
      els.T.value = 6;
      els.hbar.value = 0.5;
      if (els.pot) els.pot.value = "ho";
    } else if (p === "dwell") {
      els.N.value = 16;
      els.T.value = 8;
      els.hbar.value = 1;
      if (els.pot) els.pot.value = "doublewell";
    }
    onControlChange(true);
  }

  function bindPresets() {
    document.querySelectorAll("[data-preset]").forEach((btn) => {
      btn.addEventListener("click", () => applyPreset(btn.getAttribute("data-preset")));
    });
  }

  function activateTab(name) {
    state.tab = name;
    document.querySelectorAll("[data-tab]").forEach((b) => {
      const on = b.getAttribute("data-tab") === name;
      b.classList.toggle("active", on);
      b.setAttribute("aria-selected", on ? "true" : "false");
    });
    document.querySelectorAll(".stage-wrap").forEach((w) => {
      w.classList.toggle("active", w.id === "stage-" + name);
    });
    updateCaption();
    if (name === "bridges") drawBridges();
    else {
      if (!state.density) estimateDensity();
      drawDensity();
    }
  }

  function bindKeys() {
    window.addEventListener("keydown", (e) => {
      if (e.target && /INPUT|SELECT|TEXTAREA/.test(e.target.tagName)) return;
      const map = { "1": "milestone", "2": "fine", "3": "smallh", "4": "dwell" };
      if (map[e.key]) applyPreset(map[e.key]);
      if (e.key === "b" || e.key === "B") activateTab("bridges");
      if (e.key === "d" || e.key === "D") activateTab("density");
    });
  }

  async function loadCSV(path) {
    const res = await fetch(path);
    const text = await res.text();
    const lines = text.trim().split(/\r?\n/);
    const headers = lines[0].split(",");
    return lines.slice(1).map((line) => {
      const cols = line.split(",");
      const obj = {};
      headers.forEach((h, i) => {
        obj[h.trim()] = cols[i];
      });
      return obj;
    });
  }

  async function init() {
    bindTabs();
    bindPresets();
    bindKeys();
    ["N", "T", "hbar"].forEach((k) => {
      if (els[k]) els[k].addEventListener("input", () => onControlChange(true));
    });
    if (els.pot) els.pot.addEventListener("change", () => onControlChange(true));

    try {
      state.exactAlpha = await loadCSV("data/exact_alpha.csv");
      const summaryRows = await loadCSV("data/summary.csv");
      summaryRows.forEach((r) => {
        state.suiteSummary[r.key] = Number(r.value);
      });
      const suiteMC = document.getElementById("suite-mc");
      if (suiteMC && state.suiteSummary.direct_psi2_1) {
        suiteMC.textContent =
          Number(state.suiteSummary.direct_psi2_1).toFixed(4) +
          " ± " +
          Number(state.suiteSummary.direct_psi2_1_err).toFixed(4);
      }
    } catch (e) {
      console.warn("Suite data load failed", e);
      setStatus("Suite CSV load failed — interactive sampler still runs.");
    }

    refresh(true);
    setStatus("Ready — browser MC (indicative; suite uses more samples). Keys: 1–4 presets, B/D tabs.");
    window.addEventListener("resize", () => {
      if (state.tab === "bridges") drawBridges();
      else drawDensity();
      drawNT();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
