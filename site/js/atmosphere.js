/**
 * Subtle imaginary-time lattice / worldline atmosphere.
 * Never overpowers content; respects prefers-reduced-motion.
 */
(function () {
  "use strict";

  const canvas = document.getElementById("atmosphere");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  let w = 0;
  let h = 0;
  let dpr = 1;
  const paths = [];
  const PATH_COUNT = window.matchMedia("(max-width: 640px)").matches ? 8 : 14;

  function resize() {
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    w = window.innerWidth;
    h = window.innerHeight;
    canvas.width = Math.floor(w * dpr);
    canvas.height = Math.floor(h * dpr);
    canvas.style.width = w + "px";
    canvas.style.height = h + "px";
    canvas.style.opacity = w < 700 ? "0.38" : "0.62";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function seedPaths() {
    paths.length = 0;
    for (let i = 0; i < PATH_COUNT; i++) {
      const beads = 18 + Math.floor(Math.random() * 10);
      const amp = 40 + Math.random() * 90;
      const phase = Math.random() * Math.PI * 2;
      const drift = (Math.random() - 0.5) * 0.15;
      const x0 = (Math.random() * 0.7 + 0.15) * w;
      paths.push({ beads, amp, phase, drift, x0, weight: 0.25 + Math.random() * 0.55 });
    }
  }

  function drawLattice() {
    const N = 8;
    const Trows = 12;
    ctx.save();
    ctx.strokeStyle = "rgba(196, 163, 90, 0.045)";
    ctx.lineWidth = 1;
    for (let i = 0; i <= N; i++) {
      const x = (w * 0.12) + (i / N) * w * 0.76;
      ctx.beginPath();
      ctx.moveTo(x, h * 0.08);
      ctx.lineTo(x, h * 0.92);
      ctx.stroke();
    }
    for (let j = 0; j <= Trows; j++) {
      const y = (h * 0.08) + (j / Trows) * h * 0.84;
      ctx.beginPath();
      ctx.moveTo(w * 0.12, y);
      ctx.lineTo(w * 0.88, y);
      ctx.stroke();
    }
    ctx.restore();
  }

  function drawPaths(t) {
    for (let p = 0; p < paths.length; p++) {
      const path = paths[p];
      ctx.beginPath();
      for (let j = 0; j <= path.beads; j++) {
        const u = j / path.beads;
        const y = h * 0.1 + u * h * 0.8;
        const sway =
          Math.sin(path.phase + u * 4.2 + t * 0.00035 * (1 + path.drift * 4)) * path.amp +
          Math.sin(path.phase * 0.7 + u * 9 + t * 0.0002) * path.amp * 0.28;
        const x = path.x0 + sway + path.drift * u * 40;
        if (j === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.strokeStyle = `rgba(232, 220, 196, ${0.07 + path.weight * 0.1})`;
      ctx.lineWidth = 1.15;
      ctx.stroke();

      // beads
      for (let j = 0; j <= path.beads; j += 2) {
        const u = j / path.beads;
        const y = h * 0.1 + u * h * 0.8;
        const sway =
          Math.sin(path.phase + u * 4.2 + t * 0.00035 * (1 + path.drift * 4)) * path.amp +
          Math.sin(path.phase * 0.7 + u * 9 + t * 0.0002) * path.amp * 0.28;
        const x = path.x0 + sway + path.drift * u * 40;
        ctx.beginPath();
        ctx.arc(x, y, 1.5, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(196, 163, 90, ${0.12 + path.weight * 0.16})`;
        ctx.fill();
      }
    }
  }

  let frame = 0;
  let last = 0;

  function paint(ts) {
    if (!last) last = ts;
    const t = reduced ? 0 : ts;
    ctx.clearRect(0, 0, w, h);

    // soft vignette plane
    const g = ctx.createRadialGradient(w * 0.65, h * 0.2, 40, w * 0.5, h * 0.5, Math.max(w, h) * 0.75);
    g.addColorStop(0, "rgba(42, 34, 25, 0.35)");
    g.addColorStop(0.55, "rgba(18, 16, 14, 0.15)");
    g.addColorStop(1, "rgba(8, 7, 5, 0)");
    ctx.fillStyle = g;
    ctx.fillRect(0, 0, w, h);

    drawLattice();
    drawPaths(t);

    if (!reduced) {
      frame = requestAnimationFrame(paint);
    }
  }

  function init() {
    resize();
    seedPaths();
    if (reduced) {
      paint(0);
    } else {
      frame = requestAnimationFrame(paint);
    }
  }

  window.addEventListener("resize", () => {
    cancelAnimationFrame(frame);
    resize();
    seedPaths();
    last = 0;
    if (reduced) paint(0);
    else frame = requestAnimationFrame(paint);
  });

  init();
})();
