import json
import os
import subprocess
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXE = os.path.join(ROOT, "build", "feynman.exe")
PORT = 8765

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Path integral engine</title>
<style>
  :root { --bg:#141820; --panel:#1c2230; --ink:#f3efe6; --muted:#b7b1a6; --line:#343c4e; --teal:#1f8a84; --copper:#d7a441; }
  * { box-sizing: border-box; }
  body { margin:0; background:var(--bg); color:var(--ink); font-family: Georgia, "Palatino Linotype", serif; }
  header { padding: 22px 28px 8px; }
  h1 { font-weight: normal; font-size: 28px; margin: 0; letter-spacing: 0.01em; }
  header p { color: var(--muted); margin: 6px 0 0; font-size: 14px; }
  main { display: grid; grid-template-columns: 280px 1fr; gap: 18px; padding: 16px 28px 28px; }
  .panel { background: var(--panel); border: 1px solid var(--line); border-radius: 10px; padding: 14px; }
  label { display:block; color: var(--muted); font-size: 12px; margin: 10px 0 4px; }
  input, select { width: 100%; background: #121722; color: var(--ink); border: 1px solid var(--line); border-radius: 6px; padding: 7px 8px; font-family: Consolas, monospace; font-size: 13px; }
  .presets { display:flex; flex-wrap: wrap; gap: 6px; margin-bottom: 8px; }
  button { background: #243044; color: var(--ink); border: 1px solid var(--line); border-radius: 999px; padding: 6px 10px; cursor: pointer; font-family: Georgia, serif; font-size: 13px; }
  button.run { width: 100%; margin-top: 14px; background: var(--teal); border: none; color: #f6fffe; border-radius: 8px; padding: 10px; font-size: 16px; }
  button:disabled { opacity: 0.5; cursor: wait; }
  .nums { display:grid; grid-template-columns: 1fr 1fr; gap: 8px; margin-top: 12px; }
  .num { background:#121722; border-radius: 8px; padding: 8px; }
  .num b { display:block; font-family: Consolas, monospace; font-size: 14px; font-weight: normal; }
  .num span { color: var(--muted); font-size: 11px; }
  canvas { width: 100%; background: #10141c; border-radius: 8px; display:block; }
  h2 { font-weight: normal; font-size: 16px; margin: 0 0 8px; }
  .stack { display:flex; flex-direction: column; gap: 14px; }
  #status { min-height: 1.2em; color: var(--copper); font-size: 13px; margin-top: 8px; }
</style>
</head>
<body>
<header>
  <h1>Path integral engine</h1>
  <p>One dimension, imaginary time. The library underneath is the same C++ code as the milestone.</p>
</header>
<main>
  <section class="panel">
    <div class="presets" id="presets"></div>
    <label>potential</label>
    <select id="pot">
      <option value="ho">harmonic oscillator</option>
      <option value="anharm">anharmonic</option>
      <option value="doublewell">double well</option>
      <option value="morse">Morse</option>
    </select>
    <label>slices N</label><input id="N" value="36">
    <label>imaginary time T</label><input id="T" value="8">
    <label>mass</label><input id="mass" value="1">
    <label>hbar</label><input id="hbar" value="1">
    <label>sweeps</label><input id="sweeps" value="6000">
    <label>x range</label>
    <div style="display:flex; gap:8px">
      <input id="xmin" value="-6"><input id="xmax" value="6">
    </div>
    <button class="run" id="go">Run</button>
    <div id="status">Ready.</div>
    <div class="nums">
      <div class="num"><span>path integral at x = 1</span><b id="pi">–</b></div>
      <div class="num"><span>Schrödinger at x = 1</span><b id="se">–</b></div>
      <div class="num"><span>ground energy</span><b id="en">–</b></div>
      <div class="num"><span>acceptance</span><b id="acc">–</b></div>
    </div>
  </section>
  <section class="stack">
    <div class="panel">
      <h2>Ground-state density</h2>
      <canvas id="density" width="900" height="360"></canvas>
    </div>
    <div class="panel">
      <h2>A few paths from the chain</h2>
      <canvas id="paths" width="900" height="280"></canvas>
    </div>
  </section>
</main>
<script>
const presets = {
  "Milestone": {pot:"ho", N:7, T:4, mass:1, hbar:1, sweeps:8000, xmin:-4, xmax:4},
  "Harmonic": {pot:"ho", N:40, T:10, mass:1, hbar:1, sweeps:7000, xmin:-5, xmax:5},
  "Small hbar": {pot:"ho", N:40, T:10, mass:1, hbar:0.45, sweeps:7000, xmin:-3, xmax:3},
  "Double well": {pot:"doublewell", N:48, T:10, mass:1, hbar:1, sweeps:9000, xmin:-3.2, xmax:3.2},
  "Morse": {pot:"morse", N:48, T:10, mass:1, hbar:1, sweeps:8000, xmin:-1.5, xmax:5.5}
};
const box = document.getElementById("presets");
Object.keys(presets).forEach(name => {
  const b = document.createElement("button");
  b.textContent = name;
  b.onclick = () => apply(presets[name]);
  box.appendChild(b);
});
function apply(p) { for (const k in p) document.getElementById(k).value = p[k]; }
function num(id) { return document.getElementById(id).value; }
function drawDensity(data) {
  const c = document.getElementById("density");
  const g = c.getContext("2d");
  const w = c.width, h = c.height;
  g.clearRect(0,0,w,h);
  const pad = {l:46, r:16, t:16, b:32};
  const xs = data.x.concat(data.se_x);
  const ys = data.psi2.concat(data.se_psi2);
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const ymax = Math.max(...ys) * 1.12;
  const X = x => pad.l + (x - xmin) / (xmax - xmin) * (w - pad.l - pad.r);
  const Y = y => h - pad.b - y / ymax * (h - pad.t - pad.b);
  g.strokeStyle = "#343c4e"; g.lineWidth = 1;
  g.beginPath(); g.moveTo(pad.l, pad.t); g.lineTo(pad.l, h-pad.b); g.lineTo(w-pad.r, h-pad.b); g.stroke();
  g.fillStyle = "#8ea0b8"; g.font = "12px Georgia";
  g.fillText("x", w-pad.r-8, h-8);
  g.fillText("|phi|^2", 8, 18);
  function line(xs, ys, color, width) {
    g.beginPath();
    xs.forEach((x,i) => i ? g.lineTo(X(x), Y(ys[i])) : g.moveTo(X(x), Y(ys[i])));
    g.strokeStyle = color; g.lineWidth = width; g.stroke();
  }
  g.globalAlpha = 0.25; g.fillStyle = "#1f8a84";
  g.beginPath();
  data.x.forEach((x,i) => {
    const y0 = Y(Math.max(0, data.psi2[i] - data.err[i]));
    const y1 = Y(data.psi2[i] + data.err[i]);
    i ? g.lineTo(X(x), y1) : g.moveTo(X(x), y1);
  });
  for (let i = data.x.length - 1; i >= 0; i--) g.lineTo(X(data.x[i]), Y(Math.max(0, data.psi2[i] - data.err[i])));
  g.closePath(); g.fill(); g.globalAlpha = 1;
  line(data.se_x, data.se_psi2, "#f3efe6", 1.6);
  line(data.x, data.psi2, "#1f8a84", 2);
}
function drawPaths(data) {
  const c = document.getElementById("paths");
  const g = c.getContext("2d");
  const w = c.width, h = c.height;
  g.clearRect(0,0,w,h);
  const pad = {l:46, r:16, t:12, b:28};
  let ymin = Infinity, ymax = -Infinity;
  data.paths.forEach(p => p.forEach(v => { ymin = Math.min(ymin, v); ymax = Math.max(ymax, v); }));
  if (!isFinite(ymin)) return;
  const span = (ymax - ymin) || 1;
  ymin -= 0.08*span; ymax += 0.08*span;
  const X = (i, n) => pad.l + i / (n - 1) * (w - pad.l - pad.r);
  const Y = y => h - pad.b - (y - ymin) / (ymax - ymin) * (h - pad.t - pad.b);
  g.strokeStyle = "#343c4e";
  g.beginPath(); g.moveTo(pad.l, pad.t); g.lineTo(pad.l, h-pad.b); g.lineTo(w-pad.r, h-pad.b); g.stroke();
  data.paths.forEach((p, k) => {
    g.beginPath();
    p.forEach((v,i) => i ? g.lineTo(X(i, p.length), Y(v)) : g.moveTo(X(i, p.length), Y(v)));
    const t = data.paths.length === 1 ? 0.5 : k / (data.paths.length - 1);
    g.strokeStyle = `rgb(${80+Math.floor(160*t)}, ${150-Math.floor(40*t)}, ${90+Math.floor(40*(1-t))})`;
    g.lineWidth = 1.1; g.stroke();
  });
  g.fillStyle = "#8ea0b8"; g.font = "12px Georgia";
  g.fillText("slice", w-48, h-8);
  g.fillText("x", 8, 16);
}
document.getElementById("go").onclick = async () => {
  const btn = document.getElementById("go");
  btn.disabled = true;
  document.getElementById("status").textContent = "Running the chain...";
  const body = {
    pot: num("pot"), N: +num("N"), T: +num("T"), mass: +num("mass"),
    hbar: +num("hbar"), sweeps: +num("sweeps"), xmin: +num("xmin"), xmax: +num("xmax")
  };
  try {
    const res = await fetch("/api/run", {method:"POST", headers:{"Content-Type":"application/json"}, body: JSON.stringify(body)});
    const data = await res.json();
    if (!data.ok) throw new Error(data.error || "run failed");
    drawDensity(data);
    drawPaths(data);
    document.getElementById("pi").textContent = data.psi2_at_1.toFixed(4) + " ± " + data.psi2_at_1_err.toFixed(4);
    document.getElementById("se").textContent = data.se_psi2_at_1.toFixed(4);
    document.getElementById("en").textContent = data.energy.toFixed(4);
    document.getElementById("acc").textContent = (100*data.acceptance).toFixed(0) + "%";
    document.getElementById("status").textContent = data.label + "   N=" + data.N + "   T=" + data.T + "   eps=" + data.eps.toFixed(3);
  } catch (err) {
    document.getElementById("status").textContent = String(err);
  }
  btn.disabled = false;
};
</script>
</body>
</html>
"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        return

    def _send(self, code, body, content_type):
        raw = body if isinstance(body, bytes) else body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            self._send(200, PAGE, "text/html; charset=utf-8")
            return
        self._send(404, "not found", "text/plain")

    def do_POST(self):
        if self.path != "/api/run":
            self._send(404, "not found", "text/plain")
            return
        length = int(self.headers.get("Content-Length", "0"))
        try:
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
        except json.JSONDecodeError:
            self._send(400, json.dumps({"ok": False, "error": "bad json"}), "application/json")
            return
        pot = str(payload.get("pot", "ho"))
        if pot not in {"ho", "anharm", "doublewell", "morse", "quartic"}:
            self._send(400, json.dumps({"ok": False, "error": "unknown potential"}), "application/json")
            return
        if not os.path.isfile(EXE):
            self._send(500, json.dumps({"ok": False, "error": "build/feynman.exe is missing"}), "application/json")
            return

        def grab(key, default, lo, hi):
            try:
                value = float(payload.get(key, default))
            except (TypeError, ValueError):
                value = default
            return min(hi, max(lo, value))

        N = int(grab("N", 36, 4, 96))
        sweeps = int(grab("sweeps", 6000, 400, 20000))
        args = [
            EXE, "run",
            "--pot", pot,
            "--N", str(N),
            "--T", str(grab("T", 8, 0.5, 40)),
            "--mass", str(grab("mass", 1, 0.05, 20)),
            "--hbar", str(grab("hbar", 1, 0.05, 10)),
            "--sweeps", str(sweeps),
            "--xmin", str(grab("xmin", -6, -20, 0)),
            "--xmax", str(grab("xmax", 6, 0.2, 20)),
            "--chains", "2",
            "--bins", "120",
            "--paths", "7",
            "--seed", "7",
        ]
        try:
            proc = subprocess.run(args, capture_output=True, text=True, timeout=120, cwd=ROOT)
        except subprocess.TimeoutExpired:
            self._send(500, json.dumps({"ok": False, "error": "timed out"}), "application/json")
            return
        if proc.returncode != 0 or not proc.stdout.strip():
            err = (proc.stderr or proc.stdout or "engine failed").strip()
            self._send(500, json.dumps({"ok": False, "error": err[-400:]}), "application/json")
            return
        self._send(200, proc.stdout, "application/json")


def main():
    server = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    url = f"http://127.0.0.1:{PORT}"
    print(url, flush=True)
    webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()


if __name__ == "__main__":
    main()
