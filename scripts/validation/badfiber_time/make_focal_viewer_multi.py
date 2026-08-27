"""
make_focal_viewer_multi.py
Interactive focal-plane viewer showing multiple petals simultaneously.
Each petal gets its own accent colour; bad fibers are orange regardless of petal.

Usage:
  python make_focal_viewer_multi.py --petals 0 7
"""
import argparse, json, os, datetime

MJD_EPOCH = datetime.date(1858, 11, 17)
def date_to_mjd(yyyymmdd):
    s = str(yyyymmdd)
    d = datetime.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return (d - MJD_EPOCH).days

CCD_SWAPS_ALL = {
    0: [(date_to_mjd(20260106), 'r')],
    1: [(date_to_mjd(20220613), 'b'), (date_to_mjd(20230727), 'b'),
        (date_to_mjd(20231129), 'z'), (date_to_mjd(20241119), 'r')],
    2: [(date_to_mjd(20201113), 'r')],
    3: [(date_to_mjd(20250514), 'z')],
    4: [(date_to_mjd(20221108), 'b'), (date_to_mjd(20221108), 'r')],
    5: [(date_to_mjd(20210622), 'b'), (date_to_mjd(20221108), 'z')],
    6: [],
    7: [(date_to_mjd(20240724), 'r'), (date_to_mjd(20251007), 'r'),
        (date_to_mjd(20251007), 'z')],
    8: [(date_to_mjd(20230507), 'b')],
    9: [(date_to_mjd(20210615), 'r'), (date_to_mjd(20260429), 'r')],
}
parser = argparse.ArgumentParser()
parser.add_argument('--petals', type=int, nargs='+', default=[0, 7])
args = parser.parse_args()
PETALS = args.petals

DATADIR  = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/all-fibers-vs-time'
label    = '_'.join(str(p) for p in PETALS)
out_path = f'{DATADIR}/lrg_petals{label}_focal_viewer.html'

# load each petal JSON
petal_data = {}
for p in PETALS:
    jp = f'{DATADIR}/lrg_petal{p}_focal_data.json'
    with open(jp) as f:
        petal_data[p] = json.load(f)

# all_positions is the same for every petal — just use the first
all_positions = petal_data[PETALS[0]]['all_positions']
bins_mjd  = petal_data[PETALS[0]]['bins_mjd']
survey_fr = petal_data[PETALS[0]]['survey_fr']
survey_er = petal_data[PETALS[0]]['survey_er']

# merge fibers dict: key = str(fiberId), value = fiber record + 'petal' field
fibers_merged = {}
for p in PETALS:
    for fid, rec in petal_data[p]['fibers'].items():
        r2 = dict(rec)
        r2['petal'] = p
        fibers_merged[fid] = r2

combined = {
    'petals':       PETALS,
    'all_positions': all_positions,
    'bins_mjd':     bins_mjd,
    'survey_fr':    survey_fr,
    'survey_er':    survey_er,
    'fibers':       fibers_merged,
}
raw_json = json.dumps(combined, separators=(',', ':'))

# one accent colour per petal (index into PETALS list)
PETAL_COLORS = ['#7878a0', '#4a9a6a', '#6a9aca', '#ca7a4a',
                '#9a6aca', '#ca4a6a', '#4acaaa', '#caaa4a',
                '#aa4aca', '#4a7aca']

HTML = """<meta charset="UTF-8">
<title>DESI Focal Plane — LRG redshift failures</title>
<style>
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
body{background:#13131a;color:#c8c8d8;font-family:ui-monospace,'SF Mono',Menlo,monospace;
     display:flex;flex-direction:column;height:100vh;overflow:hidden}
header{padding:10px 18px 8px;border-bottom:1px solid #252530;display:flex;align-items:baseline;
       gap:12px;flex-shrink:0;flex-wrap:wrap}
.htitle{font-size:.82rem;letter-spacing:.08em;color:#e0e0f0;font-weight:700}
.hmeta{font-size:.70rem;color:#9898b8;font-weight:600}
.pleg{display:flex;gap:10px;margin-left:auto}
.credit{font-size:.60rem;color:#7070a0;letter-spacing:.04em;margin-left:16px;align-self:center;white-space:nowrap}
.pchip{font-size:.62rem;padding:2px 7px;border-radius:2px;border:1px solid}
main{display:flex;flex:1;overflow:hidden}
#fp-wrap{flex-shrink:0;display:flex;align-items:center;justify-content:center;
          padding:14px 10px 10px 14px;border-right:1px solid #252530}
#fp-canvas{cursor:crosshair;display:block}
#tooltip{position:fixed;pointer-events:none;background:#1e1e2c;border:1px solid #38384a;
          border-radius:3px;padding:5px 9px;font-size:.65rem;color:#dde;opacity:0;
          transition:opacity .1s;z-index:99;white-space:nowrap}
#plot-wrap{flex:1;display:flex;flex-direction:column;padding:14px 16px 10px;overflow:hidden;min-width:0}
#fiber-title{font-size:.95rem;color:#e0e0f0;margin-bottom:8px;letter-spacing:.03em;min-height:1.4em}
#plot-canvas{display:block;max-width:100%}
#legend{display:flex;gap:22px;margin-top:10px;font-size:.78rem;color:#6a6a80}
.leg{display:flex;align-items:center;gap:5px}
.lsw{width:20px;height:2px}
.lsw-sv{background:rgba(190,190,215,.45);height:6px;border-radius:1px}
#goto-box{display:flex;align-items:center;gap:8px;padding:6px 0 2px;flex-basis:100%}
#goto-input{background:#1e1e2c;border:1px solid #38384a;color:#e0e0f0;
             font-family:inherit;font-size:.85rem;padding:5px 10px;width:90px;
             border-radius:3px;outline:none}
#goto-input:focus{border-color:#5a5a7a}
#goto-input.error{border-color:#c03030}
#goto-btn{background:#252535;border:1px solid #38384a;color:#c8c8d8;
           font-family:inherit;font-size:.85rem;padding:5px 14px;
           border-radius:3px;cursor:pointer}
#goto-btn:hover{background:#32324a}
#goto-label{font-size:.72rem;color:#55556a;letter-spacing:.05em;text-transform:uppercase}
.lsw-md{background:repeating-linear-gradient(90deg,#aaa 0,#aaa 4px,transparent 4px,transparent 7px);height:2px;border:none}
.save-btn{background:#1e2830;border:1px solid #2e4050;color:#88aac0;
           font-family:inherit;font-size:.72rem;padding:4px 10px;
           border-radius:3px;cursor:pointer;letter-spacing:.03em}
.save-btn:hover{background:#253545;color:#aaccdd}
#img-overlay{display:none;position:fixed;inset:0;background:rgba(0,0,0,.85);
             z-index:200;flex-direction:column;align-items:center;justify-content:center;gap:10px}
#img-overlay.show{display:flex}
#img-overlay img{max-width:95vw;max-height:88vh;border:1px solid #38384a}
#img-overlay-hint{color:#9898b8;font-size:.72rem;letter-spacing:.05em}
#img-overlay-close{color:#6a6a80;font-size:.72rem;cursor:pointer;text-decoration:underline}
</style>

<div id="tooltip"></div>
<div id="img-overlay">
  <img id="img-overlay-img" src="" alt="saved image">
  <div id="img-overlay-hint">Right-click → Save image as…</div>
  <div id="img-overlay-close">close</div>
</div>
<header>
  <span class="htitle">DESI Focal Plane — LRG Redshift failures</span>
  <span class="hmeta">hover to identify · click for time plot · orange = bad fiber</span>
  <div class="pleg" id="pleg"></div>
  <span class="credit">J. Rohlf &amp; Claude Sonnet 4.6 (2026)</span>
  <div id="goto-box">
    <span id="goto-label">Fiber</span>
    <input id="goto-input" type="number" placeholder="number" min="0" max="4999">
    <button id="goto-btn">Go</button>
  </div>
</header>
<main>
  <div id="fp-wrap">
    <canvas id="fp-canvas"></canvas>
  </div>
  <div id="plot-wrap">
    <div id="fiber-title">← click a highlighted fiber</div>
    <canvas id="plot-canvas"></canvas>
    <div id="legend">
      <div class="leg"><div class="lsw lsw-sv"></div>survey avg ± 1σ</div>
      <div class="leg"><div class="lsw lsw-md"></div>model prediction</div>
      <div class="leg"><div class="lsw" id="fb-swatch" style="background:#c03030"></div>fiber fail rate</div>
      <div class="leg"><div class="lsw lsw-md" style="background:repeating-linear-gradient(90deg,#4169e1 0,#4169e1 4px,transparent 4px,transparent 7px)"></div>CCD b swap</div>
      <div class="leg"><div class="lsw lsw-md" style="background:repeating-linear-gradient(90deg,#dc143c 0,#dc143c 4px,transparent 4px,transparent 7px)"></div>CCD r swap</div>
      <div class="leg"><div class="lsw lsw-md" style="background:repeating-linear-gradient(90deg,#9932cc 0,#9932cc 4px,transparent 4px,transparent 7px)"></div>CCD z swap</div>
    </div>
    <div style="margin-top:8px">
      <button class="save-btn" id="save-plot-btn">Save plot</button>
    </div>
  </div>
</main>

<script>
const DATA         = __JSON__;
const PETALS       = DATA.petals;
const PETAL_COLORS = __COLORS__;
const CCD_SWAPS    = __CCD_SWAPS__;
const CCD_COLOR    = {'b': '#4169e1', 'r': '#dc143c', 'z': '#9932cc'};


const fpCanvas   = document.getElementById('fp-canvas');
const fpCtx      = fpCanvas.getContext('2d');
const plotCanvas = document.getElementById('plot-canvas');
const plotCtx    = plotCanvas.getContext('2d');
const tooltip    = document.getElementById('tooltip');
const fiberTitle = document.getElementById('fiber-title');

const FP_SIZE = Math.min(window.innerHeight - 80, 520);
fpCanvas.width = fpCanvas.height = FP_SIZE;
const SCALE = FP_SIZE / 900;
const cx = FP_SIZE / 2, cy = FP_SIZE / 2;
function mm2px(x, y) { return [cx + x * SCALE, cy - y * SCALE]; }

const fibPos = {};
for (const [fid, xy] of Object.entries(DATA.all_positions)) {
  fibPos[fid] = mm2px(xy[0], xy[1]);
}

const activeFibers = Object.keys(DATA.fibers).map(Number);
const badSet = new Set(activeFibers.filter(f => DATA.fibers[f].is_bad));

// colour per fiber — all active fibers green, bad fibers orange
function fiberColor(fib, hover, selected) {
  if (selected)        return '#ffffff';
  if (hover)           return '#ffffffcc';
  if (badSet.has(fib)) return '#e8910a';
  return '#4aaa6a';
}

let selectedFiber = null;
let plotDots = [];   // [{x, y, mjd_lo, mjd_hi, fr}] for current time plot

const MJD_EPOCH_MS = Date.UTC(1858, 10, 17);
function mjdToDate(mjd) {
  const d = new Date(MJD_EPOCH_MS + mjd * 86400000);
  return d.toLocaleDateString('en-US', { month: 'short', year: 'numeric', timeZone: 'UTC' });
}

function drawFP(hoverId) {
  const ctx = fpCtx;
  ctx.clearRect(0, 0, FP_SIZE, FP_SIZE);
  ctx.fillStyle = '#13131a';
  ctx.fillRect(0, 0, FP_SIZE, FP_SIZE);

  ctx.beginPath();
  ctx.arc(cx, cy, 410 * SCALE, 0, 2 * Math.PI);
  ctx.strokeStyle = 'rgba(255,255,255,.12)'; ctx.lineWidth = 0.8; ctx.stroke();

  ctx.strokeStyle = 'rgba(255,255,255,.06)'; ctx.lineWidth = 0.6;
  for (let p = 0; p < 10; p++) {
    const ang = (p * 36) * Math.PI / 180;
    ctx.beginPath();
    ctx.moveTo(cx, cy);
    ctx.lineTo(cx + 420 * SCALE * Math.cos(ang), cy - 420 * SCALE * Math.sin(ang));
    ctx.stroke();
  }

  // background fibers
  const activeSet = new Set(activeFibers.map(String));
  ctx.fillStyle = '#2e2e3e';
  for (const [fid, pos] of Object.entries(fibPos)) {
    if (activeSet.has(fid)) continue;
    ctx.beginPath(); ctx.arc(pos[0], pos[1], 1.5, 0, 2 * Math.PI); ctx.fill();
  }

  // active fibers
  for (const fib of activeFibers) {
    const pos = fibPos[fib];
    if (!pos) continue;
    const isSel = fib === selectedFiber;
    const isHov = fib === hoverId;
    const r = isSel ? 5 : (isHov ? 4 : (badSet.has(fib) ? 3 : 2.5));
    ctx.beginPath(); ctx.arc(pos[0], pos[1], r, 0, 2 * Math.PI);
    ctx.fillStyle = fiberColor(fib, isHov, isSel);
    ctx.fill();
  }
}
drawFP(null);

function hitTest(mx, my) {
  let best = null, bestD = 144;
  for (const fib of activeFibers) {
    const pos = fibPos[fib];
    if (!pos) continue;
    const dx = pos[0] - mx, dy = pos[1] - my;
    const d2 = dx * dx + dy * dy;
    if (d2 < bestD) { bestD = d2; best = fib; }
  }
  return best;
}

let lastHov = null;
let currentView = 'timePlot';  // 'timePlot' or 'histogram'

fpCanvas.addEventListener('mousemove', e => {
  const r = fpCanvas.getBoundingClientRect();
  const fib = hitTest(e.clientX - r.left, e.clientY - r.top);
  if (fib !== lastHov) { lastHov = fib; drawFP(fib); }
  if (fib !== null) {
    const fd = DATA.fibers[fib];
    const ns  = fd.nsig !== null ? '  σ = ' + fd.nsig.toFixed(1) : '';
    const bad = fd.is_bad ? '  ★ BAD' : '';
    tooltip.textContent   = 'fiber ' + fib + '  p' + fd.petal + '   n_obs = ' + fd.n_obs + ns + bad;
    tooltip.style.left    = (e.clientX + 14) + 'px';
    tooltip.style.top     = (e.clientY - 10) + 'px';
    tooltip.style.opacity = '1';
    fpCanvas.style.cursor = 'pointer';
  } else {
    tooltip.style.opacity = '0';
    fpCanvas.style.cursor = 'crosshair';
  }
});
fpCanvas.addEventListener('mouseleave', () => {
  tooltip.style.opacity = '0'; lastHov = null; drawFP(null);
});
fpCanvas.addEventListener('click', e => {
  const r = fpCanvas.getBoundingClientRect();
  const fib = hitTest(e.clientX - r.left, e.clientY - r.top);
  if (fib !== null) { currentView = 'timePlot'; selectedFiber = fib; drawFP(lastHov); drawTimePlot(fib); }
});

// ── go-to fiber box ───────────────────────────────────────────────────────────
function goToFiber(val) {
  const inp = document.getElementById('goto-input');
  const fib = parseInt(val, 10);
  if (isNaN(fib) || !DATA.fibers[fib]) {
    inp.classList.add('error');
    fiberTitle.innerHTML = '<span style="color:#c03030">fiber ' + fib + ' not in loaded petals</span>';
    return;
  }
  inp.classList.remove('error');
  currentView = 'timePlot';
  selectedFiber = fib;
  drawFP(null);
  drawTimePlot(fib);
}
document.getElementById('goto-btn').addEventListener('click', () => {
  goToFiber(document.getElementById('goto-input').value);
});
document.getElementById('goto-input').addEventListener('keydown', e => {
  if (e.key === 'Enter') goToFiber(e.target.value);
  e.target.classList.remove('error');
});

// ── plot canvas dot mouseover ─────────────────────────────────────────────────
plotCanvas.addEventListener('mousemove', e => {
  const r  = plotCanvas.getBoundingClientRect();
  const mx = e.clientX - r.left, my = e.clientY - r.top;
  let best = null, bestD = 12 * 12;
  for (const dot of plotDots) {
    const dx = dot.x - mx, dy = dot.y - my;
    const d2 = dx * dx + dy * dy;
    if (d2 < bestD) { bestD = d2; best = dot; }
  }
  if (best && currentView === 'timePlot') {
    tooltip.textContent = mjdToDate(best.mjd_lo) + ' – ' + mjdToDate(best.mjd_hi) +
                          '   fail = ' + (best.fr * 100).toFixed(1) + '%  · click for run histogram';
    tooltip.style.left    = (e.clientX + 14) + 'px';
    tooltip.style.top     = (e.clientY - 10) + 'px';
    tooltip.style.opacity = '1';
    plotCanvas.style.cursor = 'pointer';
  } else if (currentView === 'histogram') {
    tooltip.textContent   = 'click to return to time plot';
    tooltip.style.left    = (e.clientX + 14) + 'px';
    tooltip.style.top     = (e.clientY - 10) + 'px';
    tooltip.style.opacity = '1';
    plotCanvas.style.cursor = 'pointer';
  } else {
    tooltip.style.opacity   = '0';
    plotCanvas.style.cursor = 'default';
  }
});
plotCanvas.addEventListener('mouseleave', () => { tooltip.style.opacity = '0'; });
plotCanvas.addEventListener('click', e => {
  if (selectedFiber === null) return;
  if (currentView === 'histogram') {
    currentView = 'timePlot';
    drawTimePlot(selectedFiber);
    return;
  }
  const r  = plotCanvas.getBoundingClientRect();
  const mx = e.clientX - r.left, my = e.clientY - r.top;
  let best = null, bestD = 15 * 15;
  for (const dot of plotDots) {
    const dx = dot.x - mx, dy = dot.y - my;
    if (dx * dx + dy * dy < bestD) { bestD = dx * dx + dy * dy; best = dot; }
  }
  if (best) { currentView = 'histogram'; drawClickView(selectedFiber, best); }
});

// ── time plot ─────────────────────────────────────────────────────────────────
function drawTimePlot(fib) {
  const fd   = DATA.fibers[fib];
  const wrap = document.getElementById('plot-wrap');
  const W    = wrap.clientWidth - 32;
  const H    = Math.min(Math.round(W * 0.52), 340);
  plotCanvas.width = W; plotCanvas.height = H;

  const ctx = plotCtx;
  const PAD = { l: 52, r: 18, t: 18, b: 40 };
  const pw = W - PAD.l - PAD.r, ph = H - PAD.t - PAD.b;

  ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, W, H);

  const bins  = DATA.bins_mjd;
  const sv_fr = DATA.survey_fr;
  const sv_er = DATA.survey_er;
  const fr    = fd.fr, er = fd.er, pred = fd.pred;

  const xmin = bins[0], xmax = bins[bins.length - 1];
  const xScale = v => PAD.l + (v - xmin) / (xmax - xmin) * pw;
  const allY = [...fr, ...sv_fr, ...pred].filter(v => v !== null);
  const ymax = allY.length ? Math.max(...allY) * 1.25 : 0.4;
  const yScale = v => PAD.t + ph - (v / ymax) * ph;

  const rawStep = ymax / 4;
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const step = [1, 2, 5].map(f => f * mag).find(s => s >= rawStep) || mag;
  const yticks = [];
  for (let v = step; v <= ymax * 1.05; v = Math.round((v + step) * 1e6) / 1e6) yticks.push(v);
  ctx.strokeStyle = 'rgba(0,0,0,.08)'; ctx.lineWidth = 0.6;
  for (const y of yticks) {
    ctx.beginPath(); ctx.moveTo(PAD.l, yScale(y)); ctx.lineTo(PAD.l + pw, yScale(y)); ctx.stroke();
  }

  // survey avg fill + line
  ctx.beginPath();
  let first = true;
  for (let i = 0; i < bins.length; i++) {
    if (sv_fr[i] === null || sv_er[i] === null) { first = true; continue; }
    first ? (ctx.moveTo(xScale(bins[i]), yScale(sv_fr[i] + sv_er[i])), first = false)
          : ctx.lineTo(xScale(bins[i]), yScale(sv_fr[i] + sv_er[i]));
  }
  for (let i = bins.length - 1; i >= 0; i--) {
    if (sv_fr[i] === null || sv_er[i] === null) continue;
    ctx.lineTo(xScale(bins[i]), yScale(sv_fr[i] - sv_er[i]));
  }
  ctx.closePath(); ctx.fillStyle = 'rgba(180,180,210,.14)'; ctx.fill();
  ctx.beginPath(); first = true;
  for (let i = 0; i < bins.length; i++) {
    if (sv_fr[i] === null) { first = true; continue; }
    first ? (ctx.moveTo(xScale(bins[i]), yScale(sv_fr[i])), first = false)
          : ctx.lineTo(xScale(bins[i]), yScale(sv_fr[i]));
  }
  ctx.strokeStyle = 'rgba(180,180,210,.4)'; ctx.lineWidth = 1; ctx.stroke();

  // model prediction
  ctx.beginPath(); first = true; ctx.setLineDash([4, 3]);
  for (let i = 0; i < bins.length; i++) {
    if (pred[i] === null) { first = true; continue; }
    first ? (ctx.moveTo(xScale(bins[i]), yScale(pred[i])), first = false)
          : ctx.lineTo(xScale(bins[i]), yScale(pred[i]));
  }
  ctx.strokeStyle = 'rgba(200,200,200,.6)'; ctx.lineWidth = 1.2; ctx.stroke();
  ctx.setLineDash([]);

  // year lines
  ctx.strokeStyle = 'rgba(0,0,0,.10)'; ctx.lineWidth = 0.6;
  ctx.fillStyle   = 'rgba(0,0,0,.35)';
  ctx.font = '10px ui-monospace,Menlo,monospace'; ctx.textAlign = 'left';
  for (let yr = 2021; yr <= 2026; yr++) {
    const mjd = Math.round((Date.UTC(yr, 0, 1) - Date.UTC(1858, 10, 17)) / 86400000);
    if (mjd < xmin || mjd > xmax) continue;
    const x = xScale(mjd);
    ctx.beginPath(); ctx.moveTo(x, PAD.t); ctx.lineTo(x, PAD.t + ph); ctx.stroke();
    ctx.fillText(yr, x + 2, PAD.t + ph - 3);
  }

  // CCD swap lines
  ctx.setLineDash([4, 3]); ctx.lineWidth = 1.2;
  for (const [mjd, ch] of (CCD_SWAPS[String(fd.petal)] || [])) {
    if (mjd < xmin || mjd > xmax) continue;
    const x = xScale(mjd);
    ctx.strokeStyle = CCD_COLOR[ch] + 'bb';
    ctx.beginPath(); ctx.moveTo(x, PAD.t); ctx.lineTo(x, PAD.t + ph); ctx.stroke();
  }
  ctx.setLineDash([]);

  // time plot always red
  const fcolor = '#c03030';
  document.getElementById('fb-swatch').style.background = fcolor;

  // fail rate line
  ctx.beginPath(); first = true;
  for (let i = 0; i < bins.length; i++) {
    if (fr[i] === null) { first = true; continue; }
    first ? (ctx.moveTo(xScale(bins[i]), yScale(fr[i])), first = false)
          : ctx.lineTo(xScale(bins[i]), yScale(fr[i]));
  }
  ctx.strokeStyle = fcolor + '99'; ctx.lineWidth = 1; ctx.stroke();

  // error bars + dots
  const binW = bins.length > 1 ? (bins[1] - bins[0]) : 0;
  plotDots = [];
  for (let i = 0; i < bins.length; i++) {
    if (fr[i] === null) continue;
    const x = xScale(bins[i]), y = yScale(fr[i]);
    const ey = er[i] !== null ? (er[i] / ymax) * ph : 0;
    ctx.strokeStyle = fcolor + '88'; ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(x, y - ey); ctx.lineTo(x, y + ey);
    ctx.moveTo(x - 2, y - ey); ctx.lineTo(x + 2, y - ey);
    ctx.moveTo(x - 2, y + ey); ctx.lineTo(x + 2, y + ey);
    ctx.stroke();
    ctx.beginPath(); ctx.arc(x, y, 2.8, 0, 2 * Math.PI);
    ctx.fillStyle = fcolor; ctx.fill();
    plotDots.push({ x, y, mjd_lo: bins[i] - binW/2, mjd_hi: bins[i] + binW/2, fr: fr[i] });
  }

  // axes
  ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(PAD.l, PAD.t); ctx.lineTo(PAD.l, PAD.t + ph); ctx.lineTo(PAD.l + pw, PAD.t + ph);
  ctx.stroke();
  ctx.fillStyle = 'rgba(0,0,0,.7)'; ctx.textAlign = 'right';
  for (const y of yticks) {
    const py = yScale(y);
    ctx.beginPath(); ctx.moveTo(PAD.l - 3, py); ctx.lineTo(PAD.l, py);
    ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.lineWidth = 1; ctx.stroke();
    ctx.fillText(y.toFixed(2), PAD.l - 6, py + 3.5);
  }
  ctx.textAlign = 'center'; ctx.fillStyle = 'rgba(0,0,0,.4)';
  ctx.fillText('MJD', PAD.l + pw / 2, H - 8);

  const ns        = fd.nsig !== null ? '   σ = ' + fd.nsig.toFixed(1) : '';
  const bad       = fd.is_bad ? '   ★ BAD' : '';
  const titleColor = fd.is_bad ? '#e8910a' : '#e0e0f0';
  fiberTitle.innerHTML =
    '<span style="color:' + titleColor + '">fiber ' + fib +
    ' &nbsp;·&nbsp; petal ' + fd.petal +
    ' &nbsp;·&nbsp; n_obs = ' + fd.n_obs + ns + bad + '</span>';
}

function drawClickView(fib, dot) {
  const fd    = DATA.fibers[fib];
  const wrap  = document.getElementById('plot-wrap');
  const CW    = wrap.clientWidth - 32;
  const failM = fd.obs_fail_m || [];
  const failR = fd.obs_fail_r || [];
  const rl    = fd.rl || [];

  const totalRuns = rl.reduce((a, b) => a + b, 0);
  let meanRL = 0;
  for (let i = 0; i < rl.length; i++) meanRL += (i + 1) * rl[i];
  if (totalRuns > 0) meanRL /= totalRuns;
  const maxRL  = failR.length > 0 ? Math.max(...failR) : 3;
  const yRLmax = maxRL + 0.8;

  // Canvas layout
  const TH  = Math.min(Math.round(CW * 0.33), 220);   // top panel height
  const BH  = Math.min(Math.round(CW * 0.38), 260);   // bottom panels height
  const GY  = 26;                                       // gap between rows
  const CH  = TH + GY + BH;
  plotCanvas.width = CW; plotCanvas.height = CH;

  const ctx  = plotCtx;
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, CW, CH);

  const bins = DATA.bins_mjd;
  const xmin = bins[0], xmax = bins[bins.length - 1];

  // ── helper: draw RLE lines ────────────────────────────────────────────────
  function drawRLElines(xS, yS, x0, y0, w, h, tlo, thi, highlightDot) {
    ctx.save(); ctx.beginPath();
    ctx.rect(x0, y0, w, h + 1); ctx.clip();
    for (let i = 0; i < failM.length; i++) {
      if (failM[i] < tlo || failM[i] > thi) continue;
      const x  = xS(failM[i]);
      const inW = highlightDot && failM[i] >= highlightDot.mjd_lo && failM[i] < highlightDot.mjd_hi;
      ctx.lineWidth   = inW ? 1.3 : 0.8;
      ctx.strokeStyle = inW ? 'rgba(139,26,26,.90)' : 'rgba(139,26,26,.55)';
      ctx.beginPath(); ctx.moveTo(x, yS(0)); ctx.lineTo(x, yS(failR[i])); ctx.stroke();
    }
    ctx.restore();
  }

  // ── helper: draw integer y-ticks ─────────────────────────────────────────
  function drawIntYticks(ctx, x0, yS, maxV) {
    ctx.fillStyle = 'rgba(0,0,0,.7)'; ctx.textAlign = 'right';
    ctx.font = '9px ui-monospace,Menlo,monospace';
    for (let v = 1; v <= maxV; v++) {
      const py = yS(v);
      ctx.beginPath(); ctx.moveTo(x0 - 3, py); ctx.lineTo(x0, py);
      ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.lineWidth = 1; ctx.stroke();
      ctx.fillText(v, x0 - 5, py + 3.5);
    }
  }

  // ── TOP PANEL: full-survey RLE timeline ──────────────────────────────────
  const T  = { x0: 52, y0: 22, w: CW - 70, h: TH - 57 };
  const txS = v => T.x0 + (v - xmin) / (xmax - xmin) * T.w;
  const tyS = v => T.y0 + T.h - (v / yRLmax) * T.h;

  // gold highlight
  if (dot) {
    const hx0 = txS(dot.mjd_lo), hx1 = txS(dot.mjd_hi);
    ctx.fillStyle = 'rgba(255,200,50,.18)';
    ctx.fillRect(hx0, T.y0, hx1 - hx0, T.h);
  }

  // year grid
  ctx.strokeStyle = 'rgba(0,0,0,.10)'; ctx.lineWidth = 0.5;
  ctx.fillStyle   = 'rgba(0,0,0,.30)';
  ctx.font = '9px ui-monospace,Menlo,monospace'; ctx.textAlign = 'left';
  for (let yr = 2021; yr <= 2026; yr++) {
    const mjd = Math.round((Date.UTC(yr, 0, 1) - Date.UTC(1858, 10, 17)) / 86400000);
    if (mjd < xmin || mjd > xmax) continue;
    const x = txS(mjd);
    ctx.beginPath(); ctx.moveTo(x, T.y0); ctx.lineTo(x, T.y0 + T.h); ctx.stroke();
    ctx.fillText(yr, x + 2, T.y0 + T.h - 2);
  }

  // CCD swaps
  ctx.setLineDash([4, 3]); ctx.lineWidth = 1;
  for (const [mjd, ch] of (CCD_SWAPS[String(fd.petal)] || [])) {
    if (mjd < xmin || mjd > xmax) continue;
    ctx.strokeStyle = CCD_COLOR[ch] + 'bb';
    const x = txS(mjd);
    ctx.beginPath(); ctx.moveTo(x, T.y0); ctx.lineTo(x, T.y0 + T.h); ctx.stroke();
  }
  ctx.setLineDash([]);

  // mean line
  if (meanRL > 0) {
    const ym = tyS(meanRL);
    ctx.setLineDash([4, 3]); ctx.lineWidth = 0.8; ctx.strokeStyle = 'rgba(0,0,0,.28)';
    ctx.beginPath(); ctx.moveTo(T.x0, ym); ctx.lineTo(T.x0 + T.w, ym); ctx.stroke();
    ctx.setLineDash([]);
  }

  drawRLElines(txS, tyS, T.x0, T.y0, T.w, T.h, xmin, xmax, dot);

  // axes
  ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(T.x0, T.y0); ctx.lineTo(T.x0, T.y0 + T.h); ctx.lineTo(T.x0 + T.w, T.y0 + T.h);
  ctx.stroke();
  drawIntYticks(ctx, T.x0, tyS, maxRL);
  ctx.textAlign = 'center'; ctx.fillStyle = 'rgba(0,0,0,.4)'; ctx.font = '10px ui-monospace,Menlo,monospace';
  ctx.fillText('MJD', T.x0 + T.w / 2, TH - 5);
  ctx.save(); ctx.translate(10, T.y0 + T.h / 2); ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center'; ctx.fillText('Run length', 0, 0); ctx.restore();
  ctx.fillStyle = 'rgba(0,0,0,.5)'; ctx.textAlign = 'center';
  ctx.fillText('Run-length encoding — full survey  (' + fd.n_obs + ' obs, ' + totalRuns + ' failure runs)', T.x0 + T.w / 2, T.y0 - 6);

  // ── BOTTOM-LEFT PANEL: zoomed RLE ────────────────────────────────────────
  const BLW = Math.round(CW * 0.56);
  const BY0 = TH + GY;
  const BL  = { x0: 52, y0: BY0 + 22, w: BLW - 60, h: BH - 57 };

  let tlo = xmin, thi = xmax;
  if (dot) {
    const span = dot.mjd_hi - dot.mjd_lo;
    tlo = dot.mjd_lo - span * 0.8;
    thi = dot.mjd_hi + span * 0.8;
  }
  const blxS = v => BL.x0 + (v - tlo) / (thi - tlo) * BL.w;
  const blyS = v => BL.y0 + BL.h - (v / yRLmax) * BL.h;

  // gold highlight
  if (dot) {
    ctx.fillStyle = 'rgba(255,200,50,.18)';
    ctx.fillRect(blxS(dot.mjd_lo), BL.y0, blxS(dot.mjd_hi) - blxS(dot.mjd_lo), BL.h);
  }

  // month labels
  const MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
  ctx.strokeStyle = 'rgba(0,0,0,.10)'; ctx.lineWidth = 0.5;
  ctx.fillStyle   = 'rgba(0,0,0,.30)'; ctx.font = '8px ui-monospace,Menlo,monospace';
  for (let yr = 2020; yr <= 2027; yr++) {
    for (let mo = 0; mo < 12; mo++) {
      const mjd = Math.round((Date.UTC(yr, mo, 1) - Date.UTC(1858, 10, 17)) / 86400000);
      if (mjd < tlo || mjd > thi) continue;
      const x = blxS(mjd);
      if (x < BL.x0 || x > BL.x0 + BL.w) continue;
      ctx.beginPath(); ctx.moveTo(x, BL.y0); ctx.lineTo(x, BL.y0 + BL.h); ctx.stroke();
      ctx.textAlign = 'center'; ctx.fillText(MONTHS[mo], x, BL.y0 + BL.h + 10);
    }
  }

  // CCD swaps
  ctx.setLineDash([4, 3]); ctx.lineWidth = 1;
  for (const [mjd, ch] of (CCD_SWAPS[String(fd.petal)] || [])) {
    if (mjd < tlo || mjd > thi) continue;
    ctx.strokeStyle = CCD_COLOR[ch] + 'bb';
    const x = blxS(mjd);
    ctx.beginPath(); ctx.moveTo(x, BL.y0); ctx.lineTo(x, BL.y0 + BL.h); ctx.stroke();
  }
  ctx.setLineDash([]);

  // mean line
  if (meanRL > 0) {
    const ym = blyS(meanRL);
    ctx.setLineDash([4, 3]); ctx.lineWidth = 0.8; ctx.strokeStyle = 'rgba(0,0,0,.28)';
    ctx.beginPath(); ctx.moveTo(BL.x0, ym); ctx.lineTo(BL.x0 + BL.w, ym); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(0,0,0,.4)'; ctx.font = '8px ui-monospace,Menlo,monospace'; ctx.textAlign = 'left';
    ctx.fillText('mean=' + meanRL.toFixed(1), BL.x0 + 3, ym - 3);
  }

  let nFail = 0;
  for (let i = 0; i < failM.length; i++) {
    if (dot && failM[i] >= dot.mjd_lo && failM[i] < dot.mjd_hi) nFail++;
  }
  drawRLElines(blxS, blyS, BL.x0, BL.y0, BL.w, BL.h, tlo, thi, dot);

  // axes
  ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(BL.x0, BL.y0); ctx.lineTo(BL.x0, BL.y0 + BL.h); ctx.lineTo(BL.x0 + BL.w, BL.y0 + BL.h);
  ctx.stroke();
  drawIntYticks(ctx, BL.x0, blyS, maxRL);
  ctx.textAlign = 'center'; ctx.fillStyle = 'rgba(0,0,0,.4)'; ctx.font = '10px ui-monospace,Menlo,monospace';
  ctx.fillText('MJD', BL.x0 + BL.w / 2, BY0 + BH - 5);
  ctx.save(); ctx.translate(10, BL.y0 + BL.h / 2); ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center'; ctx.fillText('Run length', 0, 0); ctx.restore();
  const zTitle = dot ? mjdToDate(dot.mjd_lo) + ' – ' + mjdToDate(dot.mjd_hi) + '   (' + nFail + ' fails in window)'
                     : 'Full survey';
  ctx.fillStyle = 'rgba(0,0,0,.5)'; ctx.textAlign = 'center'; ctx.font = '10px ui-monospace,Menlo,monospace';
  ctx.fillText(zTitle, BL.x0 + BL.w / 2, BY0 + 12);

  // ── BOTTOM-RIGHT PANEL: histogram (for clicked window) ─────────────────────
  const BRX0 = BLW + 20;
  const BR   = { x0: BRX0 + 36, y0: BY0 + 22, w: CW - BRX0 - 48, h: BH - 57 };

  // Compute run-length counts for the window using run IDs to avoid double-counting
  const failRid = fd.obs_fail_rid || [];
  const winCounts = {};
  const seenRids  = new Set();
  for (let i = 0; i < failM.length; i++) {
    if (!dot || (failM[i] >= dot.mjd_lo && failM[i] < dot.mjd_hi)) {
      if (!seenRids.has(failRid[i])) {
        seenRids.add(failRid[i]);
        winCounts[failR[i]] = (winCounts[failR[i]] || 0) + 1;
      }
    }
  }
  const winMaxLen   = Object.keys(winCounts).length > 0 ? Math.max(...Object.keys(winCounts).map(Number)) : 0;
  const winRL       = winMaxLen > 0 ? Array.from({length: winMaxLen}, (_, i) => winCounts[i+1] || 0) : [];
  const winTotal    = winRL.reduce((a, b) => a + b, 0);

  if (winRL.length > 0) {
    const maxCount = Math.max(...winRL);
    const yHmax = maxCount * 1.25;
    const brxS  = i => BR.x0 + (i / winRL.length) * BR.w;
    const bryS  = v => BR.y0 + BR.h - (v / yHmax) * BR.h;
    const barW  = BR.w / winRL.length * 0.72;

    const rawStep3 = Math.max(yHmax / 4, 1);
    const mag3  = Math.pow(10, Math.floor(Math.log10(rawStep3)));
    const step3 = [1, 2, 5].map(f => f * mag3).find(s => s >= rawStep3) || 1;
    const yticks3 = [];
    for (let v = step3; v <= yHmax * 1.05; v = Math.round((v + step3) * 1e6) / 1e6) yticks3.push(v);

    ctx.strokeStyle = 'rgba(0,0,0,.08)'; ctx.lineWidth = 0.6;
    for (const y of yticks3) {
      ctx.beginPath(); ctx.moveTo(BR.x0, bryS(y)); ctx.lineTo(BR.x0 + BR.w, bryS(y)); ctx.stroke();
    }

    for (let i = 0; i < winRL.length; i++) {
      if (winRL[i] === 0) continue;
      const cx = brxS(i + 0.5);
      ctx.fillStyle = '#8b1a1a';
      ctx.fillRect(cx - barW / 2, bryS(winRL[i]), barW, bryS(0) - bryS(winRL[i]));
      ctx.fillStyle = 'rgba(0,0,0,.7)'; ctx.font = '10px ui-monospace,Menlo,monospace'; ctx.textAlign = 'center';
      ctx.fillText(winRL[i], cx, bryS(winRL[i]) - 4);
    }

    // mean line
    let meanH = 0;
    for (let i = 0; i < winRL.length; i++) meanH += (i + 1) * winRL[i];
    if (winTotal > 0) meanH /= winTotal;
    const mxH = brxS(meanH);
    ctx.setLineDash([4, 3]); ctx.lineWidth = 1; ctx.strokeStyle = 'rgba(0,0,0,.35)';
    ctx.beginPath(); ctx.moveTo(mxH, BR.y0); ctx.lineTo(mxH, BR.y0 + BR.h); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(0,0,0,.4)'; ctx.font = '8px ui-monospace,Menlo,monospace'; ctx.textAlign = 'left';
    ctx.fillText('μ=' + meanH.toFixed(1), mxH + 2, BR.y0 + 10);

    // axes
    ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(BR.x0, BR.y0); ctx.lineTo(BR.x0, BR.y0 + BR.h); ctx.lineTo(BR.x0 + BR.w, BR.y0 + BR.h);
    ctx.stroke();
    ctx.fillStyle = 'rgba(0,0,0,.7)'; ctx.textAlign = 'right'; ctx.font = '9px ui-monospace,Menlo,monospace';
    for (const y of yticks3) {
      const py = bryS(y);
      ctx.beginPath(); ctx.moveTo(BR.x0 - 3, py); ctx.lineTo(BR.x0, py);
      ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.stroke();
      ctx.fillText(y, BR.x0 - 5, py + 3.5);
    }
    ctx.textAlign = 'center'; ctx.fillStyle = 'rgba(0,0,0,.4)'; ctx.font = '10px ui-monospace,Menlo,monospace';
    for (let i = 0; i < winRL.length; i++) {
      const cx = brxS(i + 0.5);
      ctx.beginPath(); ctx.moveTo(cx, BR.y0 + BR.h); ctx.lineTo(cx, BR.y0 + BR.h + 4);
      ctx.strokeStyle = 'rgba(0,0,0,.4)'; ctx.stroke();
      ctx.fillStyle = 'rgba(0,0,0,.7)'; ctx.fillText(i + 1, cx, BR.y0 + BR.h + 14);
    }
    ctx.fillStyle = 'rgba(0,0,0,.4)';
    ctx.fillText('Run length', BR.x0 + BR.w / 2, BY0 + BH - 5);
    ctx.fillStyle = 'rgba(0,0,0,.5)'; ctx.textAlign = 'center';
    ctx.fillText('Run-length dist — window  (' + winTotal + ' runs)', BR.x0 + BR.w / 2, BY0 + 12);
  } else {
    ctx.fillStyle = 'rgba(0,0,0,.4)'; ctx.font = '11px ui-monospace,Menlo,monospace';
    ctx.textAlign = 'center'; ctx.fillText('No failure runs in window', BR.x0 + BR.w / 2, BR.y0 + BR.h / 2);
  }

  // ── Fiber title ───────────────────────────────────────────────────────────
  const ns    = fd.nsig !== null ? '   σ = ' + fd.nsig.toFixed(1) : '';
  const bad   = fd.is_bad ? '   ★ BAD' : '';
  const tcolor = fd.is_bad ? '#e8910a' : '#e0e0f0';
  fiberTitle.innerHTML =
    '<span style="color:' + tcolor + '">fiber ' + fib +
    ' &nbsp;·&nbsp; petal ' + fd.petal + ns + bad +
    ' &nbsp;·&nbsp; click to return</span>';
}

const overlay     = document.getElementById('img-overlay');
const overlayImg  = document.getElementById('img-overlay-img');
document.getElementById('img-overlay-close').addEventListener('click', () => overlay.classList.remove('show'));
overlay.addEventListener('click', e => { if (e.target === overlay) overlay.classList.remove('show'); });

function showCanvasOverlay(canvas) {
  overlayImg.src = canvas.toDataURL('image/png');
  overlay.classList.add('show');
}

document.getElementById('save-plot-btn').addEventListener('click', () => showCanvasOverlay(plotCanvas));
</script>
"""

colors_js    = json.dumps([PETAL_COLORS[i % len(PETAL_COLORS)] for i in range(len(PETALS))])
ccd_swaps_js = json.dumps({str(k): v for k, v in CCD_SWAPS_ALL.items()})
html = HTML.replace('__JSON__', raw_json) \
           .replace('__LABEL__', label) \
           .replace('__COLORS__', colors_js) \
           .replace('__CCD_SWAPS__', ccd_swaps_js)

with open(out_path, 'w') as f:
    f.write(html)
print(f'Saved {out_path}  ({os.path.getsize(out_path)//1024} KB)', flush=True)
