// public/script.js

// Chart instances
let zChart = null, tChart = null, chi2Chart = null, fChart = null;

// Utility: mean & variance for step breakdown
function mean(x) {
  return x.reduce((s, v) => s + v, 0) / x.length;
}
function variance(x, ddof = 0) {
  const μ = mean(x);
  return x.reduce((s, v) => s + (v - μ) ** 2, 0) / (x.length - ddof);
}

// PDF helpers (Gamma via Lanczos, plus t, χ², and F PDFs)
function gamma(z) {
  const g = 7;
  const p = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
  ];
  if (z < 0.5) {
    return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z));
  }
  z -= 1;
  let x = p[0];
  for (let i = 1; i < p.length; i++) {
    x += p[i] / (z + i);
  }
  const t = z + g + 0.5;
  return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
}

function tPdf(t, v) {
  return (gamma((v + 1) / 2) /
    (Math.sqrt(v * Math.PI) * gamma(v / 2))) *
    Math.pow(1 + (t * t) / v, -(v + 1) / 2);
}

function chi2Pdf(x, k) {
  return Math.pow(x, k / 2 - 1) * Math.exp(-x / 2) /
    (Math.pow(2, k / 2) * gamma(k / 2));
}

function fPdf(x, d1, d2) {
  const num = Math.pow(d1 * x, d1) * Math.pow(d2, d2);
  const den = Math.pow(d1 * x + d2, d1 + d2);
  return Math.sqrt(num / den) /
    (x * (gamma(d1 / 2) * gamma(d2 / 2) / gamma((d1 + d2) / 2)));
}

// ── Z-Test ──────────────────────────────────────────────────────────────
async function runZ() {
  const data = document.getElementById('z-data').value.split(',').map(Number);
  const mu0 = parseFloat(document.getElementById('z-mu0').value);
  const sigma = parseFloat(document.getElementById('z-sigma').value);
  const alt = document.getElementById('z-alt').value;
  const res = await fetch('/api/ztest', {
    method: 'POST',
    headers: {'Content-Type':'application/json'},
    body: JSON.stringify({ data, mu0, sigma, alternative: alt })
  }).then(r=>r.json());

  document.getElementById('z-result')
    .textContent = `Z = ${res.statistic.toFixed(3)}, p = ${res.pValue.toFixed(3)}`;

  const n = data.length;
  const xbar = data.reduce((s,v)=>s+v,0)/n;
  const se = sigma/Math.sqrt(n);
  const zVal = (xbar - mu0)/se;
  const steps = [
    `n = ${n}`,
    `x̄ = (${data.join(' + ')})/ ${n} = ${xbar.toFixed(3)}`,
    `σ = ${sigma}`,
    `SE = σ/√n = ${se.toFixed(3)}`,
    ``,
    `Z = (x̄−μ₀)/SE = (${xbar.toFixed(3)}−${mu0})/${se.toFixed(3)} = ${zVal.toFixed(3)}`,
    ``,
    `Alt: ${alt}`,
    alt==='two-sided'
      ? `p = 2·(1−Φ(|${zVal.toFixed(3)}|)) = ${res.pValue.toFixed(3)}`
      : alt==='greater'
        ? `p = 1−Φ(${zVal.toFixed(3)}) = ${res.pValue.toFixed(3)}`
        : `p = Φ(${zVal.toFixed(3)}) = ${res.pValue.toFixed(3)}`
  ].join('\n');
  document.getElementById('z-steps-content').textContent = steps;

  drawZChart(zVal, alt);
}

function drawZChart(z, alternative) {
  const N = 500;
  const xs = Array.from({ length: N }, (_, i) => -4 + 8 * i/(N-1));
  const ys = xs.map(x => (1/Math.sqrt(2*Math.PI))*Math.exp(-0.5*x*x));
  const crit = 1.96;
  const bg = xs.map(x => {
    if(alternative==='two-sided'&&(x<=-crit||x>=crit)) return 'rgba(255,99,132,0.4)';
    if(alternative==='greater'&& x>=crit)           return 'rgba(255,99,132,0.4)';
    if(alternative==='less'   && x<=-crit)           return 'rgba(255,99,132,0.4)';
    return 'rgba(0,0,0,0)';
  });
  if (zChart) zChart.destroy();
  const ctx = document.getElementById('z-chart').getContext('2d');
  zChart = new Chart(ctx, {
    type: 'line',
    data: {
      labels: xs,
      datasets: [
        { data: ys, borderColor:'blue', borderWidth:2, fill:false, pointRadius:0 },
        { data: ys.map((y,i)=>({x:xs[i],y})), backgroundColor:bg, type:'bar',
          barPercentage:1, categoryPercentage:1 },
        { data:[{x:z,y:0},{x:z,y:(1/Math.sqrt(2*Math.PI))*Math.exp(-0.5*z*z)}],
          borderColor:'black', borderWidth:2, fill:false, showLine:true, pointRadius:0 }
      ]
    },
    options: {
      scales: {
        x:{ type:'linear', title:{ display:true, text:'z' } },
        y:{ title:{ display:true, text:'Density' }}
      },
      plugins:{ legend:{ display:false }}
    }
  });
}

// ── T-Test ──────────────────────────────────────────────────────────────
async function runT() {
  const data = document.getElementById('t-data').value.split(',').map(Number);
  const mu0 = parseFloat(document.getElementById('t-mu0').value);
  const alt = document.getElementById('t-alt').value;
  const res = await fetch('/api/ttest', {
    method:'POST',
    headers:{'Content-Type':'application/json'},
    body:JSON.stringify({ data, mu0, alternative: alt })
  }).then(r=>r.json());

  document.getElementById('t-result')
    .textContent = `t = ${res.statistic.toFixed(3)}, p = ${res.pValue.toFixed(3)}`;

  const n = data.length;
  const meanVal = mean(data);
  const s = Math.sqrt(data.reduce((S,v)=>S+(v-meanVal)**2,0)/(n-1));
  const se = s/Math.sqrt(n);
  const tVal = (meanVal - mu0)/se;
  const steps = [
    `n = ${n}`,
    `x̄ = (${data.join(' + ')})/ ${n} = ${meanVal.toFixed(3)}`,
    `s = √[Σ(xᵢ−x̄)²/(n−1)] = ${s.toFixed(3)}`,
    `SE = s/√n = ${se.toFixed(3)}`,
    ``,
    `t = (x̄−μ₀)/SE = (${meanVal.toFixed(3)}−${mu0})/${se.toFixed(3)} = ${tVal.toFixed(3)}`,
    ``,
    `Alt: ${alt}`,
    alt==='two-sided'
      ? `p = 2·(1−Tₙ₋₁(|${tVal.toFixed(3)}|)) = ${res.pValue.toFixed(3)}`
      : alt==='greater'
        ? `p = 1−Tₙ₋₁(${tVal.toFixed(3)}) = ${res.pValue.toFixed(3)}`
        : `p = Tₙ₋₁(${tVal.toFixed(3)}) = ${res.pValue.toFixed(3)}`
  ].join('\n');
  document.getElementById('t-steps-content').textContent = steps;

  drawTChart(tVal, alt, n-1);
}

function drawTChart(tVal, alt, df) {
  const N = 500, min=-5, max=5;
  const xs = Array.from({length:N},(_,i)=>min + (max-min)*i/(N-1));
  const ys = xs.map(x=>tPdf(x,df));
  if (tChart) tChart.destroy();
  const ctx = document.getElementById('t-chart').getContext('2d');
  tChart = new Chart(ctx, {
    type:'line',
    data:{ labels: xs, datasets:[{ data:ys, borderColor:'green', borderWidth:2, fill:false }]},
    options:{
      scales:{ x:{ title:{ display:true, text:'t' }}, y:{ title:{ display:true, text:'Density' }}},
      plugins:{ legend:{ display:false }}
    }
  });
}

// ── χ²-Test ─────────────────────────────────────────────────────────────
async function runChi2() {
  const obs = document.getElementById('chi2-observed').value.split(',').map(Number);
  const exp = document.getElementById('chi2-expected').value.split(',').map(Number);
  const res = await fetch('/api/chi2', {
    method:'POST',
    headers:{'Content-Type':'application/json'},
    body:JSON.stringify({ observed: obs, expected: exp })
  }).then(r=>r.json());

  document.getElementById('chi2-result')
    .textContent = `χ² = ${res.statistic.toFixed(3)}, p = ${res.pValue.toFixed(3)}`;

  const chi2 = obs.reduce((sum,o,i)=>sum + (o-exp[i])**2/exp[i], 0);
  const df = obs.length - 1;
  const steps = [
    `Observed: [${obs.join(', ')}]`,
    `Expected: [${exp.join(', ')}]`,
    ``,
    `χ² = Σ (Oᵢ−Eᵢ)² / Eᵢ = ${chi2.toFixed(3)}`,
    `df = k−1 = ${obs.length}−1 = ${df}`,
    ``,
    `p = 1−F_χ²(χ², df) = ${res.pValue.toFixed(3)}`
  ].join('\n');
  document.getElementById('chi2-steps-content').textContent = steps;

  drawChi2Chart(chi2, df);
}

function drawChi2Chart(chi2, df) {
  const N = 500, max = Math.max(chi2*1.5, df*3);
  const xs = Array.from({length:N},(_,i)=> i*max/(N-1));
  const ys = xs.map(x=>chi2Pdf(x,df));
  if (chi2Chart) chi2Chart.destroy();
  const ctx = document.getElementById('chi2-chart').getContext('2d');
  chi2Chart = new Chart(ctx, {
    type:'line',
    data:{ labels: xs, datasets:[{ data:ys, borderColor:'orange', borderWidth:2, fill:false }]},
    options:{
      scales:{ x:{ title:{ display:true, text:'χ²' }}, y:{ title:{ display:true, text:'Density' }}},
      plugins:{ legend:{ display:false }}
    }
  });
}

// ── F-Test ──────────────────────────────────────────────────────────────
async function runF() {
  const resultEl = document.getElementById('f-result');
  const stepsEl  = document.getElementById('f-steps-content');

  try {
    // 1) read inputs
    const d1  = document.getElementById('f-data1').value.split(',').map(Number);
    const d2  = document.getElementById('f-data2').value.split(',').map(Number);
    const alt = document.getElementById('f-alt').value;

    // 2) call endpoint
    const resp = await fetch('/api/ftest', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ data1: d1, data2: d2, alternative: alt })
    });

    if (!resp.ok) {
      throw new Error(`Server returned ${resp.status} ${resp.statusText}`);
    }
    const data = await resp.json();
    console.log('F-Test response JSON:', data);

    // 3) guard against missing/null values
    if (data.statistic == null || data.pValue == null) {
      throw new Error('Invalid response from server (statistic or pValue is null)');
    }

    // 4) render results
    resultEl.textContent =
      `F = ${data.statistic.toFixed(3)}, p = ${data.pValue.toFixed(3)} ` +
      `(df1=${data.df1}, df2=${data.df2})`;

    // 5) render computation steps
    const steps = [
      `n₁ = ${d1.length}, n₂ = ${d2.length}`,
      `s₁² = ${variance(d1,1).toFixed(3)}`,
      `s₂² = ${variance(d2,1).toFixed(3)}`,
      `F = s₁²/s₂² = ${data.statistic.toFixed(3)}`,
      `df1 = ${data.df1}, df2 = ${data.df2}`,
      ``,
      `Alt: ${alt}`,
      `p-value = ${data.pValue.toFixed(3)}`
    ].join('\n');
    stepsEl.textContent = steps;

    // 6) draw chart
    drawFChart(data.statistic, alt, data.df1, data.df2);

  } catch (err) {
    console.error('runF error:', err);
    resultEl.textContent = `Error: ${err.message}`;
    stepsEl.textContent = '';      // clear old steps
  }
}

// ── Wire up forms ───────────────────────────────────────────────────────
document.getElementById('z-form').addEventListener('submit', e=>{ e.preventDefault(); runZ(); });
document.getElementById('t-form').addEventListener('submit', e=>{ e.preventDefault(); runT(); });
document.getElementById('chi2-form').addEventListener('submit', e=>{ e.preventDefault(); runChi2(); });
document.getElementById('f-form').addEventListener('submit', e=>{ e.preventDefault(); runF(); });
