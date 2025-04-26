// stats.js
// --- 1) utility functions ------------------------------------------------
function mean(x) {
    return x.reduce((s,v) => s+v, 0) / x.length;
  }
  function variance(x, ddof = 0) {
    const μ = mean(x);
    return x.reduce((s,v) => s + (v-μ)**2, 0) / (x.length - ddof);
  }
  function std(x, ddof = 0) {
    return Math.sqrt(variance(x, ddof));
  }
  
  // --- 2) error function + normal CDF --------------------------------------
  function erf(x) {
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);
    // Abramowitz & Stegun approximation
    const [a1,a2,a3,a4,a5,p] = [0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429, 0.3275911];
    const t = 1/(1 + p*x);
    const y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);
    return sign*y;
  }
  function normalCdf(z) {
    return 0.5*(1 + erf(z/Math.sqrt(2)));
  }
  
  // --- 3) numeric integrator (Simpson’s rule) ------------------------------
  function integrate(f, a, b, n = 1000) {
    if (n % 2 === 1) n++;
    const h = (b - a)/n;
    let s = f(a) + f(b);
    for (let i = 1; i < n; i++) {
      s += f(a + i*h) * (i%2 === 0 ? 2 : 4);
    }
    return s * h/3;
  }
  
  // --- 4) Gamma via Lanczos approximation ---------------------------------
  function gamma(z) {
    const g = 7;
    const p = [
      0.99999999999980993,
      676.5203681218851,  -1259.1392167224028,
      771.32342877765313, -176.61502916214059,
      12.507343278686905, -0.13857109526572012,
      9.9843695780195716e-6, 1.5056327351493116e-7
    ];
    if (z < 0.5) {
      return Math.PI / (Math.sin(Math.PI*z) * gamma(1-z));
    }
    z -= 1;
    let x = p[0];
    for (let i = 1; i < g+2; i++) {
      x += p[i]/(z + i);
    }
    const t = z + g + 0.5;
    return Math.sqrt(2*Math.PI) * Math.pow(t, z+0.5) * Math.exp(-t) * x;
  }
  
  // --- 5) t‐distribution PDF & CDF ------------------------------------------
  function tPdf(t, v) {
    return ( gamma((v+1)/2) /
             ( Math.sqrt(v*Math.PI) * gamma(v/2) ) )
           * Math.pow(1 + (t*t)/v, -(v+1)/2);
  }
  function tCdf(t, v) {
    // integrate from –∞→t; in practice from –10 to t
    return integrate(u => tPdf(u, v), -10, t, 2000);
  }
  
  // --- 6) χ² PDF & CDF ------------------------------------------------------
  function chi2Pdf(x, k) {
    return Math.pow(x, k/2 - 1) * Math.exp(-x/2)
           / ( Math.pow(2, k/2) * gamma(k/2) );
  }
  function chi2Cdf(x, k) {
    return integrate(u => chi2Pdf(u, k), 0, x, 2000);
  }
  
  // --- 7) test functions ----------------------------------------------------
  function zTest1Sample(x, mu0, sigma, alt='two-sided') {
    const n = x.length, xbar = mean(x);
    const z = (xbar - mu0)/(sigma/Math.sqrt(n));
    let p;
    if (alt==='greater')  p = 1 - normalCdf(z);
    else if (alt==='less') p = normalCdf(z);
    else                  p = 2*(1 - normalCdf(Math.abs(z)));
    return { statistic: z, pValue: p };
  }
  
  function tTest1Sample(x, mu0, alt='two-sided') {
    const n = x.length, ν = n-1, xbar = mean(x), s = std(x,1);
    const t = (xbar - mu0)/(s/Math.sqrt(n));
    let p;
    if (alt==='greater')  p = 1 - tCdf(t,ν);
    else if (alt==='less') p = tCdf(t,ν);
    else                  p = 2*(1 - tCdf(Math.abs(t),ν));
    return { statistic: t, pValue: p };
  }
  
  function chi2Gof(obs, exp) {
    const χ2 = obs.reduce((s,o,i) => s + (o-exp[i])**2/exp[i], 0);
    const df = obs.length - 1;
    return { statistic: χ2, pValue: 1 - chi2Cdf(χ2, df) };
  }

  // ── Beta function ────────────────────────────────────────────────────────
function beta(a, b) {
    return gamma(a)*gamma(b)/gamma(a+b);
  }
  
  // ── F-distribution PDF & CDF ─────────────────────────────────────────────
  function fPdf(x, d1, d2) {
    const num = Math.pow(d1*x, d1)*Math.pow(d2, d2);
    const den = Math.pow(d1*x + d2, d1 + d2);
    return Math.sqrt(num/den)/(x * beta(d1/2, d2/2));
  }
  function fCdf(x, d1, d2) {
    return integrate(u => fPdf(u, d1, d2), 0, x, 2000);
  }
  
  // ── F-test (two-sample variance) ──────────────────────────────────────────
  function fTestTwoSample(x, y, alt='two-sided') {
    const n1 = x.length, n2 = y.length;
    const v1 = variance(x,1), v2 = variance(y,1);
    const F = v1/v2, df1 = n1-1, df2 = n2-1;
    let p;
    if (alt==='greater')  p = 1 - fCdf(F, df1, df2);
    else if (alt==='less') p = fCdf(F, df1, df2);
    else                  p = 2*Math.min(fCdf(F, df1, df2), 1 - fCdf(F, df1, df2));
    return { statistic: F, pValue: p, df1, df2 };
  }
  
  module.exports = {
    zTest1Sample,
    tTest1Sample,
    chi2Gof,
    fTestTwoSample
  };
  