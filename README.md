# 786-M14 Omnibus Extended — multi-paradigm evidence synthesis dashboard

A single-file, **fully offline** browser dashboard that runs several
meta-analytic paradigms over one long-format dataset: pairwise pooling
(including 3-level multilevel and robust-variance estimation), network
meta-analysis (fixed / DL / REML / multilevel), component NMA, bivariate
diagnostic-test-accuracy, dose-response trend, a GOSH plot, BIC model-averaging
(RoBMA-style), a Metropolis–Hastings MCMC sampler, and Risk-of-Bias / GRADE
panels.

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN — Plotly is vendored locally.

## Layout

```
index.html     single-file UI (loads plotly.min.js then engine.js)
engine.js      pure statistical core — runs unchanged in Node and the browser
tests.js       Node test harness, 42 assertions, hand-derived expectations
plotly.min.js  vendored Plotly 2.24.1 (offline)
LICENSE        Apache License 2.0 (as originally shipped in this repo)
```

## Statistical core (`engine.js`)

All pure functions/objects were extracted verbatim from the inline app script so
there is a single source of truth; the page now loads `engine.js` and the inline
duplicates were deleted.

| Object | What it does |
|---|---|
| `Stat` | `pnorm` (normal CDF via `½(1+erf(x/√2))`), `qnorm` (Acklam inverse), `qt` (Cornish–Fisher t-quantile), `invLogit` |
| `Matrix` | dense `zeros` / `t` (transpose) / `dot` / `inv` (Gauss–Jordan with partial pivoting) |
| `Optim.nelderMead` | derivative-free simplex minimiser (used by the DTA bivariate fit) |
| `Pooling.pool` | log-scale pairwise pooling: DL τ²=`max(0,(Q−(k−1))/C)`, `C=ΣW−ΣW²/ΣW`; I²; 95% CI; HKSJ option with floor `√max(1,Q/df)` and `t_{k−1}` critical; prediction interval |
| `Pooling.egger` | Egger's small-study precision-weighted regression test |
| `Multilevel.fit` | 3-level (within/between cluster) random-effects IGLS fit + Cheung I² partition |
| `RVE.fit` | cluster-robust (sandwich) variance estimator, `t_{m−1}` critical |
| `NMA_RE.solve` | network meta-analysis via GLS on the contrast design (fixed / DL / REML) |
| `NMA_Multilevel.solve` | multilevel NMA with separate within/between τ² |

Effect sizes for ratio metrics are pooled on the **log scale** and
back-transformed only for display.

## Fixes applied during revival (2026-06-05)

- **Offline**: vendored Plotly 2.24.1 locally (`plotly.min.js`); switched the
  `<script>` from the `cdn.plot.ly` URL to the local file; removed the Google
  Fonts `<link>`/preconnect (system fonts fall back). The page now loads **zero
  external resources**.
- **Single source of truth**: extracted the pure statistical core into
  `engine.js` and deleted the inline duplicates; `index.html` loads `engine.js`.
- **Added `tests.js`** (42 assertions, all passing) with hand-derived expected
  values.
- **Fixed a `k=1` / zero-heterogeneity NaN bug**: for a single study (or two
  identical studies) the DL denominator `C = ΣW − ΣW²/ΣW` is 0 and `Q = 0`, so
  `τ² = (Q−df)/C = 0/0 = NaN` poisoned the pooled estimate, and `I² = (Q−df)/Q`
  was likewise `NaN`. Both now degrade to the correct fixed-effect result
  (`τ²=0`, `I²=0`). All other math was verified correct and left unchanged —
  including the normal CDF, which already correctly returns 0.5 at x=0.
- Renamed `Omni.html` → `index.html`; added `.nojekyll`, `.gitignore`, this
  README and an E156 protocol.

## Tests

```
node tests.js
# 42 passed, 0 failed
```

Checks include normal-CDF reference points (Φ(0)=0.5 not 0; Φ(1.96)=0.975), a
hand-worked DerSimonian–Laird example (τ²≈0.05333, Q≈4.6667, I²≈57.14%,
est≈0.06667, se≈0.17638), the HKSJ floor widening behaviour, single-study and
two-identical edge cases, an Egger empty/`n<3` guard, the matrix identity
`A·A⁻¹=I`, and NMA / RVE / Multilevel property checks.

## Caveats

DerSimonian–Laird under-estimates τ² for small *k* (REML/Paule–Mandel preferred
for k<10); the REML option in the NMA module is a light fixed-point iteration,
not a full profile-likelihood solver. The MCMC and GOSH panels are stochastic.
The DTA AUC is a heuristic summary, not an HSROC integral. Treat all pooled
outputs as exploratory synthesis aids, not as a substitute for a fully
specified analysis plan. Licensed under Apache 2.0 as originally shipped.
