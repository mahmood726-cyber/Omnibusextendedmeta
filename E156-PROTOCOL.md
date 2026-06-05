# E156-PROTOCOL — 786-M14 Omnibus Extended

- **Project:** Omnibusextendedmeta (GitHub repo `Omnibusextendedmeta`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `Omni.html` dump, title "786-M14 Omnibus Extended")
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: vendored Plotly 2.24.1 locally (`plotly.min.js`) and
  pointed the `<script>` at it instead of `cdn.plot.ly`; removed the Google
  Fonts `<link>` (system fonts fall back). The page loads no external resource.
- Extracted the statistical core into a pure `engine.js` (single source of
  truth; the inline duplicates were deleted and the page now loads `engine.js`).
- Added `tests.js` (42 assertions, all passing) with hand-derived expectations.
- **Fixed a `k=1` / two-identical NaN bug** in the DL pooler: `τ²=(Q−df)/C` and
  `I²=(Q−df)/Q` both went `0/0 → NaN` and poisoned the estimate; they now
  degrade to the correct fixed-effect result (`τ²=0`, `I²=0`). All other math,
  including the normal CDF (correctly 0.5 at x=0), was verified and left intact.
- Added Pages scaffold (`.nojekyll`, README, `.gitignore`); renamed
  `Omni.html` → `index.html`. LICENSE (Apache 2.0) was already present and kept.

## Body (E156 draft — CURRENT BODY)

Can one offline dashboard run many meta-analytic paradigms over a single
long-format dataset without diverging from a shared, auditable statistical core?
It ingests binary, continuous, diagnostic and dose-response arms and exposes
pairwise, multilevel and robust-variance pooling, network and component network
meta-analysis, diagnostic accuracy, GOSH, BIC model averaging and an MCMC
sampler. Every routine now lives in one pure engine that pools ratio metrics on
the log scale and computes DerSimonian–Laird τ² as max(0,(Q−(k−1))/C), I²,
confidence and prediction intervals, and a Knapp–Hartung floor. A revival audit
extracted that engine verbatim, then fixed a k=1 divide-by-zero that had
silently returned NaN for single-study and identical-study pools. Forty-two hand-derived assertions now lock the core,
including a worked DerSimonian–Laird example and a matrix inverse-identity
property. The honest read is exploratory: light REML, heuristic AUC and
stochastic panels mean outputs aid synthesis, not a pre-specified plan. The
contribution is a transparent, offline, test-backed multi-method workbench, not
a new estimator.

SUBMITTED: [ ]
