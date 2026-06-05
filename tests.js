// tests.js — pure-Node test harness for engine.js (786-M14 Omnibus Extended).
// Expected values are hand-derived INDEPENDENTLY of the engine (see comments),
// not produced by running the code under test.
//   run:  node tests.js   ->  prints "N passed, M failed"; exit(M===0?0:1)

const E = require('./engine.js');
const { Stat, Matrix, Multilevel, RVE, NMA_RE, Pooling } = E;

let pass = 0, fail = 0;
const approx = (a, b, tol = 1e-6) => Math.abs(a - b) <= tol;
function check(name, cond, got, want) {
    if (cond) { pass++; }
    else { fail++; console.log(`  FAIL: ${name}` + (got !== undefined ? `  got=${got} want=${want}` : '')); }
}
function near(name, got, want, tol = 1e-6) { check(name, approx(got, want, tol), got, want); }

// ---------------------------------------------------------------------------
// 1. normalCDF sanity (a 0 at x=0 would be a BUG). Phi(0)=0.5; Phi(1.96)=0.975.
// ---------------------------------------------------------------------------
near('pnorm(0) === 0.5 (NOT 0)', Stat.pnorm(0), 0.5, 1e-9);
near('pnorm(1.959964)≈0.975', Stat.pnorm(1.959963985), 0.975, 5e-7);
near('pnorm(-1.959964)≈0.025', Stat.pnorm(-1.959963985), 0.025, 5e-7);
check('pnorm monotone', Stat.pnorm(1) > Stat.pnorm(0) && Stat.pnorm(0) > Stat.pnorm(-1), undefined);
// qnorm is the inverse of pnorm: qnorm(0.975) ≈ 1.959964
near('qnorm(0.975)≈1.959964', Stat.qnorm(0.975), 1.959963985, 1e-4);
near('qnorm(0.5)≈0', Stat.qnorm(0.5), 0, 1e-6);

// ---------------------------------------------------------------------------
// 2. HAND-WORKED DerSimonian-Laird pooling (the named primary estimand).
//    Studies (log scale): te=[0, 0.4, -0.2], se=0.2 each => vi=0.04, w_FE=25.
//      sW=75; mu_FE=(25*0+25*0.4+25*-0.2)/75 = 5/75 = 0.0666667
//      Q = 25*(0-mu)^2 + 25*(0.4-mu)^2 + 25*(-0.2-mu)^2
//        = 0.111111 + 2.777778 + 1.777778 = 4.6666667 ; df=2
//      C = sW - sum(w^2)/sW = 75 - 3*625/75 = 75 - 25 = 50
//      tau2 = (Q-df)/C = 2.6666667/50 = 0.05333333
//      I2 = (Q-df)/Q = 2.6666667/4.6666667 = 0.5714286 (57.14%)
//      w_RE = 1/(0.04+0.05333333) = 1/0.09333333 = 10.7142857 each; sWr=32.1428571
//      est_RE = 10.7142857*(0+0.4-0.2)/32.1428571 = 2.1428571/32.1428571 = 0.0666667
//      se_RE = sqrt(1/32.1428571) = sqrt(0.03111111) = 0.17638342
// ---------------------------------------------------------------------------
{
    const eff = [0, 0.4, -0.2].map(te => ({ te, se: 0.2, w: 1 / 0.04 }));
    const r = Pooling.pool(eff, 'dl');
    near('DL tau2 = 0.05333333', r.tau2, 0.0533333333, 1e-7);
    near('DL Q = 4.6666667', r.Q, 4.6666666667, 1e-7);
    near('DL I2 = 0.5714286', r.i2, 0.5714285714, 1e-7);
    near('DL est = 0.0666667', r.est, 0.0666666667, 1e-7);
    near('DL se = 0.17638342', r.se, 0.1763834207, 1e-7);
    // CI = est ± 1.96*se = 0.0666667 ± 0.3456115  => [-0.2789448, 0.4122782]
    near('DL ci.lo', r.lo, 0.0666666667 - 1.96 * 0.1763834207, 1e-7);
    near('DL ci.hi', r.hi, 0.0666666667 + 1.96 * 0.1763834207, 1e-7);
}

// ---------------------------------------------------------------------------
// 3. HKSJ (Knapp-Hartung) floor: adj_se = se_RE * sqrt(max(1, Q/df)),
//    crit = t_{df}. Here Q/df = 4.6666667/2 = 2.3333333 > 1 => floor active.
//      adj_se = 0.17638342 * sqrt(2.3333333) = 0.17638342*1.5275252 = 0.26943
//    (we re-derive adj_se independently and assert it matches.)
// ---------------------------------------------------------------------------
{
    const eff = [0, 0.4, -0.2].map(te => ({ te, se: 0.2, w: 1 / 0.04 }));
    const r = Pooling.pool(eff, 'hk');
    const adj = 0.1763834207 * Math.sqrt(2.3333333333);
    near('HKSJ adj_se floor max(1,Q/df)', r.se, adj, 1e-6);
    // HKSJ widens vs plain RE (since Q/df>1): hi-lo larger
    const re = Pooling.pool(eff, 'dl');
    check('HKSJ wider than RE-z when Q/df>1', (r.hi - r.lo) > (re.hi - re.lo), undefined);
}

// ---------------------------------------------------------------------------
// 4. EDGE: k=1 single study. Q=0, df=0; tau2 guard -> max(0, .../0).
//    est == te, no heterogeneity. Must not throw / produce NaN est.
// ---------------------------------------------------------------------------
{
    const r = Pooling.pool([{ te: 0.5, se: 0.3, w: 1 / 0.09 }], 'dl');
    near('k=1 est == te', r.est, 0.5, 1e-9);
    check('k=1 est finite', Number.isFinite(r.est), undefined);
}

// ---------------------------------------------------------------------------
// 5. EDGE: two IDENTICAL studies => tau2 = 0 and I2 = 0 (Q = 0 exactly).
//    te=[0.3,0.3], se=0.2 each. mu=0.3; Q=0; tau2=max(0,(0-1)/C)=0; I2=0.
//    se_RE = sqrt(1/(2*25)) = sqrt(1/50) = 0.14142136
// ---------------------------------------------------------------------------
{
    const eff = [0.3, 0.3].map(te => ({ te, se: 0.2, w: 1 / 0.04 }));
    const r = Pooling.pool(eff, 'dl');
    near('two-identical tau2 = 0', r.tau2, 0, 1e-12);
    near('two-identical I2 = 0', r.i2, 0, 1e-12);
    near('two-identical est = 0.3', r.est, 0.3, 1e-12);
    near('two-identical se = sqrt(1/50)', r.se, Math.sqrt(1 / 50), 1e-9);
}

// ---------------------------------------------------------------------------
// 6. EDGE: empty input guard for Egger's test (n<3 returns neutral {int:0,p:1}).
// ---------------------------------------------------------------------------
{
    const e0 = Pooling.egger([]);
    check('egger empty -> int 0', e0.int === 0, e0.int, 0);
    check('egger empty -> p 1', e0.p === 1, e0.p, 1);
    const e2 = Pooling.egger([{ te: 0.1, se: 0.2 }, { te: 0.2, se: 0.3 }]);
    check('egger n<3 -> neutral', e2.int === 0 && e2.p === 1, undefined);
}

// ---------------------------------------------------------------------------
// 7. MATRIX property: A * A^{-1} == I  (3x3 invertible matrix).
//    A = [[2,1,1],[1,3,2],[1,0,0]] ; det = 1 (computed by hand).
// ---------------------------------------------------------------------------
{
    const A = [[2, 1, 1], [1, 3, 2], [1, 0, 0]];
    const Ai = Matrix.inv(A);
    const P = Matrix.dot(A, Ai);
    let ok = true;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
        if (!approx(P[i][j], i === j ? 1 : 0, 1e-9)) ok = false;
    }
    check('Matrix A*inv(A) == I', ok, undefined);
    // transpose property: (A^T)_{ij} = A_{ji}
    const At = Matrix.t(A);
    check('Matrix transpose', At[0][1] === A[1][0] && At[2][0] === A[0][2], undefined);
    // dot dimension correctness: [2x3]·[3x2] -> [2x2]
    const D = Matrix.dot([[1, 2, 3], [4, 5, 6]], [[1, 0], [0, 1], [1, 1]]);
    check('Matrix dot dims', D.length === 2 && D[0].length === 2, undefined);
    near('Matrix dot value D[0][0]=1*1+2*0+3*1=4', D[0][0], 4, 1e-12);
}

// ---------------------------------------------------------------------------
// 8. NMA_RE property checks: covariance Vb symmetric (PSD-ish), tau2>=0, Q>=0,
//    and the reference contrast is 0. Two-study network: S1{Plac,A}, S2{Plac,A}.
//    With identical contrasts (logOR), the A-vs-Plac estimate must be finite.
// ---------------------------------------------------------------------------
{
    // contrasts vs ref 'Plac'; t1=ref, t2=A. te=log OR, var given.
    const con = [
        { study: 'S1', t1: 'Plac', t2: 'A', te: 0.5, var: 0.1, cov: 0 },
        { study: 'S2', t1: 'Plac', t2: 'A', te: 0.5, var: 0.1, cov: 0 },
    ];
    const r = NMA_RE.solve(con, ['Plac', 'A'], 'Plac', { method: 'dl' });
    check('NMA tau2 >= 0', r.tau2 >= 0, r.tau2);
    check('NMA Q >= 0', r.Q >= 0, r.Q);
    near('NMA ref contrast est == 0', r.est['Plac'].est, 0, 1e-12);
    // Two identical contrasts (te=0.5) => pooled A estimate = 0.5
    near('NMA pooled A est = 0.5', r.est['A'].est, 0.5, 1e-6);
    // Vb symmetric (1x1 here trivially) and SE positive
    check('NMA A se > 0', r.est['A'].se > 0, undefined);
    // I2 = max(0,(Q-df)/Q)*100 in [0,100]
    check('NMA I2 in [0,100]', r.I2 >= 0 && r.I2 <= 100, r.I2);
}

// ---------------------------------------------------------------------------
// 9. RVE: with one cluster per study and 2 clusters, df = m-1 = 1; t-crit large.
//    est == inverse-variance weighted mean. te=[0,0.4], se=0.2 => w=25 each;
//    mu = (0+0.4*... )/... ; symmetric weights => mu = 0.2.
// ---------------------------------------------------------------------------
{
    const r = RVE.fit([0, 0.4], [0.04, 0.04], ['c1', 'c2']);
    near('RVE est = IVW mean 0.2', r.est, 0.2, 1e-9);
    check('RVE se finite > 0', Number.isFinite(r.se) && r.se > 0, undefined);
    check('RVE CI ordered lo<hi', r.ci[0] < r.ci[1], undefined);
}

// ---------------------------------------------------------------------------
// 10. Multilevel.fit smoke: 4 effects in 2 clusters; partitions I2 into within
//     + between; total I2 must be in [0,100] and components non-negative.
// ---------------------------------------------------------------------------
{
    const r = Multilevel.fit([0.1, 0.3, -0.1, 0.5], [0.05, 0.05, 0.05, 0.05], ['c1', 'c1', 'c2', 'c2']);
    check('Multilevel sig2_w >= 0', r.sig2_w >= 0, r.sig2_w);
    check('Multilevel sig2_b >= 0', r.sig2_b >= 0, r.sig2_b);
    check('Multilevel I2_total in [0,100]', r.I2_total >= 0 && r.I2_total <= 100, r.I2_total);
    check('Multilevel est finite', Number.isFinite(r.est), undefined);
    check('Multilevel CI ordered', r.ci[0] < r.ci[1], undefined);
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail === 0 ? 0 : 1);
