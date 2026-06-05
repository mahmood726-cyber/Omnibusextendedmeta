// engine.js — pure statistical core for 786-M14 Omnibus Extended meta-analysis.
// Extracted VERBATIM from the inline app script (single source of truth).
// Pure functions/objects only — no DOM / Plotly dependencies:
//   Math.erf            : Abramowitz-Stegun error-function polynomial
//   Stat                : pnorm (normal CDF), qnorm (inverse), qt (t-quantile approx), invLogit
//   Matrix              : dense matrix algebra (zeros/transpose/dot/inv via Gauss-Jordan)
//   Optim.nelderMead    : derivative-free simplex minimiser
//   Multilevel.fit      : 3-level (within/between cluster) random-effects IGLS fit + Cheung I2
//   RVE.fit             : robust variance (cluster sandwich) estimator with t-critical df=m-1
//   NMA_RE.solve        : network meta-analysis (fixed / DL / REML) via GLS on contrast design
//   NMA_Multilevel.solve: multilevel network meta-analysis (within+between tau2)
//   pooling             : pairwise log-scale DL/HKSJ pooler + Egger's small-study test
// Effect sizes are pooled on the log scale and back-transformed in the UI.

if(!Math.erf){Math.erf=function(x){var s=(x>=0)?1:-1;x=Math.abs(x);var t=1/(1+0.3275911*x);return s*(1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*Math.exp(-x*x))}}
const Stat={pnorm:z=>0.5*(1+Math.erf(z/Math.sqrt(2))),qnorm:p=>{var a1=-39.6968302866538,a2=220.946098424521,a3=-275.928510446969,a4=138.357751867269,a5=-30.6647980661472,a6=2.50662827745924;var b1=-54.4760987982241,b2=161.585836858041,b3=-155.698979859887,b4=66.8013118877197,b5=-13.2806815528857;var c1=-7.78489400243029e-03,c2=-3.22396458041136e-01,c3=-2.40075827716184,c4=-2.54973253934373,c5=4.37466414146497,c6=2.93816398269878;var d1=7.78469570904146e-03,d2=3.22467129070040e-01,d3=2.44513413714300,d4=3.75440866190742;var q,r;if(p<0.02425){q=Math.sqrt(-2*Math.log(p));return(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)}else if(p>1-0.02425){q=Math.sqrt(-2*Math.log(1-p));return-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)}else{q=p-0.5;r=q*q;return(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)}},qt:(p,df)=>{const x=Stat.qnorm(p),g=(x*x+1)/4;return x*(1+(g/df)+(g*g*5/96)/(df*df))},invLogit:x=>Math.exp(x)/(1+Math.exp(x))};
const Matrix={zeros:(r,c)=>Array(r).fill(0).map(()=>Array(c).fill(0)),t:m=>m[0].map((_,i)=>m.map(r=>r[i])),dot:(a,b)=>{const r1=a.length,c1=a[0].length,c2=b[0].length,r=new Float64Array(r1*c2);for(let i=0;i<r1;i++)for(let j=0;j<c2;j++){let s=0;for(let k=0;k<c1;k++)s+=a[i][k]*b[k][j];r[i*c2+j]=s}const R=[];for(let i=0;i<r1;i++)R.push(Array.from(r.slice(i*c2,(i+1)*c2)));return R},inv:m=>{const n=m.length,A=m.map(r=>[...r]),I=Matrix.zeros(n,n);for(let i=0;i<n;i++)I[i][i]=1;for(let i=0;i<n;i++){let p=i,mx=Math.abs(A[i][i]);for(let k=i+1;k<n;k++)if(Math.abs(A[k][i])>mx){p=k;mx=Math.abs(A[k][i])}[A[i],A[p]]=[A[p],A[i]];[I[i],I[p]]=[I[p],I[i]];const d=A[i][i];if(Math.abs(d)<1e-12)continue;for(let j=0;j<n;j++){A[i][j]/=d;I[i][j]/=d}for(let k=0;k<n;k++)if(k!==i){const f=A[k][i];for(let j=0;j<n;j++){A[k][j]-=f*A[i][j];I[k][j]-=f*I[i][j]}}}return I}};
const Optim={nelderMead:(f,x0)=>{let s=x0.map((v,i)=>{let x=[...x0];x[i]+=0.05;return x});s.unshift(x0);let fs=s.map(f);for(let k=0;k<300;k++){const idx=fs.map((v,i)=>[v,i]).sort((a,b)=>a[0]-b[0]).map(x=>x[1]),b=s[idx[0]],w=s[idx[s.length-1]],c=b.map((_,i)=>s.reduce((a,x,j)=>j===idx[s.length-1]?a:a+x[i],0)/(s.length-1)),r=c.map((v,i)=>v+(v-w[i])),fr=f(r);if(fr<fs[idx[0]]){const e=c.map((v,i)=>v+2*(v-w[i])),fe=f(e);if(fe<fr){s[idx[s.length-1]]=e;fs[idx[s.length-1]]=fe}else{s[idx[s.length-1]]=r;fs[idx[s.length-1]]=fr}}else if(fr<fs[idx[s.length-2]]){s[idx[s.length-1]]=r;fs[idx[s.length-1]]=fr}else{const ct=c.map((v,i)=>v+0.5*(w[i]-v)),fc=f(ct);if(fc<fs[idx[s.length-1]]){s[idx[s.length-1]]=ct;fs[idx[s.length-1]]=fc}else{for(let i=1;i<s.length;i++){s[idx[i]]=s[idx[i]].map((v,j)=>b[j]+0.5*(v-b[j]));fs[idx[i]]=f(s[idx[i]])}}}}let mv=Infinity,mi=0;for(let i=0;i<fs.length;i++)if(!isNaN(fs[i])&&fs[i]<mv){mv=fs[i];mi=i}return s[mi]}};

// --- MULTILEVEL & RVE MODULES ---
const Multilevel = {
    fit: (Y, V, cluster, opts = {}) => {
        const n = Y.length;
        const maxIter = opts.maxIter || 100;
        const tol = opts.tol || 1e-6;
        const clusters = [...new Set(cluster)];
        const m = clusters.length;

        const muInit = Y.reduce((a, b) => a + b, 0) / n;
        const ssTotal = Y.reduce((a, y) => a + (y - muInit) ** 2, 0);
        let sig2_w = Math.max(0.01, ssTotal / n - V.reduce((a, b) => a + b, 0) / n);
        let sig2_b = sig2_w / 2;

        let mu = muInit, prevLL = -Infinity;

        for (let iter = 0; iter < maxIter; iter++) {
            const Sigma = Matrix.zeros(n, n);
            for (let i = 0; i < n; i++) {
                for (let j = 0; j < n; j++) {
                    if (i === j) Sigma[i][j] = V[i] + sig2_w + sig2_b;
                    else if (cluster[i] === cluster[j]) Sigma[i][j] = sig2_b;
                }
            }
            const W = Matrix.inv(Sigma);
            let sumW = 0, sumWY = 0;
            for (let i = 0; i < n; i++) {
                for (let j = 0; j < n; j++) {
                    sumW += W[i][j];
                    sumWY += W[i][j] * Y[j];
                }
            }
            mu = sumWY / sumW;
            const r = Y.map(y => y - mu);
            const P = W.map((row, i) => row.map((w, j) => {
                let adj = 0;
                for (let a = 0; a < n; a++) adj += W[i][a];
                for (let b = 0; b < n; b++) adj *= W[b][j] / sumW;
                return w - adj;
            }));
            let tr_Pw = 0, tr_Pb = 0, tr_Pww = 0, tr_Pbb = 0;
            let rPr = 0;
            for (let i = 0; i < n; i++) {
                tr_Pw += P[i][i];
                rPr += r[i] * P[i].reduce((a, p, j) => a + p * r[j], 0);
                for (let j = 0; j < n; j++) {
                    tr_Pww += P[i][j] * P[j][i];
                    if (cluster[i] === cluster[j]) {
                        tr_Pb += P[i][j];
                        tr_Pbb += P[i].reduce((a, p, k) => a + (cluster[k] === cluster[j] ? p * P[k][j] : 0), 0);
                    }
                }
            }
            const score_w = -0.5 * tr_Pw + 0.5 * r.reduce((a, ri, i) => a + ri * P[i].reduce((s, p, j) => s + p * r[j], 0), 0);
            const score_b = -0.5 * tr_Pb + 0.5 * r.reduce((a, ri, i) => a + ri * P[i].reduce((s, p, j) => s + (cluster[i] === cluster[j] ? p * r[j] : 0), 0), 0);
            const step = 0.5;
            sig2_w = Math.max(0, sig2_w + step * score_w / Math.max(1, Math.abs(tr_Pww)));
            sig2_b = Math.max(0, sig2_b + step * score_b / Math.max(1, Math.abs(tr_Pbb)));
            const ll = -0.5 * (n * Math.log(2 * Math.PI) + rPr); // Approx LL
            if (Math.abs(ll - prevLL) < tol) break;
            prevLL = ll;
        }
        const totalVar = sig2_w + sig2_b + V.reduce((a, b) => a + b, 0) / n;
        const Sigma = Matrix.zeros(n, n);
        for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if(i===j) Sigma[i][j] = V[i] + sig2_w + sig2_b; else if(cluster[i] === cluster[j]) Sigma[i][j] = sig2_b;
        const W = Matrix.inv(Sigma);
        let sumW = 0; for(let i=0;i<n;i++)for(let j=0;j<n;j++) sumW+=W[i][j];

        const se = Math.sqrt(1/sumW);
        const pi_se = Math.sqrt(se**2 + sig2_w + sig2_b);
        const crit = 1.96;

        return {
            est: mu,
            se: se,
            ci: [mu - crit*se, mu + crit*se],
            pi: [mu - crit*pi_se, mu + crit*pi_se],
            I2_2: (sig2_w / totalVar) * 100,
            I2_3: (sig2_b / totalVar) * 100,
            I2_total: ((sig2_w+sig2_b)/totalVar)*100,
            sig2_w, sig2_b
        };
    }
};

const RVE = {
    fit: (Y, V, cluster) => {
        const n = Y.length, clusters = [...new Set(cluster)], m = clusters.length;
        const w = V.map(v => 1 / v), sumW = w.reduce((a, b) => a + b, 0);
        const mu = Y.reduce((a, y, i) => a + w[i] * y, 0) / sumW;
        const r = Y.map(y => y - mu);
        const bread = 1 / (sumW ** 2);
        let meat = 0;
        clusters.forEach(c => {
            const idx = cluster.map((cl, i) => cl === c ? i : -1).filter(i => i >= 0);
            const clusterSum = idx.reduce((a, i) => a + w[i] * r[i], 0);
            meat += clusterSum ** 2;
        });
        const varRobust = bread * meat * (m / (m - 1));
        const se = Math.sqrt(varRobust);
        const df = m - 1;
        const crit = Stat.qt(0.975, df);
        return {
            est: mu,
            se: se,
            ci: [mu - crit*se, mu + crit*se]
        };
    }
};

const NMA_RE = {
    solve: (con, txs, ref, opts = {}) => {
        const method = opts.method || 'dl';
        const lst = txs.filter(t => t !== ref), k = con.length, p = lst.length;
        const Y = con.map(c => [c.te]), X = Matrix.zeros(k, p);
        con.forEach((c, i) => {
            const i2 = lst.indexOf(c.t2), i1 = lst.indexOf(c.t1);
            if (i2 !== -1) X[i][i2] = 1; if (i1 !== -1) X[i][i1] = -1;
        });
        const sm = {}; con.forEach((c, i) => { if (!sm[c.study]) sm[c.study] = []; sm[c.study].push(i); });
        const buildV = (tau2) => {
            const V = Matrix.zeros(k, k);
            Object.values(sm).forEach(idx => {
                if (idx.length === 1) V[idx[0]][idx[0]] = con[idx[0]].var + tau2;
                else { const cv = con[idx[0]].cov; idx.forEach(i => idx.forEach(j => { V[i][j] = (i === j) ? con[i].var + tau2 : cv; })); }
            }); return V;
        };
        let tau2 = 0;
        // Simple Q calc for DL
        const V0 = buildV(0), W0 = Matrix.inv(V0), Xt = Matrix.t(X);
        const Vb0 = Matrix.inv(Matrix.dot(Matrix.dot(Xt, W0), X)), b0 = Matrix.dot(Matrix.dot(Vb0, Matrix.dot(Xt, W0)), Y);
        let Q = 0; con.forEach((c, i) => {
            const i2 = lst.indexOf(c.t2), i1 = lst.indexOf(c.t1); let pred = 0;
            if (i2 !== -1) pred += b0[i2][0]; if (i1 !== -1) pred -= b0[i1][0];
            Q += (c.te - pred) ** 2 / c.var;
        });
        const df = k - p;
        if (method === 'dl' && Q > df) {
            const sumW = con.reduce((a, c) => a + 1 / c.var, 0); const sumW2 = con.reduce((a, c) => a + 1 / c.var ** 2, 0);
            tau2 = Math.max(0, (Q - df) / (sumW - sumW2 / sumW));
        }
        if (method === 'reml') {
            for (let iter = 0; iter < 20; iter++) {
                const V = buildV(tau2), W = Matrix.inv(V);
                const Vb = Matrix.inv(Matrix.dot(Matrix.dot(Xt, W), X)), b = Matrix.dot(Matrix.dot(Vb, Matrix.dot(Xt, W)), Y);
                let rss = 0; con.forEach((c, i) => {
                    const i2 = lst.indexOf(c.t2), i1 = lst.indexOf(c.t1); let pred = 0;
                    if (i2 !== -1) pred += b[i2][0]; if (i1 !== -1) pred -= b[i1][0];
                    rss += (c.te - pred) ** 2 * W[i][i];
                });
                const newTau2 = Math.max(0, tau2 + (rss - (k - p)) / (k - p));
                if (Math.abs(newTau2 - tau2) < 1e-6) break; tau2 = newTau2;
            }
        }
        const V = buildV(tau2), W = Matrix.inv(V), Vb = Matrix.inv(Matrix.dot(Matrix.dot(Xt, W), X)), b = Matrix.dot(Matrix.dot(Vb, Matrix.dot(Xt, W)), Y);
        const est = {}; est[ref] = { est: 0, se: 0 }; lst.forEach((t, i) => { est[t] = { est: b[i][0], se: Math.sqrt(Vb[i][i]) }; });
        return { est, Vb, lst, ref, VFull: V, tau2, Q, df, I2: Math.max(0, (Q - df) / Q) * 100 };
    }
};

const NMA_Multilevel = {
    solve: (con, txs, ref) => {
        const lst = txs.filter(t => t !== ref), k = con.length, p = lst.length;
        const Y = con.map(c => [c.te]), X = Matrix.zeros(k, p);
        con.forEach((c, i) => {
            const i2 = lst.indexOf(c.t2), i1 = lst.indexOf(c.t1);
            if (i2 !== -1) X[i][i2] = 1; if (i1 !== -1) X[i][i1] = -1;
        });
        let tau2_w = 0.05, tau2_b = 0.05;
        for (let iter = 0; iter < 10; iter++) {
            const V = Matrix.zeros(k, k);
            for (let i = 0; i < k; i++) {
                for (let j = 0; j < k; j++) {
                    if (i === j) V[i][j] = con[i].var + tau2_w + tau2_b;
                    else if (con[i].study === con[j].study) V[i][j] = con[i].cov + tau2_b;
                    else if ((con[i].cluster||con[i].study) === (con[j].cluster||con[j].study)) V[i][j] = tau2_b;
                }
            }
            const W = Matrix.inv(V), Xt = Matrix.t(X), Vb = Matrix.inv(Matrix.dot(Matrix.dot(Xt, W), X)), b = Matrix.dot(Matrix.dot(Vb, Matrix.dot(Xt, W)), Y);
            let ss_w = 0, n_w = 0;
            con.forEach((c, i) => {
                const i2 = lst.indexOf(c.t2), i1 = lst.indexOf(c.t1); let pred = 0;
                if (i2 !== -1) pred += b[i2][0]; if (i1 !== -1) pred -= b[i1][0];
                ss_w += (c.te - pred) ** 2 - c.var; n_w++;
            });
            tau2_w = Math.max(0, ss_w / n_w - tau2_b);
        }
        const V = Matrix.zeros(k, k); // Final V
        for(let i=0;i<k;i++)for(let j=0;j<k;j++) if(i===j)V[i][j]=con[i].var+tau2_w+tau2_b; else if(con[i].study===con[j].study)V[i][j]=con[i].cov+tau2_b; else if((con[i].cluster||con[i].study)===(con[j].cluster||con[j].study))V[i][j]=tau2_b;
        const W = Matrix.inv(V), Xt = Matrix.t(X), Vb = Matrix.inv(Matrix.dot(Matrix.dot(Xt, W), X)), b = Matrix.dot(Matrix.dot(Vb, Matrix.dot(Xt, W)), Y);
        const est = {}; est[ref] = { est: 0, se: 0 }; lst.forEach((t, i) => { est[t] = { est: b[i][0], se: Math.sqrt(Vb[i][i]) }; });
        return { est, Vb, lst, ref, VFull: V, tau2: tau2_w+tau2_b, Q: 0, I2: 0, tau2_w, tau2_b };
    }
};

// --- PAIRWISE POOLING (pure) ---
// Random/fixed-effects pooler on the (already log-transformed) effect sizes.
// eff = [{te, se, w:1/se^2}, ...]; mod: 'fixed' | 'dl' | 'hk' (Knapp-Hartung) | 'mh'
const Pooling = {
    pool:(eff,mod)=>{
        let w=eff.map(e=>e.w);const sW=w.reduce((a,b)=>a+b,0),mu=eff.reduce((a,e,i)=>a+e.te*w[i],0)/sW;
        // k<2 guard: df=0 makes the DL denominator C=sW-sum(w^2)/sW=0 and Q=0,
        // so tau2=(Q-df)/C = 0/0 = NaN which poisons est. With one study tau2 is
        // undefined; the correct degenerate result is the fixed-effect estimate (tau2=0).
        const Q=eff.reduce((a,e,i)=>a+w[i]*(e.te-mu)**2,0),df=eff.length-1;
        const C=sW-w.reduce((a,b)=>a+b*b,0)/sW;
        const tau2=(df<1||C<=0)?0:Math.max(0,(Q-df)/C);
        if(mod.includes('h')||mod==='dl')w=eff.map(e=>1/(e.se**2+tau2));
        const swr=w.reduce((a,b)=>a+b,0),est=w.reduce((a,x,i)=>a+x*eff[i].te,0)/swr,se=Math.sqrt(1/swr);
        let crit=1.96,adj=se; if(mod==='hk'){crit=Stat.qt(0.975,df);adj=se*Math.sqrt(Math.max(1,Q/df))}
        const pi=Math.sqrt(adj**2+tau2);
        const i2=(Q>0)?Math.max(0,(Q-df)/Q):0;   // Q=0 (k=1 or identical studies) => no heterogeneity
        return {est,se:adj,lo:est-crit*adj,hi:est+crit*adj,i2,pi_lo:est-crit*pi,pi_hi:est+crit*pi,tau2,Q,df};
    },
    egger:(eff)=>{
        const n=eff.length; if(n<3) return {int:0,p:1};
        const x=eff.map(e=>1/e.se), y=eff.map(e=>e.te/e.se);
        const mx=x.reduce((a,b)=>a+b,0)/n, my=y.reduce((a,b)=>a+b,0)/n;
        const b1=x.reduce((a,xi,i)=>a+(xi-mx)*(y[i]-my),0)/x.reduce((a,xi)=>a+(xi-mx)**2,0), b0=my-b1*mx;
        const sse=y.reduce((a,yi,i)=>a+(yi-(b0+b1*x[i]))**2,0), seRes=Math.sqrt(sse/(n-2));
        const seB0=seRes*Math.sqrt(1/n + (mx**2)/x.reduce((a,xi)=>a+(xi-mx)**2,0));
        const t=b0/seB0, p=2*(1-Stat.pnorm(Math.abs(t)));
        return {int:b0, p:p};
    }
};

if (typeof module!=='undefined'&&module.exports){ module.exports = { Stat, Matrix, Optim, Multilevel, RVE, NMA_RE, NMA_Multilevel, Pooling }; }
