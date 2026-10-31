/* chem.js — pure chemistry logic (no DOM). Works in browser (window.Chem) and node (module.exports).
 *
 * Molecule model:
 *   mol = { atoms:[{id,el,x,y,charge}], bonds:[{id,a,b,order,stereo}] }
 *   a,b are atom ids. order 1|2|3. stereo 0 = plain, 1 = wedge, 2 = dash (narrow end at atom a).
 */
(function (root) {
  'use strict';

  /* ---------------------------------------------------------------- data */
  const EL = {
    H:  { m: 1.008,   rc: 0.31, vdw: 1.10, c: '#6b7280', c3: '#f1f1f1' },
    B:  { m: 10.81,   rc: 0.84, vdw: 1.92, c: '#d97706', c3: '#f4a98c' },
    C:  { m: 12.011,  rc: 0.76, vdw: 1.70, c: '#1f2937', c3: '#555b66' },
    N:  { m: 14.007,  rc: 0.71, vdw: 1.55, c: '#2563eb', c3: '#3b6bff' },
    O:  { m: 15.999,  rc: 0.66, vdw: 1.52, c: '#dc2626', c3: '#ef3b3b' },
    F:  { m: 18.998,  rc: 0.57, vdw: 1.47, c: '#16a34a', c3: '#7ddf6a' },
    Si: { m: 28.085,  rc: 1.11, vdw: 2.10, c: '#a16207', c3: '#e6c28a' },
    P:  { m: 30.974,  rc: 1.07, vdw: 1.80, c: '#ea580c', c3: '#ff9a3c' },
    S:  { m: 32.06,   rc: 1.05, vdw: 1.80, c: '#ca8a04', c3: '#f2d13b' },
    Cl: { m: 35.45,   rc: 1.02, vdw: 1.75, c: '#15803d', c3: '#3fd13f' },
    Br: { m: 79.904,  rc: 1.20, vdw: 1.85, c: '#9f1239', c3: '#b5502a' },
    I:  { m: 126.904, rc: 1.39, vdw: 1.98, c: '#7e22ce', c3: '#9c27b0' },
    Se: { m: 78.971,  rc: 1.20, vdw: 1.90, c: '#b45309', c3: '#ffa100' },
    Li: { m: 6.94,    rc: 1.28, vdw: 1.82, c: '#6d28d9', c3: '#cc80ff' },
    Na: { m: 22.990,  rc: 1.66, vdw: 2.27, c: '#6d28d9', c3: '#ab5cf2' },
    K:  { m: 39.098,  rc: 2.03, vdw: 2.75, c: '#6d28d9', c3: '#8f40d4' },
    Mg: { m: 24.305,  rc: 1.41, vdw: 1.73, c: '#047857', c3: '#8aff00' },
    Ca: { m: 40.078,  rc: 1.76, vdw: 2.31, c: '#047857', c3: '#3dff00' },
    Al: { m: 26.982,  rc: 1.21, vdw: 1.84, c: '#6b7280', c3: '#bfa6a6' }
  };

  // allowed valences by element and formal charge
  const VT = {
    H:  { 0: [1], 1: [0], '-1': [0] },
    B:  { 0: [3], 1: [2], '-1': [4] },
    C:  { 0: [4], 1: [3], '-1': [3] },
    N:  { 0: [3], 1: [4], '-1': [2] },
    O:  { 0: [2], 1: [3], '-1': [1] },
    F:  { 0: [1], 1: [2], '-1': [0] },
    Si: { 0: [4], 1: [3], '-1': [3] },
    P:  { 0: [3, 5], 1: [4], '-1': [2] },
    S:  { 0: [2, 4, 6], 1: [3, 5], '-1': [1] },
    Se: { 0: [2, 4, 6], 1: [3, 5], '-1': [1] },
    Cl: { 0: [1], 1: [2], '-1': [0] },
    Br: { 0: [1], 1: [2], '-1': [0] },
    I:  { 0: [1, 3, 5], 1: [2, 4], '-1': [0] },
    Li: { 0: [1], 1: [0], '-1': [0] },
    Na: { 0: [1], 1: [0], '-1': [0] },
    K:  { 0: [1], 1: [0], '-1': [0] },
    Mg: { 0: [2], 1: [1], '-1': [1] },
    Ca: { 0: [2], 1: [1], '-1': [1] },
    Al: { 0: [3], 1: [2], '-1': [4] }
  };

  function valences(el, q) {
    const t = VT[el];
    if (!t) return [4];
    if (t[q] !== undefined) return t[q];
    return t[0];
  }

  function hydrogens(a, bos) {
    const vals = valences(a.el, a.charge || 0);
    for (const v of vals) if (v >= bos) return { h: v - bos, ok: true };
    return { h: 0, ok: false };
  }

  /* --------------------------------------------------------------- graph */
  function buildGraph(mol) {
    const n = mol.atoms.length, idx = new Map();
    mol.atoms.forEach((a, i) => idx.set(a.id, i));
    const nbrs = Array.from({ length: n }, () => []);
    const bonds = [];
    mol.bonds.forEach((b, mi) => {
      const i = idx.get(b.a), j = idx.get(b.b);
      if (i === undefined || j === undefined || i === j) return;
      const bi = bonds.length;
      bonds.push({ i, j, order: b.order, stereo: b.stereo || 0, mi });
      nbrs[i].push({ j, bi, order: b.order });
      nbrs[j].push({ j: i, bi, order: b.order });
    });
    const bos = nbrs.map(l => l.reduce((s, e) => s + e.order, 0));
    const hc = [], ok = [];
    mol.atoms.forEach((a, i) => { const r = hydrogens(a, bos[i]); hc.push(r.h); ok.push(r.ok); });
    return { n, idx, nbrs, bonds, bos, hc, ok };
  }

  /* ------------------------------------------------------ formula & mass */
  function formula(mol) {
    const g = buildGraph(mol);
    const counts = {}; let q = 0;
    mol.atoms.forEach((a, i) => {
      counts[a.el] = (counts[a.el] || 0) + 1;
      if (g.hc[i]) counts.H = (counts.H || 0) + g.hc[i];
      q += a.charge || 0;
    });
    let keys = Object.keys(counts);
    if (counts.C) {
      keys = ['C'].concat(keys.filter(k => k === 'H')).concat(keys.filter(k => k !== 'C' && k !== 'H').sort());
    } else keys.sort();
    let plain = '', html = '', mass = 0;
    keys.forEach(k => {
      const c = counts[k];
      plain += k + (c > 1 ? c : '');
      html += k + (c > 1 ? '<sub>' + c + '</sub>' : '');
      mass += c * (EL[k] ? EL[k].m : 0);
    });
    const sup = ['', '⁺', '²⁺', '³⁺'], sdn = ['', '⁻', '²⁻', '³⁻'];
    if (q > 0) { plain += sup[Math.min(q, 3)]; html += '<sup>' + (q > 1 ? q : '') + '+</sup>'; }
    if (q < 0) { plain += sdn[Math.min(-q, 3)]; html += '<sup>' + (q < -1 ? -q : '') + '−</sup>'; }
    // degree of unsaturation (only meaningful for neutral C/H/N/O/S/P/halogen species)
    let dbe = null;
    const allowed = new Set(['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I']);
    if (q === 0 && keys.every(k => allowed.has(k)) && counts.C) {
      const X = (counts.F || 0) + (counts.Cl || 0) + (counts.Br || 0) + (counts.I || 0);
      const v = (2 * counts.C + 2 + (counts.N || 0) + (counts.P || 0) - (counts.H || 0) - X) / 2;
      if (Number.isInteger(v) && v >= 0) dbe = v;
    }
    return { counts, formula: plain, formulaHtml: html, mass, charge: q, dbe, heavy: mol.atoms.length, hydrogens: counts.H || 0 };
  }

  /* -------------------------------------------------------------- rings */
  function findRings(g) {
    const n = g.n, E = g.bonds.length;
    const comp = new Array(n).fill(-1); let C = 0;
    for (let s = 0; s < n; s++) {
      if (comp[s] >= 0) continue;
      const st = [s]; comp[s] = C;
      while (st.length) { const u = st.pop(); for (const e of g.nbrs[u]) if (comp[e.j] < 0) { comp[e.j] = C; st.push(e.j); } }
      C++;
    }
    const need = E - n + C;
    const rings = [], bondInRing = new Array(E).fill(false);
    if (need <= 0) return { rings, systems: [], ringSys: [], bondInRing };

    const cand = new Map();
    for (let bi = 0; bi < E; bi++) {
      const u = g.bonds[bi].i, v = g.bonds[bi].j;
      const prev = new Array(n).fill(-1), pe = new Array(n).fill(-1);
      const q = [u]; prev[u] = u; let found = false;
      while (q.length && !found) {
        const x = q.shift();
        for (const e of g.nbrs[x]) {
          if (e.bi === bi || prev[e.j] !== -1) continue;
          prev[e.j] = x; pe[e.j] = e.bi;
          if (e.j === v) { found = true; break; }
          q.push(e.j);
        }
      }
      if (!found) continue;
      const atoms = [], bonds = [bi];
      let x = v;
      while (x !== u) { atoms.push(x); bonds.push(pe[x]); x = prev[x]; }
      atoms.push(u); atoms.reverse();
      const key = bonds.slice().sort((a, b) => a - b).join(',');
      if (!cand.has(key)) cand.set(key, { atoms, bonds });
    }
    const list = Array.from(cand.values()).sort((a, b) => a.atoms.length - b.atoms.length);
    const basis = [];
    for (const c of list) {
      if (rings.length >= need) break;
      const v = new Uint8Array(E); c.bonds.forEach(b => { v[b] = 1; });
      for (const b of basis) if (v[b.pivot]) for (let k = 0; k < E; k++) v[k] ^= b.vec[k];
      const p = v.indexOf(1);
      if (p < 0) continue;
      basis.push({ vec: v, pivot: p });
      rings.push({ atoms: c.atoms, bonds: c.bonds, size: c.atoms.length });
    }
    rings.forEach(r => r.bonds.forEach(b => { bondInRing[b] = true; }));
    // fused systems: rings sharing at least one bond
    const par = rings.map((_, i) => i);
    const find = x => par[x] === x ? x : (par[x] = find(par[x]));
    for (let a = 0; a < rings.length; a++) for (let b = a + 1; b < rings.length; b++) {
      if (rings[a].bonds.some(x => rings[b].bonds.includes(x))) par[find(a)] = find(b);
    }
    const groups = new Map();
    rings.forEach((_, i) => { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); });
    const systems = Array.from(groups.values());
    const ringSys = new Array(rings.length);
    systems.forEach((s, si) => s.forEach(ri => { ringSys[ri] = si; }));
    return { rings, systems, ringSys, bondInRing };
  }

  /* --------------------------------------------------------- aromaticity */
  // pi-electron contribution of one ring atom, in the context of a (fused) ring system
  function piInfo(mol, g, i, sysBonds) {
    const a = mol.atoms[i], el = a.el, q = a.charge || 0, nb = g.nbrs[i];
    const dbl = nb.filter(x => x.order === 2), trp = nb.filter(x => x.order === 3);
    const sigma = nb.length + g.hc[i];
    if (trp.length) return trp.some(x => sysBonds.has(x.bi)) ? { e: 1 } : { block: 'sp carbon' };
    if (dbl.length > 1) return { block: 'cumulated double bonds' };
    if (dbl.length === 1) {
      const d = dbl[0];
      if (sysBonds.has(d.bi)) return { e: 1 };
      const p = mol.atoms[d.j].el;
      if (el === 'C' && (p === 'O' || p === 'S' || p === 'N' || p === 'Se')) return { e: 0, note: 'exocyclic C=' + p + ' (empty p orbital)' };
      return { block: 'exocyclic double bond (cross-conjugated)' };
    }
    if (el === 'C' || el === 'Si') {
      if (q === -1 && sigma === 3) return { e: 2, note: 'carbanion lone pair' };
      if (q === 1 && sigma === 3) return { e: 0, note: 'carbocation (empty p orbital)' };
      return { block: 'sp3 atom' };
    }
    if (el === 'B') return (q === 0 && sigma === 3) ? { e: 0, note: 'empty p orbital' } : { block: 'sp3 atom' };
    if (el === 'N' || el === 'P') return (q <= 0 && sigma <= 3) ? { e: 2, note: 'lone pair' } : { block: 'sp3 atom' };
    if (el === 'O' || el === 'S' || el === 'Se') return (q === 0 && sigma === 2) ? { e: 2, note: 'lone pair' } : { block: 'sp3 atom' };
    return { block: 'atom cannot take part in conjugation' };
  }

  function classifyRing(size, e, blocks, lab) {
    if (blocks.length) {
      const why = blocks.map(b => lab(b.a) + ' (' + b.why + ')').join(', ');
      return { type: 'nonaromatic', reason: 'Not fully conjugated: ' + why + ' interrupts the cyclic π system.' };
    }
    if (e === 0) return { type: 'nonaromatic', reason: 'No π electrons in the ring.' };
    if (e % 4 === 2) {
      const n = (e - 2) / 4;
      return { type: 'aromatic', reason: e + ' π electrons = 4n+2 with n = ' + n + '; cyclic, planar and fully conjugated → aromatic (Hückel).' };
    }
    if (e % 4 === 0) {
      if (size === 6) return { type: 'nonaromatic', reason: e + ' π electrons, but a 6-membered ring with empty p orbitals is quinoid / cross-conjugated, so it is not antiaromatic.' };
      if (size >= 8) return { type: 'nonaromatic', reason: e + ' π electrons (4n) but a ring of ' + size + ' atoms is not planar (tub shape, e.g. cyclooctatetraene), so it avoids antiaromaticity → non-aromatic.' };
      return { type: 'antiaromatic', reason: e + ' π electrons = 4n with n = ' + (e / 4) + '; cyclic, planar and fully conjugated → antiaromatic.' };
    }
    return { type: 'nonaromatic', reason: e + ' π electrons (odd): not 4n+2 and not 4n → non-aromatic.' };
  }

  function aromaticity(mol, g, rg) {
    g = g || buildGraph(mol);
    rg = rg || findRings(g);
    const lab = i => mol.atoms[i].el + (i + 1);
    const sysInfo = rg.systems.map(sys => {
      const sb = new Set(), atoms = new Set();
      sys.forEach(ri => { rg.rings[ri].bonds.forEach(b => sb.add(b)); rg.rings[ri].atoms.forEach(a => atoms.add(a)); });
      const pi = new Map(); atoms.forEach(a => pi.set(a, piInfo(mol, g, a, sb)));
      return { sys, atoms, pi };
    });
    const rings = rg.rings.map((r, ri) => {
      const S = sysInfo[rg.ringSys[ri]];
      let e = 0; const blocks = [], contrib = [];
      r.atoms.forEach(a => {
        const p = S.pi.get(a);
        if (p.block) blocks.push({ a, why: p.block });
        else { e += p.e; contrib.push({ a, label: lab(a), e: p.e, note: p.note || '' }); }
      });
      const c = classifyRing(r.size, e, blocks, lab);
      return { atoms: r.atoms, bonds: r.bonds, size: r.size, e, blocks, contrib, type: c.type, reason: c.reason, system: rg.ringSys[ri] };
    });
    // two-ring fused systems whose single rings have odd electron counts (azulene, pentalene …)
    sysInfo.forEach(S => {
      if (S.sys.length !== 2) return;
      const rs = S.sys.map(ri => rings[ri]);
      if (!rs.every(r => !r.blocks.length && r.e % 2 === 1)) return;
      let e = 0; S.atoms.forEach(a => { e += S.pi.get(a).e; });
      if (e % 4 === 2) rs.forEach(r => { r.type = 'aromatic'; r.reason = 'Rings share a conjugated bicyclic π system of ' + e + ' electrons = 4n+2 (n = ' + ((e - 2) / 4) + ') → aromatic as a fused system (azulene-type).'; });
      else if (e % 4 === 0 && e > 0) rs.forEach(r => { r.type = 'antiaromatic'; r.reason = 'Rings share a conjugated bicyclic π system of ' + e + ' electrons = 4n (n = ' + (e / 4) + ') → antiaromatic as a fused system (pentalene-type).'; });
    });
    const has = t => rings.some(r => r.type === t);
    let overall, message;
    if (!rings.length) { overall = 'acyclic'; message = 'No rings in this molecule, so it cannot be aromatic.'; }
    else if (has('aromatic') && has('antiaromatic')) { overall = 'mixed'; message = 'Contains both aromatic and antiaromatic rings.'; }
    else if (has('aromatic')) { overall = 'aromatic'; message = 'Contains ' + rings.filter(r => r.type === 'aromatic').length + ' aromatic ring(s).'; }
    else if (has('antiaromatic')) { overall = 'antiaromatic'; message = 'Contains antiaromatic ring(s) (destabilised).'; }
    else { overall = 'nonaromatic'; message = 'Rings are present but none is aromatic or antiaromatic.'; }
    return { rings, overall, message };
  }

  /* ----------------------------------------------------- geometry helpers */
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const len = a => Math.sqrt(dot(a, a));
  const unit = a => { const l = len(a) || 1e-9; return [a[0] / l, a[1] / l, a[2] / l]; };
  const rad = d => d * Math.PI / 180;

  function lcg(seed) { let s = seed >>> 0; return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296); }

  function ringAngleMap(rg) {
    const m = new Map();
    rg.rings.forEach(r => {
      const k = r.atoms.length;
      r.atoms.forEach((c, t) => {
        const p = r.atoms[(t - 1 + k) % k], nx = r.atoms[(t + 1) % k];
        const key = c + '|' + Math.min(p, nx) + '|' + Math.max(p, nx);
        const old = m.get(key);
        if (old === undefined || k < old) m.set(key, k);
      });
    });
    return m;
  }

  function topoDist(adj, maxd) {
    const n = adj.length, D = [];
    for (let s = 0; s < n; s++) {
      const d = new Int8Array(n).fill(-1); d[s] = 0; const q = [s];
      for (let h = 0; h < q.length; h++) {
        const u = q[h]; if (d[u] >= maxd) continue;
        for (const v of adj[u]) if (d[v] < 0) { d[v] = d[u] + 1; q.push(v); }
      }
      D.push(d);
    }
    return D;
  }

  /* generic damped-dynamics minimiser.  dist: [i,j,d0,k]   rep: [i,j,rmin,k]   chir: {c,a,b,d,t,k} */
  function relax(P, dist, rep, chir, o) {
    const n = P.length, dim = o.dim || 3, dt = o.dt || 0.01, iters = o.iters || 2000;
    const V = P.map(() => [0, 0, 0]), G = P.map(() => [0, 0, 0]);
    const tri = (c, a, b, d) => dot(unit(sub(a, c)), cross(unit(sub(b, c)), unit(sub(d, c))));
    for (let it = 0; it < iters; it++) {
      for (let i = 0; i < n; i++) { G[i][0] = 0; G[i][1] = 0; G[i][2] = 0; }
      for (let t = 0; t < dist.length; t++) {
        const [i, j, d0, k] = dist[t];
        const dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1], dz = dim === 3 ? P[i][2] - P[j][2] : 0;
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-6;
        const f = 2 * k * (d - d0) / d;
        G[i][0] += f * dx; G[i][1] += f * dy; G[i][2] += f * dz;
        G[j][0] -= f * dx; G[j][1] -= f * dy; G[j][2] -= f * dz;
      }
      for (let t = 0; t < rep.length; t++) {
        const [i, j, r0, k] = rep[t];
        const dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1], dz = dim === 3 ? P[i][2] - P[j][2] : 0;
        const d2 = dx * dx + dy * dy + dz * dz;
        if (d2 >= r0 * r0) continue;
        const d = Math.sqrt(d2) || 1e-6;
        const f = 2 * k * (d - r0) / d;
        G[i][0] += f * dx; G[i][1] += f * dy; G[i][2] += f * dz;
        G[j][0] -= f * dx; G[j][1] -= f * dy; G[j][2] -= f * dz;
      }
      if (dim === 3) {
        for (const ch of chir) {
          const ids = [ch.c, ch.a, ch.b, ch.d];
          const E = () => { const x = tri(P[ids[0]], P[ids[1]], P[ids[2]], P[ids[3]]) - ch.t; return ch.k * x * x; };
          for (const id of ids) for (let ax = 0; ax < 3; ax++) {
            const old = P[id][ax];
            P[id][ax] = old + 1e-4; const e1 = E();
            P[id][ax] = old - 1e-4; const e2 = E();
            P[id][ax] = old;
            G[id][ax] += (e1 - e2) / 2e-4;
          }
        }
      }
      let gmax = 0;
      for (let i = 0; i < n; i++) {
        for (let ax = 0; ax < dim; ax++) {
          V[i][ax] = V[i][ax] * 0.85 - dt * G[i][ax];
          if (V[i][ax] > 0.12) V[i][ax] = 0.12; else if (V[i][ax] < -0.12) V[i][ax] = -0.12;
          P[i][ax] += V[i][ax];
          const ga = Math.abs(G[i][ax]); if (ga > gmax) gmax = ga;
        }
      }
      if (it > 50 && gmax < 2e-4) break;
    }
  }

  /* ------------------------------------------------------ 2D clean-up */
  function clean2D(mol, Lpx) {
    const L = Lpx || 40, g = buildGraph(mol), n = g.n;
    if (n < 2) return mol.atoms.map(a => [a.x, a.y]);
    const rg = findRings(g), ram = ringAngleMap(rg);
    const P = mol.atoms.map(a => [a.x / L, a.y / L, 0]);
    const dist = [], rep = [];
    g.bonds.forEach(b => dist.push([b.i, b.j, 1, 5]));
    const adj = g.nbrs.map(l => l.map(e => e.j));
    for (let c = 0; c < n; c++) {
      const N = adj[c]; const d = N.length;
      if (d < 2 || d > 4) continue;
      const hasT = g.nbrs[c].some(e => e.order === 3), nd = g.nbrs[c].filter(e => e.order === 2).length;
      const sorted = N.slice().sort((p, q) => Math.atan2(P[p][1] - P[c][1], P[p][0] - P[c][0]) - Math.atan2(P[q][1] - P[c][1], P[q][0] - P[c][0]));
      for (let x = 0; x < d; x++) for (let y = x + 1; y < d; y++) {
        const i = sorted[x], j = sorted[y];
        const steps = Math.min(y - x, d - (y - x));
        let th;
        if (d === 2) th = (hasT || nd >= 2) ? 180 : 120;
        else th = 360 / d * steps;
        const key = c + '|' + Math.min(i, j) + '|' + Math.max(i, j);
        const rs = ram.get(key);
        if (rs && rs <= 8) th = (rs - 2) * 180 / rs;
        const d13 = Math.sqrt(2 - 2 * Math.cos(rad(th)));
        dist.push([i, j, d13, 2]);
      }
    }
    const D = topoDist(adj, 3);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (D[i][j] === -1 || D[i][j] >= 3) rep.push([i, j, 0.85, 0.4]);
    relax(P, dist, rep, [], { dim: 2, dt: 0.01, iters: 3000 });
    // keep centroid
    let ox = 0, oy = 0, nx = 0, ny = 0;
    mol.atoms.forEach((a, i) => { ox += a.x; oy += a.y; nx += P[i][0] * L; ny += P[i][1] * L; });
    ox = (ox - nx) / n; oy = (oy - ny) / n;
    return P.map(p => [p[0] * L + ox, p[1] * L + oy]);
  }

  /* ------------------------------------------------------ 3D embedding */
  function bondLength(e1, e2, order, arom) {
    let l = EL[e1].rc + EL[e2].rc;
    if (e1 === 'H' || e2 === 'H') l += 0.02;
    if (arom) l *= 0.92; else if (order === 2) l *= 0.88; else if (order === 3) l *= 0.79;
    return l;
  }

  function embed3D(mol, opt) {
    opt = opt || {};
    const g = buildGraph(mol), nh = mol.atoms.length;
    const rg = findRings(g), ar = aromaticity(mol, g, rg), ram = ringAngleMap(rg);
    const aromAtom = new Set(), aromBond = new Set();
    ar.rings.forEach(r => { if (r.type === 'aromatic') { r.atoms.forEach(a => aromAtom.add(a)); r.bonds.forEach(b => aromBond.add(b)); } });

    let Ls = 0; g.bonds.forEach(b => { const A = mol.atoms[b.i], B = mol.atoms[b.j]; Ls += Math.hypot(A.x - B.x, A.y - B.y); });
    const Lpx = g.bonds.length ? Ls / g.bonds.length : 40, S = 1.5 / Lpx;
    let cx = 0, cy = 0; mol.atoms.forEach(a => { cx += a.x; cy += a.y; }); cx /= nh; cy /= nh;
    const rnd = lcg(12345);

    const at = mol.atoms.map(a => ({ el: a.el, p: [(a.x - cx) * S, -(a.y - cy) * S, (rnd() - 0.5) * 0.2], src: -1, charge: a.charge || 0 }));
    g.bonds.forEach(b => { if (b.stereo) at[b.j].p[2] += (b.stereo === 1 ? 0.9 : -0.9); });

    const hyb = at.map((_, i) => {
      const nb = g.nbrs[i];
      const hasT = nb.some(e => e.order === 3), nd = nb.filter(e => e.order === 2).length;
      const el = mol.atoms[i].el, q = mol.atoms[i].charge || 0, sig = nb.length + g.hc[i];
      const flat = sig === 3 && ((el === 'C' && q === 1) || (el === 'B' && q === 0));   // carbocation / borane: planar, empty p orbital
      return (hasT || nd >= 2) ? 1 : (nd === 1 || aromAtom.has(i) || flat) ? 2 : 3;
    });

    // explicit bond list (heavy bonds + H bonds)
    const B3 = g.bonds.map((b, bi) => ({ a: b.i, b: b.j, order: b.order, arom: aromBond.has(bi) }));

    // place hydrogens
    for (let i = 0; i < nh; i++) {
      const k = g.hc[i]; if (!k) continue;
      const c = at[i].p;
      const U = g.nbrs[i].map(e => unit(sub(at[e.j].p, c)));
      let sum = [0, 0, 0]; U.forEach(u => { sum = add(sum, u); });
      let f;
      if (!U.length) f = [1, 0, 0];
      else { f = mul(sum, -1); f = len(f) < 0.25 ? [0, 0, (i % 2 ? 1 : -1)] : unit(f); }
      const hl = EL[at[i].el].rc + 0.33;
      const dirs = [];
      const h = hyb[i];
      if (!U.length && k === 4) {
        [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].forEach(d => dirs.push(unit(d)));
      } else if (U.length <= 1) {
        const alpha = rad(h === 1 ? 0 : h === 2 ? (U.length === 0 && k === 3 ? 90 : 60) : 70.5);
        let p1 = unit(cross(f, [0, 0, 1])); if (len(cross(f, [0, 0, 1])) < 0.2) p1 = unit(cross(f, [1, 0, 0]));
        let p2 = cross(f, p1);
        if (U.length === 1 && k === 1 && at[i].el === 'O') {      // carboxylic acid O–H: keep the H syn to the C=O
          const c0 = g.nbrs[i][0].j, co = g.nbrs[c0].find(e => e.j !== i && e.order === 2 && mol.atoms[e.j].el === 'O');
          if (co && mol.atoms[c0].el === 'C') {
            const v = sub(at[co.j].p, at[c0].p), pr = sub(v, mul(f, dot(v, f)));
            if (len(pr) > 0.1) { p1 = unit(pr); p2 = cross(f, p1); }
          }
        }
        const phis = k === 1 ? [0] : k === 2 ? [0, Math.PI] : [0, 2 * Math.PI / 3, 4 * Math.PI / 3];
        phis.forEach(ph => dirs.push(add(mul(f, Math.cos(alpha)), mul(add(mul(p1, Math.cos(ph)), mul(p2, Math.sin(ph))), Math.sin(alpha)))));
      } else if (U.length === 2) {
        let nrm = cross(U[0], U[1]);
        if (len(nrm) < 0.2) nrm = cross(U[0], [0, 0, 1]);
        if (len(nrm) < 0.2) nrm = cross(U[0], [1, 0, 0]);
        nrm = unit(nrm);
        const half = rad(54.75);
        if (k === 1) {
          if (h === 3) {
            const sgn = Math.abs(f[2]) > 0.05 ? Math.sign(f[2] * nrm[2]) || 1 : (i % 2 ? 1 : -1);
            dirs.push(add(mul(f, Math.cos(half)), mul(nrm, Math.sin(half) * sgn)));
          } else dirs.push(f);
        } else {
          dirs.push(add(mul(f, Math.cos(half)), mul(nrm, Math.sin(half))));
          dirs.push(add(mul(f, Math.cos(half)), mul(nrm, -Math.sin(half))));
        }
      } else {
        dirs.push(f);
      }
      dirs.slice(0, k).forEach(d => {
        const idx = at.length;
        at.push({ el: 'H', p: add(c, mul(unit(d), hl)), src: i, charge: 0 });
        B3.push({ a: i, b: idx, order: 1, arom: false });
      });
    }

    const N = at.length;
    const adj = Array.from({ length: N }, () => []);
    B3.forEach(b => { adj[b.a].push(b.b); adj[b.b].push(b.a); });
    const P = at.map(a => a.p);
    const bl = {};
    B3.forEach(b => { bl[b.a + '|' + b.b] = bl[b.b + '|' + b.a] = bondLength(at[b.a].el, at[b.b].el, b.order, b.arom); });

    const dist = [], rep = [], chir = [];
    B3.forEach(b => dist.push([b.a, b.b, bl[b.a + '|' + b.b], 6]));

    const hybOf = i => i < nh ? hyb[i] : 3;
    const ringAng = rs => (rs - 2) * 180 / rs;
    const angleTarget = (c, i, j) => {
      const h = hybOf(c);
      let th = h === 1 ? 180 : h === 2 ? 120 : 109.5;
      if (c < nh) {
        const rsOf = (x, y) => (x < nh && y < nh) ? ram.get(c + '|' + Math.min(x, y) + '|' + Math.max(x, y)) : undefined;
        const rs = rsOf(i, j);
        if (h === 2 && adj[c].length === 3) {
          // sp2 atom: the three angles must sum to 360°. Ring angles (3-7 membered) are fixed by the ring,
          // the remaining angle(s) share what is left (this also keeps fused 5/7 systems such as azulene flat).
          const nb = adj[c]; let sum = 0, nRing = 0;
          for (let x = 0; x < 3; x++) for (let y = x + 1; y < 3; y++) { const r2 = rsOf(nb[x], nb[y]); if (r2 && r2 <= 7) { sum += ringAng(r2); nRing++; } }
          if (rs && rs <= 7) th = ringAng(rs);
          else if (!rs && nRing > 0 && nRing < 3) th = (360 - sum) / (3 - nRing);
        } else if (rs && rs <= 5) th = ringAng(rs);
      }
      return th;
    };
    for (let c = 0; c < N; c++) {
      const nb = adj[c]; if (nb.length < 2 || nb.length > 4) continue;
      for (let x = 0; x < nb.length; x++) for (let y = x + 1; y < nb.length; y++) {
        const i = nb[x], j = nb[y];
        const l1 = bl[c + '|' + i], l2 = bl[c + '|' + j];
        const th = rad(angleTarget(c, i, j));
        dist.push([i, j, Math.sqrt(l1 * l1 + l2 * l2 - 2 * l1 * l2 * Math.cos(th)), 3]);
      }
    }

    // planar 1-4 restraints (double bonds, aromatic bonds) → preserves E/Z and flatness
    const dihedral = (a, b, c, d) => {
      const b1 = sub(b, a), b2 = sub(c, b), b3 = sub(d, c);
      const n1 = cross(b1, b2), n2 = cross(b2, b3);
      return Math.atan2(dot(cross(n1, n2), unit(b2)), dot(n1, n2));
    };
    B3.forEach(bd => {
      const b = bd.a, c = bd.b;
      const planar = (bd.order === 2 || bd.arom) && b < nh && c < nh && hyb[b] === 2 && hyb[c] === 2;
      if (!planar) return;
      adj[b].forEach(a => { if (a === c) return; adj[c].forEach(d => {
        if (d === b || d === a) return;
        const cis = Math.abs(dihedral(P[a], P[b], P[c], P[d])) < Math.PI / 2;
        const lab = bl[a + '|' + b], lbc = bl[b + '|' + c], lcd = bl[c + '|' + d];
        const t1 = rad(angleTarget(b, a, c)), t2 = rad(angleTarget(c, b, d));
        const ax = lab * Math.cos(t1), ay = lab * Math.sin(t1);
        const dx = lbc - lcd * Math.cos(t2), dy = (cis ? 1 : -1) * lcd * Math.sin(t2);
        dist.push([a, d, Math.hypot(ax - dx, ay - dy), 2]);
      }); });
    });

    // chirality restraints at sp3 centres with 4 neighbours
    for (let c = 0; c < nh; c++) {
      if (hyb[c] !== 3 || adj[c].length !== 4) continue;
      const nHs = adj[c].filter(x => x >= nh || at[x].el === 'H').length;
      if (nHs >= 2) continue;
      const nb = adj[c];
      let done = false;
      for (let x = 0; x < 4 && !done; x++) for (let y = x + 1; y < 4 && !done; y++) for (let z = y + 1; z < 4 && !done; z++) {
        const T = dot(unit(sub(P[nb[x]], P[c])), cross(unit(sub(P[nb[y]], P[c])), unit(sub(P[nb[z]], P[c]))));
        if (Math.abs(T) > 0.1) { chir.push({ c, a: nb[x], b: nb[y], d: nb[z], t: Math.sign(T) * 0.77, k: 4 }); done = true; }
      }
    }

    const D = topoDist(adj, 4);
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      const d = D[i][j];
      if (d !== -1 && d < 3) continue;
      const f = d === 3 ? 0.62 : 0.8;
      rep.push([i, j, (EL[at[i].el].vdw + EL[at[j].el].vdw) * f, 0.5]);
    }

    relax(P, dist, rep, chir, { dim: 3, dt: 0.01, iters: opt.iters || 3000 });

    // centre
    const cen = [0, 0, 0]; P.forEach(p => { cen[0] += p[0]; cen[1] += p[1]; cen[2] += p[2]; });
    cen[0] /= N; cen[1] /= N; cen[2] /= N;
    P.forEach(p => { p[0] -= cen[0]; p[1] -= cen[1]; p[2] -= cen[2]; });

    return {
      atoms: at.map((a, i) => ({ el: a.el, x: P[i][0], y: P[i][1], z: P[i][2], heavy: i < nh, charge: a.charge })),
      bonds: B3.map(b => ({ a: b.a, b: b.b, order: b.order, arom: b.arom })),
      ar
    };
  }

  /* ------------------------------------------- dipole moment (rough estimate) */
  // Pauling electronegativities
  const EN = { H: 2.20, B: 2.04, C: 2.55, N: 3.04, O: 3.44, F: 3.98, Si: 1.90, P: 2.19, S: 2.58, Cl: 3.16, Br: 2.96, I: 2.66, Se: 2.55, Li: 0.98, Na: 0.93, K: 0.82, Mg: 1.31, Ca: 1.00, Al: 1.61 };
  // empirical bond moments in Debye (values that reproduce small reference molecules).
  // key = "El1|El2|order" with the two symbols sorted; value = [atom at the negative end, moment]
  const BM = {
    'F|H|1': ['F', 1.82], 'Cl|H|1': ['Cl', 1.08], 'Br|H|1': ['Br', 0.82], 'H|I|1': ['I', 0.44],
    'H|O|1': ['O', 1.51], 'H|N|1': ['N', 1.31], 'H|S|1': ['S', 0.65],
    'C|F|1': ['F', 1.41], 'C|Cl|1': ['Cl', 1.46], 'C|Br|1': ['Br', 1.38], 'C|I|1': ['I', 1.19],
    'C|O|1': ['O', 0.74], 'C|O|2': ['O', 2.00],
    'C|N|1': ['N', 0.22], 'C|N|2': ['N', 0.90], 'C|N|3': ['N', 3.20],
    'C|S|1': ['S', 0.90], 'C|S|2': ['S', 1.80], 'O|S|2': ['O', 2.70], 'O|P|2': ['O', 3.50]
  };
  // halogen on an sp2 carbon (aryl / vinyl halides) is less polar than on sp3
  const BM_SP2 = { 'C|F|1': ['F', 1.30], 'C|Cl|1': ['Cl', 1.40], 'C|Br|1': ['Br', 1.30], 'C|I|1': ['I', 1.20] };
  const CH_M = [0.40, 0.30, 1.00];                       // C–H moment (C is the negative end) for sp3 / sp2 / sp carbon
  const CC_M = { '0|1': 0.40, '0|2': 0.50, '1|2': 0.30 };  // C–C between different hybrids (more s-character carbon is negative)
  const PI_DONOR = { N: 0.98, O: 0.45, S: 2.08 };                 // π-donation of 5-ring heteroatoms towards the ring (Debye)
  const DEBYE = 4.803;

  function dipole(mol, s3) {
    const nh = mol.atoms.length;
    if (!nh) return { supported: false, reason: 'Empty drawing.' };
    const chg = mol.atoms.map(a => a.charge || 0);
    if (chg.reduce((s, q) => s + q, 0) !== 0) return { supported: false, reason: 'The molecule carries a net charge, so its dipole depends on the chosen origin and is not defined.' };
    const g = buildGraph(mol);
    for (let i = 0; i < nh; i++) if (chg[i]) {
      const el = mol.atoms[i].el;
      const ok = chg[i] === 1
        ? (['N', 'S', 'P'].includes(el) && g.nbrs[i].some(e => chg[e.j] === -1 && mol.atoms[e.j].el === 'O'))
        : (el === 'O' && g.nbrs[i].some(e => chg[e.j] === 1 && ['N', 'S', 'P'].includes(mol.atoms[e.j].el)));
      if (!ok) return { supported: false, reason: 'Only nitro / N-oxide / S–O type charge-separated groups are supported; other zwitterions depend on charge distribution this model cannot estimate.' };
    }
    s3 = s3 || embed3D(mol);
    const A = s3.atoms, P3 = i => [A[i].x, A[i].y, A[i].z];
    const arom = new Set(); s3.bonds.forEach(b => { if (b.arom) { arom.add(b.a); arom.add(b.b); } });
    // hybridisation rank: 0 = sp3, 1 = sp2, 2 = sp
    const hyb = mol.atoms.map((_, i) => {
      const nb = g.nbrs[i], t = nb.some(e => e.order === 3), nd = nb.filter(e => e.order === 2).length;
      return (t || nd >= 2) ? 2 : (nd === 1 || arom.has(i)) ? 1 : 0;
    });
    const hy = i => i < nh ? hyb[i] : 0;
    const q = i => i < nh ? chg[i] : 0;

    function moment(ia, ib, ord, isArom, d) {
      const ea = A[ia].el, eb = A[ib].el;
      if ((q(ia) > 0 && eb === 'O') || (q(ib) > 0 && ea === 'O')) {            // N+–O−, N+=O, S+–O− …
        const xi = q(ia) > 0 ? ia : ib, oi = xi === ia ? ib : ia;
        let m = 3.3;
        if (A[xi].el === 'N') {
          const term = g.nbrs[xi].filter(e => mol.atoms[e.j].el === 'O' && g.nbrs[e.j].length === 1).length;
          const cn = g.nbrs[xi].find(e => mol.atoms[e.j].el === 'C');
          m = term >= 2 ? ((cn && hy(cn.j) >= 1) ? 3.7 : 2.8) : 4.3;         // nitro (conjugated or not) vs N-oxide
        }
        return { m, neg: oi };
      }
      const key = [ea, eb].sort().join('|');
      if (key === 'C|H') { const ci = ea === 'C' ? ia : ib; return { m: CH_M[hy(ci)], neg: ci }; }
      if (key === 'C|C') {
        const ra = hy(ia), rb = hy(ib);
        if (ra === rb) return null;
        const lo = Math.min(ra, rb), hi = Math.max(ra, rb);
        return { m: CC_M[lo + '|' + hi], neg: ra > rb ? ia : ib };
      }
      if (key === 'C|N' && isArom && g.nbrs[ea === 'N' ? ia : ib].length === 2 && g.hc[ea === 'N' ? ia : ib] === 0) return { m: 1.9, neg: ea === 'N' ? ia : ib };
      const ci = ea === 'C' ? ia : eb === 'C' ? ib : -1;
      const k1 = key + '|' + ord;
      let ent = (ci >= 0 && hy(ci) >= 1 && BM_SP2[k1]) || BM[k1];
      if (ent) {
        let m = ent[1];
        if (ci >= 0 && ['F', 'Cl', 'Br', 'I'].includes(ent[0])) {              // geminal halogens weaken each other
          const nX = g.nbrs[ci].filter(e => ['F', 'Cl', 'Br', 'I'].includes(mol.atoms[e.j].el)).length;
          m *= Math.max(0.3, 1 - 0.28 * (nX - 1));
          if (hy(ci) >= 1) g.nbrs[ci].forEach(e => {                           // halogens on the other end of a C=C
            if (e.order !== 2) return;
            const nY = g.nbrs[e.j].filter(f => ['F', 'Cl', 'Br', 'I'].includes(mol.atoms[f.j].el)).length;
            m *= Math.max(0.4, 1 - 0.24 * nY);
          });
        }
        if (key === 'O|S' && ord === 2) {                                      // sulfoxide vs sulfone
          const si = ea === 'S' ? ia : ib;
          m = g.nbrs[si].filter(e => e.order === 2 && mol.atoms[e.j].el === 'O').length >= 2 ? 3.6 : 2.7;
        }
        return { m, neg: ea === ent[0] ? ia : ib };
      }
      if (ea === eb) return null;
      const dEN = (EN[ea] || 2) - (EN[eb] || 2), ad = Math.abs(dEN);                // Hannay–Smyth fallback
      return { m: DEBYE * d * (0.16 * ad + 0.035 * ad * ad), neg: dEN > 0 ? ia : ib };
    }

    let D = [0, 0, 0];
    const parts = [];
    s3.bonds.forEach(b => {
      const pa = P3(b.a), pb = P3(b.b), v = sub(pb, pa), d = len(v);
      if (d < 1e-6) return;
      const r = moment(b.a, b.b, b.order, !!b.arom, d);
      if (!r || !r.m) return;
      const dir = unit(r.neg === b.b ? v : mul(v, -1));                         // points to the δ− end
      D = add(D, mul(dir, r.m));
      parts.push({ a: b.a, b: b.b, m: r.m });
    });
    // lone-pair donation of pyrrole / furan / thiophene type atoms into a 5-membered aromatic ring
    s3.ar.rings.forEach(r => {
      if (r.type !== 'aromatic' || r.size !== 5) return;
      let c = [0, 0, 0]; r.atoms.forEach(i => { c = add(c, P3(i)); }); c = mul(c, 1 / 5);
      r.atoms.forEach(i => {
        const k = PI_DONOR[mol.atoms[i].el];
        if (!k || g.nbrs[i].some(e => e.order > 1)) return;
        D = add(D, mul(unit(sub(c, P3(i))), k));
      });
    });
    const mu = len(D);
    const cls = mu < 0.1 ? 'non-polar' : mu < 1 ? 'weakly polar' : mu < 2.5 ? 'polar' : 'strongly polar';
    return { supported: true, mu, vec: D, cls };
  }

  /* ------------------------------------------------ E/Z (CIP priorities) */
  const ZNUM = { H: 1, Li: 3, B: 5, C: 6, N: 7, O: 8, F: 9, Na: 11, Mg: 12, Al: 13, Si: 14, P: 15, S: 16, Cl: 17, K: 19, Ca: 20, Se: 34, Br: 35, I: 53 };

  // CIP ranking helpers shared by E/Z and R/S. cmpSub(p, q, depth) > 0 when substituent p outranks q.
  // A substituent node is {z, i, kind:'real'|'dup'|'H'|'lp', anc:[atoms from the stereo centre]}.
  function cipTools(mol, g) {
    const Zof = i => ZNUM[mol.atoms[i].el] || 0;
    const kids = n => {
      if (n.kind !== 'real') return [];
      const out = [], par = n.anc[n.anc.length - 1];
      g.nbrs[n.i].forEach(e => {
        if (e.j !== par) {
          if (n.anc.includes(e.j)) out.push({ z: Zof(e.j), i: -1, kind: 'dup' });            // ring closure
          else out.push({ z: Zof(e.j), i: e.j, kind: 'real', anc: n.anc.concat(n.i) });
        }
        for (let k = 1; k < e.order; k++) out.push({ z: Zof(e.j), i: -1, kind: 'dup' });       // duplicate for a multiple bond
      });
      for (let k = 0; k < g.hc[n.i]; k++) out.push({ z: 1, i: -1, kind: 'H' });
      return out;
    };
    const sortKids = (n, depth) => kids(n).sort((a, b) => depth > 0 ? cmpSub(b, a, depth) : b.z - a.z);
    function cmpSub(p, q, depth) {                      // > 0 when p outranks q
      if (p.z !== q.z) return p.z - q.z;
      let LP = [p], LQ = [q];
      for (let s = 0; s < 8 && (LP.length || LQ.length); s++) {
        const GP = LP.map(n => sortKids(n, depth - 1)), GQ = LQ.map(n => sortKids(n, depth - 1));
        const m = Math.max(GP.length, GQ.length);
        for (let k = 0; k < m; k++) {
          const a = GP[k] || [], b = GQ[k] || [], l = Math.max(a.length, b.length);
          for (let t = 0; t < l; t++) { const za = a[t] ? a[t].z : 0, zb = b[t] ? b[t].z : 0; if (za !== zb) return za - zb; }
        }
        LP = [].concat(...GP); LQ = [].concat(...GQ);
      }
      return 0;
    }
    return { Zof, cmpSub };
  }

  // E/Z label for every acyclic C=C, C=N or N=N bond whose two ends each carry two different substituents.
  // Priorities follow the CIP rules (atomic number, then sphere by sphere outwards, with duplicated atoms for
  // multiple bonds and ring closures). Geometry is read from the 2D drawing.
  function ezStereo(mol, g, rg) {
    g = g || buildGraph(mol); rg = rg || findRings(g);
    const { Zof, cmpSub } = cipTools(mol, g);
    const bonds = []; let doubles = 0;
    g.bonds.forEach((bd, bi) => {
      if (bd.order !== 2 || rg.bondInRing[bi]) return;
      const A = bd.i, B = bd.j, ok = ['C', 'N'].includes(mol.atoms[A].el) && ['C', 'N'].includes(mol.atoms[B].el);
      if (!ok) return;
      doubles++;
      if ([A, B].some(x => g.nbrs[x].some(e => e.bi !== bi && e.order >= 2))) return;           // allene / conjugated cumulene / alkyne end
      const ends = [[A, B], [B, A]].map(([X, O]) => {
        const subs = [];
        g.nbrs[X].forEach(e => { if (e.j !== O) subs.push({ z: Zof(e.j), i: e.j, kind: 'real', anc: [X] }); });
        for (let k = 0; k < g.hc[X]; k++) subs.push({ z: 1, i: -1, kind: 'H' });
        if (mol.atoms[X].el === 'N' && subs.length === 1) subs.push({ z: 0, i: -1, kind: 'lp' });   // lone pair ranks last
        return subs;
      });
      if (ends.some(s => s.length !== 2)) return;
      const top = ends.map(s => { const c = cmpSub(s[0], s[1], 3); return c === 0 ? null : (c > 0 ? s[0] : s[1]); });
      if (!top[0] || !top[1] || top[0].i < 0 || top[1].i < 0) return;                           // identical groups, or highest group not drawn
      const pa = mol.atoms[A], pb = mol.atoms[B], ax = pb.x - pa.x, ay = pb.y - pa.y;
      const side = i => Math.sign(ax * (mol.atoms[i].y - pa.y) - ay * (mol.atoms[i].x - pa.x));
      const s1 = side(top[0].i), s2 = side(top[1].i);
      if (!s1 || !s2) return;
      const heavy = ends.map(s => s.filter(x => x.z > 1));
      let cistrans = null;
      if (heavy[0].length === 1 && heavy[1].length === 1 && heavy[0][0].i >= 0 && heavy[1][0].i >= 0)
        cistrans = side(heavy[0][0].i) === side(heavy[1][0].i) ? 'cis' : 'trans';
      bonds.push({ mi: bd.mi, a: A, b: B, label: s1 === s2 ? 'Z' : 'E', cistrans, ha: top[0].i, hb: top[1].i });
    });
    return { bonds, doubles };
  }

  // Fischer projection: a centre drawn as a cross (neighbours along the four axis directions, no wedge / dash of
  // its own). Horizontal bonds point towards the viewer, vertical bonds away. 3 drawn neighbours + implicit H: the H
  // sits on the empty arm. Returns the 4 vectors in rank order, or null when the centre is not a Fischer cross.
  function fischerVecs(mol, g, c, subs) {
    const A = mol.atoms[c], arms = [[1, 0, 0.7], [0, 1, -0.7], [-1, 0, 0.7], [0, -1, -0.7]], used = new Set();
    const arm = new Map();
    for (const s of subs) {
      if (s.i < 0) continue;
      const ang = Math.atan2(-(mol.atoms[s.i].y - A.y), mol.atoms[s.i].x - A.x) * 180 / Math.PI;
      let best = -1, bd = 1e9;
      for (let k = 0; k < 4; k++) { let d = Math.abs(((ang - 90 * k) % 360 + 540) % 360 - 180); if (d < bd) { bd = d; best = k; } }
      if (bd > 28 || used.has(best)) return null;
      used.add(best); arm.set(s.i, best);
    }
    if (used.size < 3) return null;
    const free = [0, 1, 2, 3].filter(k => !used.has(k));
    return subs.map(s => s.i >= 0 ? arms[arm.get(s.i)] : arms[free[0]]);
  }

  /* --------------------------------------- R/S (CIP) and optical activity */
  // Finds every sp3 stereocentre (C, Si, N+, P+ with four different substituents under the CIP rules) and assigns R/S
  // from the wedge / dash bonds that START at the centre (narrow end on the centre), or, when a centre is drawn as a
  // Fischer cross without wedges, from the Fischer convention. Then decides whether the molecule
  // is chiral, using an internal-mirror (meso) test. Axial and planar chirality are not assessed.
  function rsStereo(mol, g, rg) {
    g = g || buildGraph(mol); rg = rg || findRings(g);
    const { Zof, cmpSub } = cipTools(mol, g), n = mol.atoms.length;
    const okEl = (el, q) => ((el === 'C' || el === 'Si') && q === 0) || ((el === 'N' || el === 'P') && q === 1);
    const centres = [];
    for (let c = 0; c < n; c++) {
      const A = mol.atoms[c], nb = g.nbrs[c];
      if (!okEl(A.el, A.charge || 0) || nb.some(e => e.order > 1) || nb.length + g.hc[c] !== 4 || g.hc[c] > 1) continue;
      const subs = nb.map(e => ({ z: Zof(e.j), i: e.j, kind: 'real', anc: [c] }));
      if (g.hc[c] === 1) subs.push({ z: 1, i: -1, kind: 'H' });
      subs.sort((p, q) => cmpSub(q, p, 4));
      let distinct = true;
      for (let k = 0; k < 3; k++) if (cmpSub(subs[k], subs[k + 1], 4) === 0) distinct = false;
      if (!distinct) continue;
      // geometry: x right, y up (canvas y is flipped), z towards the viewer set by wedge (+) / dash (−)
      const zOf = new Map(); let anyZ = false;
      g.bonds.forEach(b => { if (b.stereo && b.i === c) { zOf.set(b.j, b.stereo === 1 ? 1 : -1); anyZ = true; } });
      let V, fischer = false;
      const fz = anyZ ? null : fischerVecs(mol, g, c, subs);
      if (fz) { V = fz; fischer = true; }
      else {
        const vec = s => {
          const dx = mol.atoms[s.i].x - A.x, dy = -(mol.atoms[s.i].y - A.y), l = Math.hypot(dx, dy) || 1;
          return [dx / l, dy / l, (zOf.get(s.i) || 0) * 0.9];
        };
        V = subs.map(s => s.i >= 0 ? vec(s) : null);
        if (V[3] === null) { let sm = [0, 0, 0]; for (let k = 0; k < 3; k++) sm = add(sm, unit(V[k])); V[3] = mul(sm, -1); }
      }
      let label = null;
      if (anyZ || fischer) {
        const d = [sub(V[0], V[3]), sub(V[1], V[3]), sub(V[2], V[3])], det = dot(d[0], cross(d[1], d[2]));
        if (Math.abs(det) > 1e-3) label = det < 0 ? 'R' : 'S';
      }
      centres.push({ atom: c, label, order: subs.map(s => s.i), hasWedge: anyZ, fischer });
    }
    const out = { centres, chiral: null, meso: false, verdict: '' };
    if (!centres.length) { out.chiral = false; out.verdict = 'none'; return out; }
    if (centres.some(x => !x.label)) { out.verdict = 'undefined'; return out; }
    if (centres.length === 1) { out.chiral = true; out.verdict = 'chiral'; return out; }
    const mir = mirrorSymmetric(mol, g, centres, ezStereo(mol, g, rg).bonds);
    if (mir === null) { out.verdict = 'unknown'; return out; }
    out.meso = mir; out.chiral = !mir; out.verdict = mir ? 'meso' : 'chiral';
    return out;
  }

  // true when some graph automorphism maps every stereocentre onto one of opposite label (the molecule equals its
  // own mirror image, e.g. meso-tartaric acid); null when the search was cut off.
  function mirrorSymmetric(mol, g, centres, ez) {
    const n = mol.atoms.length, lab = new Map(centres.map(c => [c.atom, c.label]));
    const flip = l => l === 'R' ? 'S' : 'R';
    const ezOf = new Map(); ez.forEach(b => { ezOf.set(b.a + '|' + b.b, b.label); ezOf.set(b.b + '|' + b.a, b.label); });
    let col = mol.atoms.map((a, i) => a.el + '|' + (a.charge || 0) + '|' + g.hc[i] + '|' + g.nbrs[i].length);
    for (let r = 0; r < 6; r++) {
      const nx = col.map((c, i) => c + '(' + g.nbrs[i].map(e => e.order + col[e.j]).sort().join(',') + ')');
      const ids = new Map(); col = nx.map(s => { if (!ids.has(s)) ids.set(s, ids.size); return String(ids.get(s)); });
    }
    const order = [], seen = new Array(n).fill(false);
    for (let s = 0; s < n; s++) {
      if (seen[s]) continue; seen[s] = true; const q = [s];
      for (let h = 0; h < q.length; h++) { order.push(q[h]); g.nbrs[q[h]].forEach(e => { if (!seen[e.j]) { seen[e.j] = true; q.push(e.j); } }); }
    }
    const map = new Array(n).fill(-1), used = new Array(n).fill(false), orderOf = (a, b) => { const e = g.nbrs[a].find(x => x.j === b); return e ? e.order : 0; };
    let steps = 0, cut = false;
    function fin() {
      for (const [k, l] of ezOf) { const [a, b] = k.split('|').map(Number); if (ezOf.get(map[a] + '|' + map[b]) !== l) return false; }
      return true;
    }
    function rec(k) {
      if (cut) return false;
      if (k === n) return fin();
      if (++steps > 200000) { cut = true; return false; }
      const u = order[k];
      for (let v = 0; v < n; v++) {
        if (used[v] || col[v] !== col[u]) continue;
        const lu = lab.get(u), lv = lab.get(v);
        if ((lu === undefined) !== (lv === undefined) || (lu !== undefined && lv !== flip(lu))) continue;
        let ok = true;
        for (const e of g.nbrs[u]) if (map[e.j] >= 0 && orderOf(v, map[e.j]) !== e.order) { ok = false; break; }
        if (!ok) continue;
        map[u] = v; used[v] = true;
        if (rec(k + 1)) return true;
        map[u] = -1; used[v] = false;
      }
      return false;
    }
    const res = rec(0);
    return cut ? null : res;
  }

  /* ------------------------------------ Hansen solubility parameters (estimate) */
  // Hoftyzer–Van Krevelen group contributions: [Fd, Fp, Eh, V]
  //   Fd, Fp in J^0.5 cm^1.5 / mol, Eh in J / mol, V (molar volume increment) in cm^3 / mol.
  //   dD = ΣFd / V     dP = sqrt(ΣFp²) / V     dH = sqrt(ΣEh / V)      (all in MPa^0.5)
  const HVK = {
    'CH3': [420, 0, 0, 33.5], 'CH2': [270, 0, 0, 16.1], 'CH': [80, 0, 0, -1.0], 'C': [-70, 0, 0, -19.2],
    '=CH2': [400, 0, 0, 28.5], '=CH': [200, 0, 0, 13.5], '=C<': [70, 0, 0, -5.5],
    'C6H5': [1430, 110, 0, 71.4], 'C6H4': [1270, 110, 0, 52.4],
    // unsubstituted benzene: fitted so the ring reproduces benzene itself (dD 18.4, dP 0, dH 2.0, V 89.4)
    'C6H6': [1645, 0, 360, 89.4],
    // tri- / tetra-substituted rings: extrapolated from the C6H5 → C6H4 step (about −160 in Fd and −19 cm³/mol per H removed)
    'C6H3': [1110, 110, 0, 33.4], 'C6H2': [950, 110, 0, 14.4],
    'F': [220, 0, 0, 18.0], 'Cl': [450, 550, 400, 24.0], 'Br': [550, 0, 0, 30.0],
    'CN': [430, 1100, 2500, 24.0], 'OH': [210, 500, 20000, 10.0], 'O': [100, 400, 3000, 3.8],
    'COOH': [530, 420, 10000, 28.5], 'CO': [290, 770, 2000, 10.8], 'COO': [390, 490, 7000, 18.0],
    'NH2': [280, 0, 8400, 19.2], 'NH': [160, 210, 3100, 4.5], 'NO2': [500, 1070, 1500, 32.0]
  };
  const HVK_PARTIAL = { 'F': 1, 'Br': 1 };   // polar / H-bond terms are not tabulated for these groups
  // common solvents: [name, dD, dP, dH] in MPa^0.5 (Hansen)
  const SOLV = [
    ['Water', 15.5, 16.0, 42.3, 18.0], ['Methanol', 14.7, 12.3, 22.3, 40.7], ['Ethanol', 15.8, 8.8, 19.4, 58.5], ['1-Propanol', 16.0, 6.8, 17.4, 75.2],
    ['2-Propanol', 15.8, 6.1, 16.4, 76.8], ['1-Butanol', 16.0, 5.7, 15.8, 91.5], ['Ethylene glycol', 17.0, 11.0, 26.0, 55.8], ['Glycerol', 17.4, 12.1, 29.3, 73.3],
    ['Acetone', 15.5, 10.4, 7.0, 74.0], ['2-Butanone (MEK)', 16.0, 9.0, 5.1, 90.1], ['Ethyl acetate', 15.8, 5.3, 7.2, 98.5], ['Tetrahydrofuran', 16.8, 5.7, 8.0, 81.7],
    ['Diethyl ether', 14.5, 2.9, 5.1, 104.8], ['1,4-Dioxane', 17.5, 1.8, 9.0, 85.7], ['Dichloromethane', 18.2, 6.3, 6.1, 63.9], ['Chloroform', 17.8, 3.1, 5.7, 80.7],
    ['Carbon tetrachloride', 17.8, 0.0, 0.6, 97.1], ['Toluene', 18.0, 1.4, 2.0, 106.8], ['Benzene', 18.4, 0.0, 2.0, 89.4], ['n-Hexane', 14.9, 0.0, 0.0, 131.6],
    ['n-Heptane', 15.3, 0.0, 0.0, 147.4], ['Cyclohexane', 16.8, 0.0, 0.2, 108.7], ['Acetonitrile', 15.3, 18.0, 6.1, 52.6], ['DMF', 17.4, 13.7, 11.3, 77.0],
    ['DMSO', 18.4, 16.4, 10.2, 71.3], ['NMP', 18.0, 12.3, 7.2, 96.5], ['Acetic acid', 14.5, 8.0, 13.5, 57.1], ['Pyridine', 19.0, 8.8, 5.9, 80.9], ['Nitromethane', 15.8, 18.8, 5.1, 54.3]
  ];
  const RT = 8.314 * 298.15;   // J/mol at 25 °C

  const SP1 = {
    'CH3': [-0.9714, -1.6448, -0.7813],
    'CH2': [-0.0269, -0.3045, -0.4119],
    'CH': [0.645, 0.6491, -0.2018],
    'C': [1.2686, 2.0838, 0.0866],
    'CH2=CH': [-1.0585, -2.0035, -1.2985],
    'CH=CH': [0.0048, -0.2984, -0.04],
    'CH2=C': [-0.4829, -0.7794, -0.826],
    'CH=C': [0.5372, -0.9024, -1.8872],
    'C=C': [0.3592, 1.0526, -15.4659],
    'CH2=C=CH': [-1.6518, null, -0.998],
    'CH#C': [0.232, -1.3294, 1.0736],
    'C#C': [-0.2028, -0.7598, -1.1083],
    'ACH': [0.1105, -0.5303, -0.4305],
    'AC': [0.8446, 0.6187, 0.0084],
    'ACCH3': [0.2174, -0.5705, -1.1473],
    'ACCH2': [0.6933, 0.6517, -0.1375],
    'CH3CO': [-0.3551, 2.3192, -1.3078],
    'CH2CO': [0.6527, 3.7328, -0.5344],
    'CHO': [-0.403, 3.4734, 0.1687],
    'COOH': [-0.291, 0.9042, 3.7391],
    'CH3COO': [-0.5401, -0.397, 1.5826],
    'CH2COO': [0.2913, 3.6462, 1.2523],
    'HCOO': [null, 1.9308, 2.1202],
    'COO': [0.2039, 3.4637, 1.1389],
    'OH': [-0.3462, 1.1404, 7.1908],
    'ACOH': [0.5288, 1.101, 6.958],
    'CH3O': [-0.5828, 0.1764, 0.146],
    'CH2O': [0.031, 0.8826, -0.1528],
    'CH-O': [0.8833, 1.6853, 0.447],
    'C2H5O2': [-0.1249, 3.6422, 8.3579],
    'CH2O-ring': [0.2753, 0.1994, -0.161],
    'CH2NH2': [-0.5828, 1.4084, 2.592],
    'CHNH2': [0.0112, -1.1989, 0.3818],
    'CH3NH': [null, 0.6777, 5.6646],
    'CH2NH': [0.8116, 0.9412, 1.34],
    'CH3N': [0.8769, 1.2046, 1.6062],
    'CH2N': [1.4681, 2.8345, 1.2505],
    'ACNH2': [1.6987, 1.6761, 4.5274],
    'CONH2': [-0.0689, 6.0694, 5.228],
    'CON(CH3)2': [0.4482, 5.7899, 3.002],
    'CH2SH': [1.2797, -0.8223, 4.4646],
    'CH3S': [null, 0.4944, -1.4861],
    'CH2S': [1.0595, 0.753, -0.2287],
    'I': [0.7797, 0.6777, 0.2646],
    'Br': [0.5717, 0.6997, -1.0722],
    'CH2Cl': [0.2623, 0.597, -0.5364],
    'CHCl': [0.4462, 2.806, -1.4125],
    'CCl': [2.7576, 2.0406, 0.1101],
    'CHCl2': [1.1797, 1.8361, -3.2861],
    'CCl2': [0.3653, 0.1696, -1.4334],
    'CCl3': [null, 1.2777, -2.6354],
    'ACCl': [0.8475, -0.0339, -0.784],
    'ACF': [0.117, 0.1856, -0.7182],
    'Cl-vinyl': [0.2289, 2.3444, 3.8893],
    'CF3': [-0.2293, -1.9735, -1.4665],
    'CH2NO2': [null, 6.8944, -1.2861],
    'CHNO2': [null, 8.0347, -2.3167],
    'ACNO2': [1.4195, 4.4838, -0.7167],
    'CH2CN': [-0.3392, 6.5341, -0.8892],
    'CF2': [-0.9729, null, null],
    'CF': [0.1707, null, null],
    'F': [-0.7069, null, null],
    'CH2=C=C': [-0.2804, null, -1.9167],
    'O': [0.0472, 3.3432, 0.0256],
    'Cl': [0.2256, 1.8711, -0.3295],
    'C=N': [-0.3074, -0.0012, -5.3956],
    'CH=N': [0.9672, 1.9728, 0.7668],
    'NH': [null, 0.0103, 2.2086],
    'CN': [0.0861, 6.5331, -0.6849],
    'NCO': [-0.1306, 1.6102, 4.0461],
    'SH': [1.0427, 1.9813, 4.8181],
    'S': [1.4899, 9.2072, -0.625],
    'SO2': [1.5502, 11.1758, 0.1055],
    'C=S': [0.7747, 0.0683, 3.408],
    'C=O': [-0.4343, 0.7905, 1.8147],
    'N': [1.5438, 2.578, 1.1189]
  };
  const SP2 = {
    'iPr': [0.046, 0.0019, 0.3149],
    'tBu': [-0.0738, 1.1881, -0.2966],
    'ring5': [-0.6681, -2.343, -0.3079],
    'ring6': [-0.3874, -3.6432, null],
    'diene': [-0.1355, -3.5085, -1.0795],
    'CH3-C=': [-0.0785, 0.3316, 0.3875],
    'CH2-C=': [-0.3236, -2.3179, -0.5836],
    'CHorC-C=': [-0.2798, null, -1.1164],
    'cycstring': [-0.1945, null, null],
    'CH3COCH2': [-0.0451, -0.3383, -0.4083],
    'Ccyc=O': [-0.2981, 0.4497, -0.4794],
    'ACCOOH': [-0.2293, -0.6349, -0.903],
    'CHorC-COOH': [null, -0.2187, 1.146],
    'CH3COO-CHorC': [-0.522, -0.0652, 0.3085],
    'COCH2COO': [null, -2.3792, 0.8412],
    'anhydride': [-0.2707, -1.0562, 1.6335],
    'ACHO': [0.3772, -1.811, -1.0096],
    'CHOH': [0.1123, 0.2564, -0.1928],
    'C(OH)': [-0.068, 0.1075, 1.2931],
    'diol': [null, 0.6419, 0.387],
    'aminoalcohol': [-0.0809, 0.5683, -0.6326],
    'Ccyc-OH': [-0.0876, -3.522, 0.5914],
    'C-O-C=C': [0.2063, 0.608, 1.1344],
    'AC-O-C': [0.2568, 0.8153, 0.6092],
    'Ncyc': [0.2218, -2.2018, -0.0452],
    'Scyc': [0.4892, 0.304, 0.2297],
    'ACBr': [0.1234, -0.4495, 0.3397],
    'C=C-Br': [-0.4059, -0.0024, -1.1304],
    'ring3': [0.02, 1.8288, -0.8073],
    'ACCOO': [-0.1847, 0.4059, -0.1921],
    'naph': [-0.3751, -1.298, 0.6844],
    'OcycCcycO': [0.2468, 2.7501, 0.122],
    'AC-O-AC': [-0.5646, -3.4329, 2.083],
    'pyr': [0.7002, 0.0691, -2.7661],
    'lactam': [0.2956, 2.8958, 1.3125],
    'acetal': [0.0839, 0.3451, 0.3767],
    'dicarbonyl': [-0.4862, -0.4888, 1.2482]
  };
  const SP1L = {
    'CH3': [-0.72412, 0.29901],
    'CH2': [-0.1403, -0.1161],
    'CH': [0.58978, 0.1386],
    'CH2=CH': [-0.29774, 1.35521],
    'CH=CH': [-0.22864, 0.48189],
    'CH2=C': [0.64816, 0.11148],
    'CH=C': [1.22566, -0.03066],
    'C=C': [null, -0.12117],
    'CH2=C=CH': [-0.32258, null],
    'CH#C': [-0.74895, 0.43846],
    'C#C': [null, -0.35107],
    'ACH': [-0.19313, 0.13532],
    'AC': [0.16369, -0.17405],
    'ACCH3': [-0.47724, -0.28733],
    'ACCH2': [-0.33086, -0.88084],
    'ACCH<': [0.86718, -1.44666],
    'CH3COO': [1.71923, null],
    'CH2COO': [2.16274, null],
    'COO': [1.60913, 0.37204],
    'OH': [1.84013, null],
    'CH3O': [-0.4032, null],
    'CH2O': [null, null],
    'CH-O': [null, -0.40667],
    'CH2O-ring': [-0.33305, null],
    'CH2NH2': [null, null],
    'CH2NH': [0.83214, null],
    'CHNH': [1.25999, null],
    'CH3N': [null, -0.17004],
    'CH2N': [0.65229, -1.03686],
    'CH2S': [null, 0.14606],
    'CH2Cl': [null, 0.48952],
    'CHCl': [null, 0.12996],
    'CHCl2': [null, 0.52541],
    'ACCl': [-0.10778, 0.44238],
    'CCl2F': [null, null],
    'ACF': [null, -0.37183],
    'Cl-vinyl': [null, 0.66062],
    'CF3': [null, -0.08871],
    'CH2=C=C': [1.20154, null],
    'O': [-0.48942, null],
    'Cl': [null, 1.12515],
    'S': [0.11058, null],
    'C=O': [null, -0.05529]
  };
  const SP2L = {
    'iPr': [0.200694, 1e-08],
    'ring5': [-0.897912, 0.19438],
    'ring6': [-0.956852, 2e-08],
    'diene': [0.648793, null],
    'CH3-C=': [-0.008375, -0.06137],
    'CH2-C=': [0.011009, 0.06599],
    'CHorC-C=': [-0.3972, 0.342229],
    'cycstring': [null, -0.280859],
    'ACCOO': [0.491153, 1e-08],
    'naph': [0.013012, 0.086424],
    'acetal': [0.0, null]
  };

  // Flory–Huggins: chi = Vs·Ra² / (4RT); two liquids / a solute and a solvent mix when chi < chi_c = 0.5·(1 + 1/√n)²,
  // n = molar volume of the solute / molar volume of the solvent (small molecules tolerate a much larger Ra than polymers).
  function rankSolvents(dD, dP, dH, V) {
    return SOLV.map(s => {
      const ra = Math.sqrt(4 * (dD - s[1]) * (dD - s[1]) + (dP - s[2]) * (dP - s[2]) + (dH - s[3]) * (dH - s[3]));
      const chi = s[4] * ra * ra / (4 * RT), nn = Math.max(V / s[4], 0.2), chic = 0.5 * Math.pow(1 + 1 / Math.sqrt(nn), 2), red = chi / chic;
      return { name: s[0], ra, chi, chic, red, verdict: red < 1 ? 'good' : red < 1.5 ? 'partial' : 'poor' };
    }).sort((a, b) => a.red - b.red);
  }

  function hansenVK(mol, g, rg) {
    g = g || buildGraph(mol); rg = rg || findRings(g);
    const n = mol.atoms.length;
    if (!n) return { supported: false, reason: 'Empty drawing.' };
    if (mol.atoms.some(a => (a.charge || 0) !== 0 && !(a.el === 'N' || a.el === 'O'))) return { supported: false, reason: 'Charged species are not supported.' };
    const ar = aromaticity(mol, g, rg), used = new Set(), groups = [], bad = [], extra = [];
    const el = i => mol.atoms[i].el, lab = i => el(i) + (i + 1);
    const add1 = (name, atoms) => { groups.push({ name, atoms }); atoms.forEach(a => used.add(a)); };
    // 1. benzene rings with 0–4 substituents (not fused)
    ar.rings.forEach(r => {
      if (r.type !== 'aromatic' || r.size !== 6 || r.atoms.some(a => el(a) !== 'C')) return;
      if (rg.systems[rg.ringSys[ar.rings.indexOf(r)]].length > 1) return;
      const inRing = new Set(r.atoms);
      const k = r.atoms.reduce((s, a) => s + g.nbrs[a].filter(e => !inRing.has(e.j)).length, 0);
      const nm = ['C6H6', 'C6H5', 'C6H4', 'C6H3', 'C6H2'][k];
      if (nm) {
        add1(nm, r.atoms.slice());
        if (k >= 3) extra.push('Benzene ring with ' + k + ' substituents: its group values are extrapolated, not tabulated, so treat the result as less certain.');
      } else bad.push('benzene ring with ' + k + ' substituents (up to 4 are supported)');
    });
    for (let i = 0; i < n; i++) if (!used.has(i) && ar.rings.some(r => r.type === 'aromatic' && r.atoms.includes(i))) { bad.push('aromatic or heteroaromatic ring other than a benzene ring (fused rings and heteroaromatics are not tabulated)'); used.add(i); }
    // 2. carbon-centred functional groups
    for (let c = 0; c < n; c++) {
      if (used.has(c) || el(c) !== 'C') continue;
      const nb = g.nbrs[c], dO = nb.find(e => e.order === 2 && el(e.j) === 'O' && !used.has(e.j));
      const nitrile = nb.find(e => e.order === 3 && el(e.j) === 'N');
      if (nitrile) { add1('CN', [c, nitrile.j]); continue; }
      if (!dO) continue;
      const others = nb.filter(e => e !== dO), singles = others.filter(e => e.order === 1);
      const oh = singles.find(e => el(e.j) === 'O' && g.hc[e.j] === 1 && g.nbrs[e.j].length === 1);
      const oe = singles.find(e => el(e.j) === 'O' && g.hc[e.j] === 0 && g.nbrs[e.j].length === 2);
      const nn = singles.find(e => el(e.j) === 'N');
      if (nn) {
        // amide R–C(=O)–NH2 / R–C(=O)–NHR': built from the tabulated CO + NH2 / NH groups.
        // Not accepted: tertiary amides, ureas / carbamates (extra O or N on the carbonyl carbon), formamides, charged or unsaturated N.
        const N = nn.j, nN = g.nbrs[N], simple = !(mol.atoms[N].charge) && nN.every(e => e.order === 1) && g.hc[c] === 0 && singles.length === 2 && singles.every(e => e === nn || el(e.j) === 'C') && (g.hc[N] === 1 || g.hc[N] === 2);
        if (simple) {
          add1('CO', [c, dO.j]); add1(g.hc[N] === 2 ? 'NH2' : 'NH', [N]);
          if (!extra.some(t => t.indexOf('Amide') === 0)) extra.push('Amide group: estimated as a carbonyl (CO) plus an amine (NH / NH2) increment. Real amides hydrogen-bond more strongly, so δH is probably underestimated.');
        } else { bad.push('amide / urea / carbamate-type group at ' + lab(c) + ' (only primary and secondary amides are supported)'); [c, dO.j, N].forEach(a => used.add(a)); }
      }
      else if (oh && singles.every(e => e === oh || el(e.j) === 'C') && others.length <= 2) add1('COOH', [c, dO.j, oh.j]);
      else if (oe && others.length === 2 && singles.every(e => el(e.j) === 'C' || e === oe) ) add1('COO', [c, dO.j, oe.j]);
      else if (g.hc[c] === 1 && others.length === 1) { bad.push('aldehyde at ' + lab(c) + ' (not reliably tabulated)'); [c, dO.j].forEach(a => used.add(a)); }
      else if (others.length === 2 && singles.length === 2 && singles.every(e => el(e.j) === 'C')) add1('CO', [c, dO.j]);
      else { bad.push('carbonyl group at ' + lab(c) + ' (not tabulated)'); [c, dO.j].forEach(a => used.add(a)); }
    }
    // 3. N and O groups, halogens, other elements (nitro groups first, so their oxygens are not read as stray O−)
    for (let i = 0; i < n; i++) {
      if (used.has(i) || el(i) !== 'N' || (mol.atoms[i].charge || 0) !== 1) continue;
      const os = g.nbrs[i].filter(x => el(x.j) === 'O');
      if (os.length === 2 && g.nbrs[i].length === 3 && g.nbrs[i].filter(x => el(x.j) === 'C').length === 1) add1('NO2', [i].concat(os.map(x => x.j)));
    }
    for (let i = 0; i < n; i++) {
      if (used.has(i)) continue;
      const e = el(i), nb = g.nbrs[i], h = g.hc[i];
      if (e === 'N') {
        if (nb.some(x => x.order > 1) || (mol.atoms[i].charge || 0)) { bad.push('unsaturated / charged nitrogen ' + lab(i)); used.add(i); continue; }
        if (h === 2) add1('NH2', [i]); else if (h === 1) add1('NH', [i]); else { bad.push('tertiary amine ' + lab(i) + ' (not tabulated)'); used.add(i); }
      } else if (e === 'O') {
        if ((mol.atoms[i].charge || 0) || nb.some(x => x.order > 1)) { bad.push('unusual oxygen ' + lab(i)); used.add(i); }
        else if (h === 1) add1('OH', [i]); else if (nb.length === 2) add1('O', [i]); else { bad.push(nb.length ? 'oxygen ' + lab(i) : 'water / free oxygen (not tabulated)'); used.add(i); }
      } else if (e === 'F' || e === 'Cl' || e === 'Br') add1(e, [i]);
      else if (e !== 'C') { bad.push('element ' + e + ' (' + lab(i) + ')'); used.add(i); }
    }
    // 4. remaining carbons: CH3 / CH2 / CH / C and alkene carbons
    for (let c = 0; c < n; c++) {
      if (used.has(c) || el(c) !== 'C') continue;
      const nb = g.nbrs[c], h = g.hc[c];
      if (nb.some(e => e.order === 3) || nb.filter(e => e.order === 2).length > 1) { bad.push('alkyne / allene carbon ' + lab(c) + ' (not tabulated)'); used.add(c); continue; }
      if (nb.some(e => e.order === 2)) { add1(h === 2 ? '=CH2' : h === 1 ? '=CH' : '=C<', [c]); continue; }
      if (h === 4) { bad.push('methane (single carbon)'); used.add(c); continue; }
      add1(h === 3 ? 'CH3' : h === 2 ? 'CH2' : h === 1 ? 'CH' : 'C', [c]);
    }
    if (bad.length) return { supported: false, reason: 'Contains parts without group-contribution data: ' + Array.from(new Set(bad)).slice(0, 4).join('; ') + '.', unsupported: bad };
    let Fd = 0, Fp2 = 0, Eh = 0, V = 0; const partial = [], counts = {};
    groups.forEach(gr => { const p = HVK[gr.name]; Fd += p[0]; Fp2 += p[1] * p[1]; Eh += p[2]; V += p[3]; counts[gr.name] = (counts[gr.name] || 0) + 1; if (HVK_PARTIAL[gr.name] && partial.indexOf(gr.name) < 0) partial.push(gr.name); });
    if (V < 5) return { supported: false, reason: 'The molecule is too small or too branched for the molar-volume increments to give a sensible result.' };
    // several chlorines on one carbon partly cancel (empirical factors fitted to CH2Cl2, CHCl3, CCl4)
    const notes = extra.slice(), ps = [0, 1, 0.73, 0.4, 0];
    const clOn = new Map(); groups.forEach(x => { if (x.name === 'Cl') { const c = g.nbrs[x.atoms[0]][0]; if (c) clOn.set(c.j, (clOn.get(c.j) || 0) + 1); } });
    clOn.forEach(k => { if (k >= 2) { const fp = HVK.Cl[1]; Fp2 += Math.pow(fp * ps[Math.min(k, 4)], 2) - k * fp * fp; notes.push('Several Cl on one carbon: polar term reduced with an empirical factor.'); } });
    let dD = Fd / V, dP = Math.sqrt(Math.max(Fp2, 0)) / V, dH = Math.sqrt(Math.max(Eh, 0) / V), dT;
    if (dP > 3 && !groups.some(x => ['OH', 'NH2', 'NH', 'COOH'].includes(x.name))) {   // symmetric molecules: polar groups cancel (Van Krevelen symmetry rule)
      let mu = null; try { const d = dipole(mol); if (d.supported) mu = d.mu; } catch (e) { mu = null; }
      if (mu !== null && mu < 0.5) { dP *= 0.25; notes.push('Polar groups cancel by symmetry (dipole ≈ 0): δP scaled by 0.25.'); }
    }
    dT = Math.sqrt(dD * dD + dP * dP + dH * dH);
    // Flory–Huggins: chi = Vs·Ra² / (4RT); two liquids / a solute and a solvent mix when chi < chi_c = 0.5·(1 + 1/√n)²,
    // n = molar volume of the solute / molar volume of the solvent (small molecules tolerate a much larger Ra than polymers).
    const rows = rankSolvents(dD, dP, dH, V);
    return { supported: true, method: 'Hoftyzer–Van Krevelen', dD, dP, dH, dT, V, groups: counts, partial, notes, rows };
  }

  const SP_C = { d: 17.3231, p: 7.3548, h: 7.9793, pLow: 2.7467, hLow: 1.3720 };   // Eqs. 24–28
  const SP_HAL = ['F', 'Cl', 'Br', 'I'];

  // Splits the drawing into the paper's first-order groups (UNIFAC-type: a group may span several atoms) and
  // detects the second-order groups. Returns { first, second, bad } with group counts.
  function spFragment(mol, g, rg, ar) {
    const n = mol.atoms.length, el = i => mol.atoms[i].el, chg = i => mol.atoms[i].charge || 0;
    const nb = i => g.nbrs[i], hc = i => g.hc[i], lab = i => el(i) + (i + 1);
    const arom = new Set(); ar.rings.forEach(r => { if (r.type === 'aromatic') r.atoms.forEach(a => arom.add(a)); });
    const inRing = i => rg.rings.some(r => r.atoms.includes(i));
    const used = new Array(n).fill(false), gof = new Array(n).fill(null), first = {}, bad = [];
    const take = (key, atoms) => { first[key] = (first[key] || 0) + 1; atoms.forEach(a => { used[a] = true; gof[a] = key; }); };
    const skip = atoms => atoms.forEach(a => { used[a] = true; });
    const sp3 = i => !arom.has(i) && nb(i).every(e => e.order === 1);
    const freeC = i => el(i) === 'C' && !used[i];
    const pick = (list, h) => { const r = list.find(c => hc(c) === h); return r === undefined ? -1 : r; };
    // carbon that joins a hetero-atom group: CH3 first, then CH2 (then CH when allowed)
    const pairC = (cands, ch) => { const ok = cands.filter(c => freeC(c) && sp3(c)); let r = pick(ok, 3); if (r < 0) r = pick(ok, 2); if (r < 0 && ch) r = pick(ok, 1); return r; };

    // 1. nitro, nitrile, sulfone / sulfoxide, thiocarbonyl
    for (let N = 0; N < n; N++) {
      if (used[N] || el(N) !== 'N') continue;
      const os = nb(N).filter(e => el(e.j) === 'O' && nb(e.j).length === 1), cs = nb(N).filter(e => el(e.j) === 'C');
      if (os.length === 2 && cs.length === 1 && nb(N).length === 3) {
        const c = cs[0].j, ids = [N, os[0].j, os[1].j];
        if (arom.has(c)) take('ACNO2', ids.concat(c)); else if (hc(c) === 2) take('CH2NO2', ids.concat(c)); else if (hc(c) === 1) take('CHNO2', ids.concat(c));
        else { bad.push('nitro group on ' + lab(c) + ' (needs a CH2, CH or aromatic carbon)'); skip(ids); }
      }
    }
    for (let c = 0; c < n; c++) {
      if (used[c] || el(c) !== 'C') continue;
      const t = nb(c).find(e => e.order === 3 && el(e.j) === 'N'); if (!t) continue;
      const a = nb(c).find(e => e.j !== t.j);
      if (a && freeC(a.j) && sp3(a.j) && hc(a.j) === 2) take('CH2CN', [a.j, c, t.j]); else take('CN', [c, t.j]);
    }
    for (let s = 0; s < n; s++) {
      if (used[s] || el(s) !== 'S') continue;
      const os = nb(s).filter(e => el(e.j) === 'O' && e.order === 2 && !used[e.j]);
      if (os.length === 2) take('SO2', [s, os[0].j, os[1].j]);
      else if (os.length === 1) { bad.push('sulfoxide at ' + lab(s)); skip([s, os[0].j]); }
    }

    // 1b. structures the paper's tables cannot describe: exocyclic C=O / C=S on an aromatic ring (pyridone, xanthine, uracil type) and leftover charges
    for (let c = 0; c < n; c++) {
      if (!arom.has(c) || used[c]) continue;
      const x = nb(c).find(e => e.order === 2 && (el(e.j) === 'O' || el(e.j) === 'S') && !arom.has(e.j));
      if (x) { bad.push('C=' + el(x.j) + ' on an aromatic ring at ' + lab(c) + ' (pyridone / xanthine / uracil type)'); skip([c, x.j]); }
    }
    for (let i = 0; i < n; i++) if (!used[i] && chg(i) !== 0) { bad.push('charged atom ' + lab(i)); skip([i]); }

    // 2. carbonyl-type carbons: acid, ester, amide, aldehyde, ketone, thioketone
    for (let c = 0; c < n; c++) {
      if (used[c] || el(c) !== 'C' || arom.has(c)) continue;
      const dS = nb(c).find(e => e.order === 2 && el(e.j) === 'S' && !used[e.j]);
      if (dS) { take('C=S', [c, dS.j]); continue; }
      const dO = nb(c).find(e => e.order === 2 && el(e.j) === 'O' && !used[e.j]); if (!dO) continue;
      const others = nb(c).filter(e => e !== dO), carbons = others.filter(e => el(e.j) === 'C').map(e => e.j);
      const sO = others.filter(e => el(e.j) === 'O' && e.order === 1 && !used[e.j]);
      const acid = sO.find(e => nb(e.j).length === 1 && hc(e.j) === 1), ester = sO.find(e => nb(e.j).length === 2 && hc(e.j) === 0);
      const amN = others.find(e => el(e.j) === 'N' && e.order === 1 && !used[e.j]);
      const het = others.filter(e => (el(e.j) === 'N' || el(e.j) === 'O') && e.order === 1 && !used[e.j]).map(e => e.j);
      if (het.length >= 2) { bad.push('urea / carbamate / carbonate-type group at ' + lab(c) + ' (not in the paper’s group list)'); skip([c, dO.j].concat(het)); continue; }
      if (acid) { take('COOH', [c, dO.j, acid.j]); continue; }
      if (ester) {
        const a = carbons.length ? carbons[0] : -1, ok = a >= 0 && freeC(a) && sp3(a);
        if (hc(c) === 1) take('HCOO', [c, dO.j, ester.j]);
        else if (ok && hc(a) === 3) take('CH3COO', [a, c, dO.j, ester.j]);
        else if (ok && hc(a) === 2) take('CH2COO', [a, c, dO.j, ester.j]);
        else take('COO', [c, dO.j, ester.j]);
        continue;
      }
      if (amN) {
        const N = amN.j, ns = nb(N).filter(e => e.j !== c);
        if (hc(N) === 2 && !ns.length) { take('CONH2', [c, dO.j, N]); continue; }
        const me = ns.filter(e => el(e.j) === 'C' && freeC(e.j) && sp3(e.j) && hc(e.j) === 3 && nb(e.j).length === 1);
        if (hc(N) === 0 && ns.length === 2 && me.length === 2) { take('CON(CH3)2', [c, dO.j, N, me[0].j, me[1].j]); continue; }
        take('C=O', [c, dO.j]); continue;                       // other amides: generic carbonyl, nitrogen handled below
      }
      if (hc(c) === 1) { take('CHO', [c, dO.j]); continue; }
      const a3 = carbons.find(x => freeC(x) && sp3(x) && hc(x) === 3), a2 = carbons.find(x => freeC(x) && sp3(x) && hc(x) === 2);
      if (a3 !== undefined) take('CH3CO', [a3, c, dO.j]); else if (a2 !== undefined) take('CH2CO', [a2, c, dO.j]); else take('C=O', [c, dO.j]);
    }

    // 3. alcohols / phenols, ethers
    for (let o = 0; o < n; o++) {
      if (used[o] || el(o) !== 'O' || hc(o) !== 1 || nb(o).length !== 1) continue;
      const c = nb(o)[0].j;
      if (el(c) === 'C' && arom.has(c) && !used[c]) take('ACOH', [c, o]); else take('OH', [o]);
    }
    for (let o = 0; o < n; o++) {
      if (used[o] || el(o) !== 'O' || hc(o) !== 0 || nb(o).length !== 2 || nb(o).some(e => e.order !== 1)) continue;
      const c = pairC(nb(o).map(e => e.j), true);
      if (c < 0) take('O', [o]);
      else { const h = hc(c); take(h === 3 ? 'CH3O' : h === 2 ? (inRing(o) && inRing(c) ? 'CH2O-ring' : 'CH2O') : 'CH-O', [c, o]); }
    }

    // 4. nitrogen
    for (let N = 0; N < n; N++) {
      if (used[N] || el(N) !== 'N') continue;
      const ns = nb(N), h = hc(N);
      if (chg(N) !== 0) { bad.push('charged nitrogen ' + lab(N)); skip([N]); continue; }
      if (arom.has(N)) {
        if (ns.length === 2 && h === 0) {                        // pyridine-type: N is grouped with one neighbouring ring carbon
          const cs = ns.map(e => e.j).filter(c => !used[c]);
          let c = cs.find(x => hc(x) > 0); if (c === undefined) c = cs[0];
          if (c === undefined) { bad.push('aromatic nitrogen ' + lab(N)); skip([N]); } else take(hc(c) > 0 ? 'CH=N' : 'C=N', [c, N]);
        } else take(h === 1 ? 'NH' : 'N', [N]);
        continue;
      }
      const dbl = ns.find(e => e.order === 2 && el(e.j) === 'C' && !used[e.j]);
      if (dbl) { take(hc(dbl.j) > 0 ? 'CH=N' : 'C=N', [dbl.j, N]); continue; }
      if (ns.some(e => e.order > 1)) { bad.push('N=X / N≡X group at ' + lab(N)); skip([N]); continue; }
      const cs = ns.map(e => e.j);
      if (h === 2) {
        const c = cs[0];
        if (arom.has(c) && !used[c]) take('ACNH2', [c, N]);
        else if (freeC(c) && sp3(c) && hc(c) === 2) take('CH2NH2', [c, N]);
        else if (freeC(c) && sp3(c) && hc(c) === 1) take('CHNH2', [c, N]);
        else { bad.push('primary amine on ' + lab(c) + ' (needs a CH2, CH or aromatic carbon)'); skip([N]); }
      } else if (h === 1) { const c = pairC(cs, false); if (c < 0) take('NH', [N]); else take(hc(c) === 3 ? 'CH3NH' : 'CH2NH', [c, N]); }
      else if (h === 0 && cs.length === 3) { const c = pairC(cs, false); if (c < 0) take('N', [N]); else take(hc(c) === 3 ? 'CH3N' : 'CH2N', [c, N]); }
      else { bad.push('nitrogen ' + lab(N)); skip([N]); }
    }

    // 5. sulfur
    for (let s = 0; s < n; s++) {
      if (used[s] || el(s) !== 'S') continue;
      const ns = nb(s);
      if (ns.some(e => e.order > 1) || chg(s)) { bad.push('sulfur group at ' + lab(s)); skip([s]); continue; }
      if (hc(s) === 1 && ns.length === 1) { const c = ns[0].j; if (freeC(c) && sp3(c) && hc(c) === 2) take('CH2SH', [c, s]); else take('SH', [s]); }
      else if (ns.length === 2 && hc(s) === 0) {
        if (arom.has(s)) take('S', [s]);
        else { const c = pairC(ns.map(e => e.j), false); if (c < 0) take('S', [s]); else take(hc(c) === 3 ? 'CH3S' : 'CH2S', [c, s]); }
      } else { bad.push('sulfur ' + lab(s)); skip([s]); }
    }

    // 6. halogens
    for (let c = 0; c < n; c++) {
      if (el(c) !== 'C') continue;
      const hs = nb(c).filter(e => SP_HAL.includes(el(e.j)) && !used[e.j]); if (!hs.length) continue;
      if (arom.has(c)) {
        hs.forEach(e => { const X = el(e.j); if (X === 'Br' || X === 'I') take(X, [e.j]); else if (!used[c]) take(X === 'Cl' ? 'ACCl' : 'ACF', [c, e.j]); else take(X === 'Cl' ? 'Cl' : 'F', [e.j]); });
        continue;
      }
      if (used[c]) { hs.forEach(e => take(el(e.j) === 'Cl' || el(e.j) === 'F' ? el(e.j) : el(e.j), [e.j])); continue; }
      const X = k => hs.filter(e => el(e.j) === k).map(e => e.j), Cl = X('Cl'), F = X('F'), h = hc(c), vinyl = nb(c).some(e => e.order === 2);
      X('Br').concat(X('I')).forEach(x => take(el(x), [x]));
      if (Cl.length && !F.length) {
        if (vinyl && Cl.length === 1) take('Cl-vinyl', [Cl[0]]);
        else if (Cl.length === 1 && h <= 2) take(h === 2 ? 'CH2Cl' : h === 1 ? 'CHCl' : 'CCl', [c, Cl[0]]);
        else if (Cl.length === 2 && h <= 1) take(h === 1 ? 'CHCl2' : 'CCl2', [c].concat(Cl));
        else if (Cl.length === 3 && h === 0) take('CCl3', [c].concat(Cl));
        else Cl.forEach(x => take('Cl', [x]));
      } else if (F.length && !Cl.length) {
        if (F.length === 3 && h === 0) take('CF3', [c].concat(F));
        else if (F.length === 2 && h === 0) take('CF2', [c].concat(F));
        else if (F.length === 1 && h === 0) take('CF', [c, F[0]]);
        else F.forEach(x => take('F', [x]));
      } else { Cl.forEach(x => take('Cl', [x])); F.forEach(x => take('F', [x])); }
    }

    // 7. hydrocarbon skeleton: alkynes, allenes, aromatic carbons (with a CH3 / CH2 substituent), alkenes, saturated carbons
    for (let c = 0; c < n; c++) {
      if (!freeC(c) || arom.has(c)) continue;
      const t = nb(c).find(e => e.order === 3 && el(e.j) === 'C' && freeC(e.j));
      if (t && c < t.j) take(hc(c) === 1 || hc(t.j) === 1 ? 'CH#C' : 'C#C', [c, t.j]);
    }
    for (let m = 0; m < n; m++) {
      if (!freeC(m) || arom.has(m)) continue;
      const d = nb(m).filter(e => e.order === 2 && freeC(e.j));
      if (d.length === 2 && nb(m).length === 2) {
        const ends = d.map(e => e.j).sort((a, b) => hc(b) - hc(a));
        if (hc(ends[0]) === 2 && hc(ends[1]) === 1) take('CH2=C=CH', [ends[0], m, ends[1]]);
        else if (hc(ends[0]) === 2 && hc(ends[1]) === 0) take('CH2=C=C', [ends[0], m, ends[1]]);
        else { bad.push('allene at ' + lab(m)); skip([m].concat(ends)); }
      }
    }
    for (let c = 0; c < n; c++) {
      if (!arom.has(c) || used[c]) continue;
      const s = nb(c).map(e => e.j).find(j => !arom.has(j) && freeC(j) && sp3(j) && (hc(j) === 3 || hc(j) === 2));
      if (s !== undefined) take(hc(s) === 3 ? 'ACCH3' : 'ACCH2', [c, s]); else take(hc(c) > 0 ? 'ACH' : 'AC', [c]);
    }
    for (let c = 0; c < n; c++) {
      if (!freeC(c) || arom.has(c)) continue;
      const d = nb(c).find(e => e.order === 2 && freeC(e.j) && !arom.has(e.j));
      if (!d || c > d.j) continue;
      const hs = [hc(c), hc(d.j)].sort((a, b) => b - a), key = hs[0] === 2 && hs[1] === 1 ? 'CH2=CH' : hs[0] === 1 && hs[1] === 1 ? 'CH=CH' : hs[0] === 2 && hs[1] === 0 ? 'CH2=C' : hs[0] === 1 && hs[1] === 0 ? 'CH=C' : hs[0] === 0 && hs[1] === 0 ? 'C=C' : null;
      if (key) take(key, [c, d.j]); else { bad.push('C=C double bond at ' + lab(c) + ' (ethylene-type)'); skip([c, d.j]); }
    }
    for (let c = 0; c < n; c++) {
      if (!freeC(c)) continue;
      if (!sp3(c)) { bad.push('carbon ' + lab(c) + ' (unrecognised bonding)'); skip([c]); continue; }
      const h = hc(c);
      if (h === 4) { bad.push('methane (single carbon)'); skip([c]); continue; }
      take(h === 3 ? 'CH3' : h === 2 ? 'CH2' : h === 1 ? 'CH' : 'C', [c]);
    }
    for (let i = 0; i < n; i++) {                                   // leftovers: lone halogens / oxygens / other elements
      if (used[i]) continue;
      const e = el(i);
      if (SP_HAL.includes(e)) take(e, [i]);
      else if (e === 'O' && nb(i).length) take('O', [i]);
      else { bad.push(e === 'O' ? 'water / free oxygen' : 'element ' + e + ' (' + lab(i) + ')'); skip([i]); }
    }

    // second-order groups (structural patterns)
    const second = {}, add = k => { second[k] = (second[k] || 0) + 1; };
    const cN = i => nb(i).filter(e => el(e.j) === 'C').map(e => e.j);
    const alkC = c => el(c) === 'C' && !arom.has(c) && nb(c).some(e => e.order === 2 && el(e.j) === 'C');
    const plain = ['CH3', 'CH2', 'CH', 'C'];
    const dO = c => nb(c).some(e => e.order === 2 && el(e.j) === 'O');
    const isKetone = c => el(c) === 'C' && !arom.has(c) && dO(c) && hc(c) === 0 && nb(c).filter(e => e.order === 1).every(e => el(e.j) === 'C');
    const isEster = c => el(c) === 'C' && dO(c) && nb(c).some(e => e.order === 1 && el(e.j) === 'O' && nb(e.j).length === 2);
    for (let c = 0; c < n; c++) {
      const e = el(c), k = gof[c];
      if (e === 'C') {
        if (k === 'CH' && sp3(c) && cN(c).filter(x => sp3(x) && hc(x) === 3).length >= 2) add('iPr');
        if (k === 'C' && sp3(c) && cN(c).filter(x => sp3(x) && hc(x) === 3).length >= 4 - 1) add('tBu');
        if (plain.includes(k) && sp3(c) && cN(c).some(alkC)) add(k === 'CH3' ? 'CH3-C=' : k === 'CH2' ? 'CH2-C=' : 'CHorC-C=');
        if (alkC(c)) cN(c).forEach(x => { if (x > c && alkC(x) && nb(c).some(f => f.j === x && f.order === 1)) add('diene'); });
        if (k === 'COOH') cN(c).forEach(x => { if (arom.has(x)) add('ACCOOH'); else if (sp3(x) && hc(x) <= 1) add('CHorC-COOH'); });
        if (isEster(c) && !arom.has(c) && !nb(c).some(f => el(f.j) === 'N')) {
          const a = cN(c).find(x => arom.has(x)); if (a !== undefined) add('ACCOO');
          if (k === 'CH3COO') { const o = nb(c).find(f => f.order === 1 && el(f.j) === 'O' && nb(f.j).length === 2), x = o && nb(o.j).find(f => f.j !== c); if (x && el(x.j) === 'C' && sp3(x.j) && hc(x.j) <= 1) add('CH3COO-CHorC'); }
        }
        if (k === 'CHO' && cN(c).some(x => arom.has(x))) add('ACHO');
        if (isKetone(c) && inRing(c)) add('Ccyc=O');
        if (isKetone(c)) { const cs = cN(c); if (cs.length === 2 && cs.some(x => sp3(x) && hc(x) === 3) && cs.some(x => sp3(x) && hc(x) === 2)) add('CH3COCH2'); }
        if (sp3(c) && hc(c) === 2 && cN(c).length === 2) {
          const cs = cN(c);
          if (cs.some(isKetone) && cs.some(isEster)) add('COCH2COO');
        }
        if (sp3(c) && cN(c).filter(isKetone).length === 2) add('dicarbonyl');
        if (sp3(c) && nb(c).filter(f => el(f.j) === 'O' && nb(f.j).length === 2 && f.order === 1 && !arom.has(f.j)).length === 2) add('acetal');
        if (sp3(c)) cN(c).forEach(x => {
          if (x < c || !sp3(x)) return;
          const oh = i => nb(i).some(f => el(f.j) === 'O' && hc(f.j) === 1 && nb(f.j).length === 1), am = i => nb(i).some(f => el(f.j) === 'N' && f.order === 1 && !arom.has(f.j) && !nb(f.j).some(q => dO(q.j)));
          if (oh(c) && oh(x)) add('diol'); else if ((oh(c) && am(x)) || (oh(x) && am(c))) add('aminoalcohol');
        });
      } else if (e === 'O') {
        if (k === 'OH' || k === 'ACOH') { if (k === 'OH') { const x = nb(c)[0].j; if (el(x) === 'C' && sp3(x)) add(inRing(x) ? 'Ccyc-OH' : hc(x) === 1 ? 'CHOH' : hc(x) === 0 ? 'C(OH)' : 'x'); } }
        else if (nb(c).length === 2 && hc(c) === 0 && !arom.has(c)) {
          const cs = nb(c).map(f => f.j);
          if (cs.every(x => el(x) === 'C' && arom.has(x))) add('AC-O-AC');
          else if (cs.some(x => el(x) === 'C' && arom.has(x)) && cs.some(x => el(x) === 'C' && sp3(x))) add('AC-O-C');
          else if (cs.some(alkC) && cs.some(x => el(x) === 'C' && sp3(x))) add('C-O-C=C');
          else if (cs.every(x => isEster(x) || (el(x) === 'C' && dO(x)))) add('anhydride');
        }
      } else if (e === 'N') {
        if (!arom.has(c) && inRing(c) && nb(c).every(f => f.order === 1) && !nb(c).some(f => dO(f.j))) add('Ncyc');
        if (!arom.has(c) && inRing(c) && nb(c).some(f => dO(f.j) && rg.rings.some(r => r.atoms.includes(c) && r.atoms.includes(f.j)))) add('lactam');
      } else if (e === 'S') { if (!arom.has(c) && inRing(c) && nb(c).every(f => f.order === 1) && nb(c).length === 2) add('Scyc'); }
      else if (e === 'Br') { const x = nb(c)[0] && nb(c)[0].j; if (x !== undefined && el(x) === 'C') { if (arom.has(x)) add('ACBr'); else if (alkC(x)) add('C=C-Br'); } }
    }
    rg.rings.forEach((r, ri) => {
      if (r.atoms.every(a => el(a) === 'C' && !arom.has(a))) { if (r.size === 3) add('ring3'); else if (r.size === 5) add('ring5'); else if (r.size === 6) add('ring6'); }
      if (r.size === 6 && r.atoms.every(a => arom.has(a)) && r.atoms.filter(a => el(a) === 'N').length === 1) add('pyr');
      for (let rj = ri + 1; rj < rg.rings.length; rj++) {
        const q = rg.rings[rj];
        if (r.size === 6 && q.size === 6 && [r, q].every(x => x.atoms.every(a => arom.has(a) && el(a) === 'C')) && r.bonds.some(b => q.bonds.includes(b))) add('naph');
      }
    });
    delete second.x;
    return { first, second, bad: Array.from(new Set(bad)) };
  }

  // Le Bas molar volume (cm³/mol): used only to size the Flory–Huggins mixing limit, accuracy about ±10 %
  function leBasVolume(mol, g, rg) {
    let V = 0;
    mol.atoms.forEach((a, i) => {
      const nb = g.nbrs[i]; V += 3.7 * g.hc[i];
      if (a.el === 'C') V += 14.8;
      else if (a.el === 'N') V += nb.some(e => e.order > 1) ? 15.6 : nb.length + g.hc[i] === 3 ? (nb.length === 1 ? 10.5 : nb.length === 2 ? 12.0 : 15.6) : 15.6;
      else if (a.el === 'O') V += nb.length === 2 && nb.every(e => e.order === 1) ? 9.1 : 7.4;
      else V += { S: 25.6, F: 8.7, Cl: 24.6, Br: 27.0, I: 37.0, P: 27.0 }[a.el] || 20;
    });
    rg.rings.forEach(r => { V -= r.size === 3 ? 6 : r.size === 4 ? 8.5 : r.size === 5 ? 11.5 : 15; });
    return Math.max(V, 10);
  }

  function hansenSP(mol, g, rg) {
    const n = mol.atoms.length, M = 'Stefanis–Panayiotou';
    if (!n) return { supported: false, reason: 'Empty drawing.' };
    const ar = aromaticity(mol, g, rg), fr = spFragment(mol, g, rg, ar);
    if (fr.bad.length) return { supported: false, reason: M + ' method: no data for ' + fr.bad.slice(0, 4).join('; ') + '.' };
    const W = Object.keys(fr.second).length ? 1 : 0, notes = [];
    const sum = (idx, c0, T1, T2, strict) => {
      let s = c0;
      for (const k of Object.keys(fr.first)) { const row = T1[k], v = row ? row[idx] : null; if (v === null || v === undefined) { if (strict) return { s: null, miss: k }; return null; } s += fr.first[k] * v; }
      if (W) for (const k of Object.keys(fr.second)) { const row = T2[k], v = row ? row[idx] : null; if (v !== null && v !== undefined) s += fr.second[k] * v; }
      return strict ? { s, miss: null } : s;
    };
    const dd = sum(0, SP_C.d, SP1, SP2, true);
    const WHY = { NH: ', which is what an amide N–H, a pyrrole / indole / imidazole N–H or a diarylamine N–H is classed as' };
    if (dd.s === null) return { supported: false, reason: M + ' method: the paper gives no δD contribution for the group “' + dd.miss + '”' + (WHY[dd.miss] || '') + '.' };
    const out = {};
    [['p', 1, SP_C.p, SP_C.pLow, 0], ['h', 2, SP_C.h, SP_C.hLow, 1]].forEach(([nm, idx, c0, cLow, li]) => {
      const main = sum(idx, c0, SP1, SP2, true);
      if (main.s === null) { out[nm] = { miss: main.miss }; return; }
      let v = main.s, usedLow = false;
      // low-value equations (Eqs. 27 / 28, Tables 5 / 6); only when Table 6 has an entry for every second-order group present
      const lw = Object.keys(fr.second).every(k => SP2L[k]) ? sum(li, cLow, SP1L, SP2L, false) : null;
      if (lw !== null && lw < 3) { v = lw; usedLow = true; }
      out[nm] = { v: Math.max(0, v), usedLow, main: main.s, low: lw };
    });
    if (out.p.miss || out.h.miss) return { supported: false, reason: M + ' method: the paper gives no ' + (out.p.miss ? 'δP' : 'δH') + ' contribution for the group “' + (out.p.miss || out.h.miss) + '”.' };
    const dD = dd.s, dP = out.p.v, dH = out.h.v, dT = Math.sqrt(dD * dD + dP * dP + dH * dH);
    const nC = mol.atoms.filter(a => a.el === 'C').length;
    if (nC < 3) notes.push('The Stefanis–Panayiotou method was fitted to molecules with three or more carbon atoms; expect larger errors here.');
    if (out.p.usedLow || out.h.usedLow) notes.push('Low-value equation used for ' + [out.p.usedLow ? 'δP' : '', out.h.usedLow ? 'δH' : ''].filter(Boolean).join(' and ') + ' (estimate below 3 MPa½).');
    const counts = {}; Object.keys(fr.first).forEach(k => { counts[k] = fr.first[k]; }); Object.keys(fr.second).forEach(k => { counts['2nd: ' + k] = fr.second[k]; });
    const V = leBasVolume(mol, g, rg);
    return { supported: true, method: M, dD, dP, dH, dT, V, groups: counts, partial: [], notes, rows: rankSolvents(dD, dP, dH, V), detail: { p: out.p, h: out.h } };
  }

  // Method switch: 'auto' (Stefanis–Panayiotou, the broader published method, falling back to Hoftyzer–Van Krevelen
  // where the paper has no value), 'sp' or 'vk'. In 'auto' mode a result also carries `alt`: the other method's
  // numbers when it can handle the structure too, so the two estimates can be compared.
  function hansen(mol, g, rg, opt) {
    g = g || buildGraph(mol); rg = rg || findRings(g);
    const m = (opt && opt.method) || 'auto';
    if (m === 'sp') return hansenSP(mol, g, rg);
    if (m === 'vk') return hansenVK(mol, g, rg);
    const sp = hansenSP(mol, g, rg), vk = hansenVK(mol, g, rg);
    const brief = x => ({ method: x.method, dD: x.dD, dP: x.dP, dH: x.dH, dT: x.dT });
    if (sp.supported) { sp.alt = vk.supported ? brief(vk) : null; return sp; }
    if (vk.supported) { vk.notes.unshift('The Stefanis–Panayiotou method could not be used: ' + sp.reason.replace(/^Stefanis–Panayiotou method: /, '').replace(/\.$/, '') + '. Hoftyzer–Van Krevelen was used instead.'); vk.alt = null; return vk; }
    return { supported: false, reason: sp.reason + ' Hoftyzer–Van Krevelen method: ' + vk.reason };
  }

  /* --------------------------------------------------------- file export */
  function toMolBlock(mol, Lpx) {
    const L = Lpx || 40, f = x => x.toFixed(4).padStart(10);
    const g = buildGraph(mol);
    const out = ['drawn molecule', '  chem.js', ''];
    out.push(String(mol.atoms.length).padStart(3) + String(g.bonds.length).padStart(3) + '  0  0  0  0  0  0  0  0999 V2000');
    mol.atoms.forEach(a => out.push(f(a.x / L * 1.5) + f(-a.y / L * 1.5) + f(0) + ' ' + a.el.padEnd(3) + ' 0  0  0  0  0  0  0  0  0  0  0  0'));
    g.bonds.forEach(b => out.push(String(b.i + 1).padStart(3) + String(b.j + 1).padStart(3) + String(b.order).padStart(3) + String(b.stereo === 1 ? 1 : b.stereo === 2 ? 6 : 0).padStart(3)));
    mol.atoms.forEach((a, i) => { if (a.charge) out.push('M  CHG  1 ' + String(i + 1).padStart(3) + ' ' + String(a.charge).padStart(3)); });
    out.push('M  END');
    return out.join('\n');
  }

  function toXYZ(s) {
    const out = [String(s.atoms.length), 'generated by chem.js'];
    s.atoms.forEach(a => out.push(a.el.padEnd(2) + ' ' + a.x.toFixed(5).padStart(11) + ' ' + a.y.toFixed(5).padStart(11) + ' ' + a.z.toFixed(5).padStart(11)));
    return out.join('\n');
  }

  root.Chem = { EL, valences, buildGraph, formula, findRings, aromaticity, clean2D, embed3D, dipole, ezStereo, rsStereo, hansen, toMolBlock, toXYZ, _params: { PI_DONOR } };
  if (typeof module !== 'undefined' && module.exports) module.exports = root.Chem;
})(typeof window !== 'undefined' ? window : globalThis);
