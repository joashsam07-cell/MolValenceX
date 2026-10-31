/* smiles.js — SMILES import / export for MolValenceX. Needs chem.js loaded first (window.Chem). */
(function (root) {
  'use strict';
  const C = root.Chem;
  if (!C) return;
  const ORG = ['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'];
  const AROM = ['b', 'c', 'n', 'o', 'p', 's', 'se', 'as'];
  const DEFV = { B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6], F: [1], Cl: [1], Br: [1], I: [1] };
  const cap = s => s[0].toUpperCase() + s.slice(1);

  /* ------------------------------------------------------------ parse */
  function parse(str) {
    const s = String(str).trim().split(/\s+/)[0];
    if (!s) throw new Error('empty string');
    const atoms = [], bonds = [], stack = [], open = {}, nbo = [];   // nbo[i]: neighbours of atom i in SMILES order (-1 = the implicit H)
    let prev = -1, sym = null, i = 0;
    const link = (a, b, sy) => bonds.push({ a, b, sy });
    const addAtom = a => {
      const k = atoms.length; atoms.push(a); nbo.push([]);
      if (prev >= 0) { link(prev, k, sym); nbo[prev].push(k); nbo[k].push(prev); }
      if (a.tag && a.h === 1) nbo[k].push(-1);                        // [C@H]: the H comes right after the preceding atom
      prev = k; sym = null;
    };
    while (i < s.length) {
      const ch = s[i];
      if (ch === '(') { if (prev < 0) throw new Error('branch without atom'); stack.push(prev); i++; continue; }
      if (ch === ')') { if (!stack.length) throw new Error('unbalanced )'); prev = stack.pop(); sym = null; i++; continue; }
      if ('-=#:/\\'.indexOf(ch) >= 0) { sym = ch; i++; continue; }
      if (ch === '.') { prev = -1; sym = null; i++; continue; }
      if (ch === '%' || /\d/.test(ch)) {
        let d; if (ch === '%') { d = s.substr(i + 1, 2); i += 3; } else { d = ch; i++; }
        if (!/^\d+$/.test(d)) throw new Error('bad ring label');
        if (prev < 0) throw new Error('ring digit without atom');
        if (open[d]) { const o = open[d]; delete open[d]; if (o.a === prev) throw new Error('ring closes on itself'); link(o.a, prev, sym || o.sy); nbo[o.a][o.pos] = prev; nbo[prev].push(o.a); }
        else { open[d] = { a: prev, sy: sym, pos: nbo[prev].length }; nbo[prev].push(-2); }   // -2 = ring-closure placeholder, filled in when the ring closes
        sym = null; continue;
      }
      if (ch === '[') {
        const j = s.indexOf(']', i); if (j < 0) throw new Error('missing ]');
        const m = /^\[(\d*)([A-Za-z][a-z]?)(@{0,2})(?:H(\d*))?([+-]+\d*|[+-]*)(?::\d+)?\]$/.exec(s.slice(i, j + 1));
        if (!m) throw new Error('cannot read ' + s.slice(i, j + 1));
        const e = m[2], arom = /^[a-z]/.test(e), El = cap(e);
        if (arom && AROM.indexOf(e) < 0) throw new Error('bad aromatic atom ' + e);
        if (!C.EL[El]) throw new Error('element ' + El + ' is not supported');
        let q = 0; const cs = m[5];
        if (/^[+-]\d+$/.test(cs)) q = (cs[0] === '+' ? 1 : -1) * +cs.slice(1); else q = (cs.split('+').length - 1) - (cs.split('-').length - 1);
        addAtom({ el: El, arom, charge: q, h: m[4] === undefined ? 0 : (m[4] === '' ? 1 : +m[4]), br: true, tag: m[3] || '' });
        i = j + 1; continue;
      }
      let el = null, arom = false;
      if (ch === 'C' && s[i + 1] === 'l') { el = 'Cl'; i += 2; }
      else if (ch === 'B' && s[i + 1] === 'r') { el = 'Br'; i += 2; }
      else if (ORG.indexOf(ch) >= 0) { el = ch; i++; }
      else if ('bcnops'.indexOf(ch) >= 0) { el = ch.toUpperCase(); arom = true; i++; }
      else throw new Error('unexpected "' + ch + '"');
      addAtom({ el, arom, charge: 0, h: 0, br: false });
    }
    if (stack.length) throw new Error('unbalanced (');
    if (Object.keys(open).length) throw new Error('unclosed ring digit');
    if (!atoms.length) throw new Error('no atoms');
    bonds.forEach(b => {
      b.order = b.sy === '=' ? 2 : b.sy === '#' ? 3 : b.sy === ':' ? 'ar' : b.sy ? 1 : (atoms[b.a].arom && atoms[b.b].arom ? 'ar' : 1);
    });
    // kekulise aromatic bonds
    const nbr = atoms.map(() => []); bonds.forEach((b, k) => { nbr[b.a].push(k); nbr[b.b].push(k); });
    const total = i2 => nbr[i2].reduce((t, k) => t + (bonds[k].order === 'ar' ? 1 : bonds[k].order), 0);
    const canDouble = i2 => {
      const a = atoms[i2], t = total(i2);
      if (a.br) return C.valences(a.el, a.charge).some(v => v === t + a.h + 1);
      if (a.el === 'C' || a.el === 'B') return t <= 3;
      if (a.el === 'N' || a.el === 'P') return t === 2;
      return false;
    };
    const need = new Set();
    atoms.forEach((a, i2) => { if (a.arom && nbr[i2].some(k => bonds[k].order === 'ar') && canDouble(i2)) need.add(i2); });
    const list = Array.from(need), matched = new Set(), dbl = new Set();
    (function () { })();
    const solve = () => {
      const u = list.find(x => !matched.has(x)); if (u === undefined) return true;
      for (const k of nbr[u]) {
        const b = bonds[k]; if (b.order !== 'ar') continue;
        const v = b.a === u ? b.b : b.a; if (matched.has(v) || !need.has(v)) continue;
        matched.add(u); matched.add(v); dbl.add(k);
        if (solve()) return true;
        matched.delete(u); matched.delete(v); dbl.delete(k);
      }
      return false;
    };
    if (!solve()) throw new Error('could not kekulise the aromatic rings');
    bonds.forEach((b, k) => { if (b.order === 'ar') b.order = dbl.has(k) ? 2 : 1; });
    return { atoms, bonds, nbo };
  }

  /* ------------------------------------------------ tetrahedral stereo (used by both import and export) */
  // SMILES '@' = looking from the first neighbour, the other three run anticlockwise; '@@' = clockwise.
  // `xy` = atom positions on the canvas (y down), `order` = the neighbours of centre c in SMILES order (-1 = the implicit H),
  // `zmap` = z for neighbours joined by a wedge (+1) / dash (-1) that starts at c. Same geometry model as Chem.rsStereo.
  const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const sub3 = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const det4 = V => dot3(sub3(V[1], V[0]), cross3(sub3(V[2], V[0]), sub3(V[3], V[0])));
  const REF = det4([[0, 0, 1], [1, 0, -0.33], [-0.5, 0.866, -0.33], [-0.5, -0.866, -0.33]]);   // viewed from the first atom, the rest anticlockwise = '@'
  function chirTag(xy, c, order, zmap) {
    const V = order.map(j => {
      if (j < 0) return null;
      const dx = xy[j][0] - xy[c][0], dy = -(xy[j][1] - xy[c][1]), l = Math.hypot(dx, dy) || 1;
      return [dx / l, dy / l, (zmap.get(j) || 0) * 0.9];
    });
    const hi = V.indexOf(null);
    if (hi >= 0) {                                                  // implicit H points away from the three drawn substituents
      let s = [0, 0, 0];
      V.forEach(v => { if (v) { const l = Math.hypot(v[0], v[1], v[2]) || 1; s = [s[0] + v[0] / l, s[1] + v[1] / l, s[2] + v[2] / l]; } });
      V[hi] = [-s[0], -s[1], -s[2]];
    }
    const d = det4(V);
    return { tag: d * REF > 0 ? '@' : '@@', mag: Math.abs(d) };
  }

  // Turns the @ / @@ marks of a parsed SMILES into wedge / dash bonds on the finished 2D drawing.
  function assignWedges(p, lay) {
    const xy = lay.atoms.map(a => [a.x, a.y]), idx = new Map(), st = lay.bonds.map(() => 0), from = lay.bonds.map(() => -1);
    lay.bonds.forEach((b, k) => { idx.set(b.a + ',' + b.b, k); idx.set(b.b + ',' + b.a, k); });
    const m0 = { atoms: lay.atoms.map((a, i) => ({ id: i, el: a.el, charge: a.charge, x: a.x, y: a.y })), bonds: lay.bonds.map((b, k) => ({ id: k, a: b.a, b: b.b, order: b.order, stereo: 0 })) };
    const g = C.buildGraph(m0), rg = C.findRings(g), inRing = new Set();
    rg.rings.forEach(r => r.bonds.forEach(bi => inRing.add(g.bonds[bi].mi)));
    const heavyDeg = j => p.nbo[j].filter(x => x >= 0).length;
    let skipped = 0;
    p.atoms.forEach((A, c) => {
      if (!A.tag) return;
      const order = p.nbo[c];
      if (order.length !== 4 || order.indexOf(-2) >= 0) { skipped++; return; }      // needs three drawn neighbours + H, or four
      const cands = order.filter(j => j >= 0).map(j => ({ j, k: idx.get(c + ',' + j) }))
        .filter(x => x.k !== undefined && !st[x.k] && lay.bonds[x.k].order === 1)
        .sort((u, v) => (inRing.has(u.k) - inRing.has(v.k)) || (heavyDeg(u.j) - heavyDeg(v.j)));   // prefer a non-ring, terminal substituent
      for (const x of cands) {
        let best = null;
        [1, -1].forEach(z => { const r = chirTag(xy, c, order, new Map([[x.j, z]])); if (r.tag === A.tag && r.mag > 0.15 && (!best || r.mag > best.mag)) best = { z, mag: r.mag }; });
        if (best) { st[x.k] = best.z > 0 ? 1 : 2; from[x.k] = c; return; }
      }
      skipped++;
    });
    lay.bonds.forEach((b, k) => { if (st[k]) { b.stereo = st[k]; if (b.a !== from[k]) { const t = b.a; b.a = b.b; b.b = t; } } });
    lay.skippedStereo = skipped;
  }

  /* ----------------------------------------------------------- layout */
  function layout(atoms, bonds, BL) {
    const n = atoms.length, adj = atoms.map(() => []);
    bonds.forEach((b, k) => { adj[b.a].push({ j: b.b, k }); adj[b.b].push({ j: b.a, k }); });
    const m0 = {
      atoms: atoms.map((a, i) => ({ id: i, el: a.el, charge: a.charge, x: 0, y: 0 })),
      bonds: bonds.map((b, k) => ({ id: k, a: b.a, b: b.b, order: b.order, stereo: 0 }))
    };
    const g = C.buildGraph(m0), rg = C.findRings(g), rs = bonds.map(() => 0);
    rg.rings.forEach(r => r.bonds.forEach(bi => { const k = g.bonds[bi].mi; if (!rs[k] || r.size < rs[k]) rs[k] = r.size; }));
    const pos = Array(n).fill(null), ts = Array(n).fill(1);
    function place(c, h) {
      const list = adj[c].slice().sort((p, q) => (rs[q.k] > 0) - (rs[p.k] > 0));
      let nr = 0, rr = 0;
      for (const e of list) {
        if (pos[e.j]) continue;
        let hh;
        if (rs[e.k]) { hh = h + (rr++ === 0 ? 1 : -1) * 2 * Math.PI / rs[e.k]; }
        else { const m = nr++; hh = h + (m === 0 ? ts[c] : m === 1 ? -ts[c] : 2) * Math.PI / 3; }
        ts[e.j] = -ts[c];
        pos[e.j] = [pos[c][0] + BL * Math.cos(hh), pos[c][1] + BL * Math.sin(hh)];
        place(e.j, hh);
      }
    }
    let ox = 0;
    for (let s0 = 0; s0 < n; s0++) {
      if (pos[s0]) continue;
      pos[s0] = [ox, 0]; place(s0, -Math.PI / 3);
      ox = Math.max.apply(null, pos.filter(Boolean).map(p => p[0])) + 2 * BL;
    }
    m0.atoms.forEach((a, i) => { a.x = pos[i][0]; a.y = pos[i][1]; });
    for (let t = 0; t < 2; t++) {
      try {
        const P = C.clean2D(m0, BL);
        if (P.every(p => isFinite(p[0]) && isFinite(p[1]))) m0.atoms.forEach((a, i) => { a.x = P[i][0]; a.y = P[i][1]; });
      } catch (e) { break; }
    }
    return { atoms: m0.atoms.map(a => ({ el: a.el, charge: a.charge, x: a.x, y: a.y })), bonds: bonds.map(b => ({ a: b.a, b: b.b, order: b.order, stereo: 0 })) };
  }

  function fromSMILES(str, BL) { const p = parse(str), lay = layout(p.atoms, p.bonds, BL || 40); assignWedges(p, lay); return lay; }

  /* ----------------------------------------------------------- export */
  function toSMILES(mol) {
    const g = C.buildGraph(mol), n = g.n;
    if (!n) return '';
    const rg = C.findRings(g), ar = C.aromaticity(mol, g, rg), aa = new Set(), ab = new Set();
    ar.rings.forEach(r => { if (r.type === 'aromatic') { r.atoms.forEach(a => aa.add(a)); r.bonds.forEach(b => ab.add(b)); } });
    const seen = Array(n).fill(false), tree = new Set(), done = new Set(), roots = [], clos = Array.from({ length: n }, () => []);
    function d1(u, pb) {
      seen[u] = true;
      for (const e of g.nbrs[u]) {
        if (e.bi === pb) continue;
        if (!seen[e.j]) { tree.add(e.bi); d1(e.j, e.bi); }
        else if (!tree.has(e.bi) && !done.has(e.bi)) { done.add(e.bi); clos[u].push(e.bi); clos[e.j].push(e.bi); }
      }
    }
    for (let s = 0; s < n; s++) if (!seen[s]) { roots.push(s); d1(s, -1); }
    const bsym = bi => {
      const b = g.bonds[bi];
      if (ab.has(bi)) return '';
      if (b.order === 2) return '='; if (b.order === 3) return '#';
      return (aa.has(b.i) && aa.has(b.j) && b.order === 1) ? '-' : '';
    };
    // stereocentres that carry a wedge / dash drawn from them get @ / @@
    const stereoAt = new Set();
    try { C.rsStereo(mol, g, rg).centres.forEach(c => { if (c.hasWedge) stereoAt.add(c.atom); }); } catch (e) { /* CIP centres only */ }
    // also any sp3 carbon with a wedge / dash drawn from it: pseudo-asymmetric ring carbons (cis / trans cyclohexanes) are not CIP centres but still need the tag
    g.bonds.forEach(b => { const u = b.i; if (b.stereo && mol.atoms[u].el === 'C' && !(mol.atoms[u].charge) && g.nbrs[u].length + g.hc[u] === 4 && g.nbrs[u].every(e => e.order === 1)) stereoAt.add(u); });
    const xy = mol.atoms.map(a => [a.x, a.y]);
    const atomTok = (i, tag) => {
      const el = mol.atoms[i].el, q = mol.atoms[i].charge || 0, h = g.hc[i], arom = aa.has(i), sy = arom ? el.toLowerCase() : el;
      let bare = false;
      if (tag) return '[' + sy + tag + (h ? 'H' + (h > 1 ? h : '') : '') + (q ? (q > 0 ? '+' : '-') + (Math.abs(q) > 1 ? Math.abs(q) : '') : '') + ']';
      if (q === 0 && DEFV[el]) {
        if (arom) bare = (el === 'C' || h === 0) && ['C', 'N', 'O', 'S', 'P', 'B'].indexOf(el) >= 0;
        else { const bos = g.bos[i], v = DEFV[el].find(x => x >= bos); bare = v !== undefined && v - bos === h; }
      }
      if (bare) return sy;
      return '[' + sy + (h ? 'H' + (h > 1 ? h : '') : '') + (q ? (q > 0 ? '+' : '-') + (Math.abs(q) > 1 ? Math.abs(q) : '') : '') + ']';
    };
    const num = new Map(), used = new Set();
    const lab = d => d > 9 ? '%' + d : String(d);
    function emit(u, pb) {
      let tag = '';
      if (stereoAt.has(u)) {
        // neighbours in the order they will appear in the string: parent, implicit H, ring-closure partners, then the branches
        const other = bi => (g.bonds[bi].i === u ? g.bonds[bi].j : g.bonds[bi].i), order = [];
        if (pb >= 0) order.push(other(pb));
        if (g.hc[u] === 1) order.push(-1);
        clos[u].forEach(bi => order.push(other(bi)));
        g.nbrs[u].filter(e => tree.has(e.bi) && e.bi !== pb).forEach(e => order.push(e.j));
        if (order.length === 4) {
          const z = new Map(); g.bonds.forEach(b => { if (b.stereo && b.i === u) z.set(b.j, b.stereo === 1 ? 1 : -1); });
          tag = chirTag(xy, u, order, z).tag;
        }
      }
      let s = atomTok(u, tag);
      clos[u].forEach(bi => {
        if (num.has(bi)) { const d = num.get(bi); s += lab(d); used.delete(d); num.delete(bi); }
        else { let d = 1; while (used.has(d)) d++; used.add(d); num.set(bi, d); s += bsym(bi) + lab(d); }
      });
      const kids = g.nbrs[u].filter(e => tree.has(e.bi) && e.bi !== pb);
      kids.forEach((e, k) => { const part = bsym(e.bi) + emit(e.j, e.bi); s += k < kids.length - 1 ? '(' + part + ')' : part; });
      return s;
    }
    return roots.map(r => emit(r, -1)).join('.');
  }

  C.fromSMILES = fromSMILES; C.toSMILES = toSMILES;
  if (typeof module !== 'undefined' && module.exports) module.exports = C;
})(typeof window !== 'undefined' ? window : globalThis);
