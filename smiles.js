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
    const atoms = [], bonds = [], stack = [], open = {};
    let prev = -1, sym = null, i = 0;
    const link = (a, b, sy) => bonds.push({ a, b, sy });
    const addAtom = a => { const k = atoms.length; atoms.push(a); if (prev >= 0) link(prev, k, sym); prev = k; sym = null; };
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
        if (open[d]) { const o = open[d]; delete open[d]; if (o.a === prev) throw new Error('ring closes on itself'); link(o.a, prev, sym || o.sy); }
        else open[d] = { a: prev, sy: sym };
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
        addAtom({ el: El, arom, charge: q, h: m[4] === undefined ? 0 : (m[4] === '' ? 1 : +m[4]), br: true });
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
    return { atoms, bonds };
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
    return { atoms: m0.atoms.map(a => ({ el: a.el, charge: a.charge, x: a.x, y: a.y })), bonds: bonds.map(b => ({ a: b.a, b: b.b, order: b.order })) };
  }

  function fromSMILES(str, BL) { const p = parse(str); return layout(p.atoms, p.bonds, BL || 40); }

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
    const atomTok = i => {
      const el = mol.atoms[i].el, q = mol.atoms[i].charge || 0, h = g.hc[i], arom = aa.has(i), sy = arom ? el.toLowerCase() : el;
      let bare = false;
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
      let s = atomTok(u);
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
