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
    const link = (a, b, sy) => bonds.push({ a, b, sy, dsym: (sy === '/' || sy === '\\') ? sy : '', dfrom: a });   // dsym / dfrom: direction of a / or \ mark (read as 'from dfrom towards the other atom')
    const addAtom = a => {
      if (sym && prev < 0) throw new Error('a bond symbol needs an atom on both sides');
      const k = atoms.length; atoms.push(a); nbo.push([]);
      if (prev >= 0) { link(prev, k, sym); nbo[prev].push(k); nbo[k].push(prev); }
      if (a.tag && a.h === 1) nbo[k].push(-1);                        // [C@H]: the H comes right after the preceding atom
      prev = k; sym = null;
    };
    while (i < s.length) {
      const ch = s[i];
      if (ch === '(') { if (prev < 0) throw new Error('branch without atom'); stack.push(prev); i++; continue; }
      if (ch === ')') { if (!stack.length) throw new Error('unbalanced )'); if (sym) throw new Error('bond symbol before )'); prev = stack.pop(); sym = null; i++; continue; }
      if ('-=#:/\\'.indexOf(ch) >= 0) { if (sym) throw new Error('two bond symbols in a row'); sym = ch; i++; continue; }
      if (ch === '.') { if (prev < 0 || sym) throw new Error('empty fragment or dangling bond before "."'); prev = -1; sym = null; i++; continue; }
      if (ch === '%' || /\d/.test(ch)) {
        let d; if (ch === '%') { d = s.substr(i + 1, 2); i += 3; } else { d = ch; i++; }
        if (!/^\d+$/.test(d)) throw new Error('bad ring label');
        if (prev < 0) throw new Error('ring digit without atom');
        if (open[d]) {
          const o = open[d]; delete open[d];
          if (o.a === prev) throw new Error('ring closes on itself');
          if (bonds.some(x => (x.a === o.a && x.b === prev) || (x.a === prev && x.b === o.a))) throw new Error('two bonds between the same pair of atoms');
          link(o.a, prev, sym || o.sy);
          const nb = bonds[bonds.length - 1], isDir = z => z === '/' || z === '\\';
          if (isDir(o.sy)) { nb.dsym = o.sy; nb.dfrom = o.a; } else if (isDir(sym)) { nb.dsym = sym; nb.dfrom = prev; } else nb.dsym = '';
          nbo[o.a][o.pos] = prev; nbo[prev].push(o.a);
        }
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
    if (sym) throw new Error('the string ends with a bond symbol');
    if (prev < 0) throw new Error('empty fragment after "."');
    if (stack.length) throw new Error('unbalanced (');
    if (Object.keys(open).length) throw new Error('unclosed ring digit');
    if (!atoms.length) throw new Error('no atoms');
    /* explicit hydrogens such as [H]C(F)(Cl)Br or [2H]: the editor keeps hydrogens implicit, so a terminal H atom is folded
       into its neighbour (its slot in the neighbour list becomes the implicit-H marker -1, so @ / @@ still reads correctly) */
    (function () {
      const deg = atoms.map(() => 0); bonds.forEach(b => { deg[b.a]++; deg[b.b]++; });
      const gone = atoms.map((a, i) => a.el === 'H' && !a.charge && deg[i] === 1);
      bonds.forEach(b => { if (gone[b.a] && gone[b.b]) { gone[b.a] = false; gone[b.b] = false; } });      // H-H stays
      if (!gone.some(Boolean)) return;
      bonds.forEach(b => {
        if (gone[b.a] || gone[b.b]) {
          const h = gone[b.a] ? b.a : b.b, o = h === b.a ? b.b : b.a;
          nbo[o] = nbo[o].map(x => x === h ? -1 : x);
          if (atoms[o].br) atoms[o].h++;
        }
      });
      const map = atoms.map(() => -1); let m = 0; atoms.forEach((a, i) => { if (!gone[i]) map[i] = m++; });
      const na = [], nn = [], nb = [];
      atoms.forEach((a, i) => { if (!gone[i]) { na.push(a); nn.push(nbo[i].map(x => x >= 0 ? map[x] : x)); } });
      bonds.forEach(b => { if (!gone[b.a] && !gone[b.b]) { b.dfrom = map[b.dfrom]; b.a = map[b.a]; b.b = map[b.b]; nb.push(b); } });
      atoms.length = 0; na.forEach(a => atoms.push(a)); bonds.length = 0; nb.forEach(b => bonds.push(b)); nbo.length = 0; nn.forEach(x => nbo.push(x));
    })();
    if (!atoms.length) throw new Error('no atoms');
    // an aromatic bond is only a ring bond: between two separate aromatic rings (biphenyl, c1ccccc1c1ccccc1) it is single
    (function () {
      const adjB = atoms.map(() => []); bonds.forEach((b, k) => { adjB[b.a].push([b.b, k]); adjB[b.b].push([b.a, k]); });
      bonds.forEach((b, k) => {
        if (b.sy === ':' || (!b.sy && atoms[b.a].arom && atoms[b.b].arom)) {
          const seen = new Set([b.a]), q = [b.a]; let ring = false;
          for (let h = 0; h < q.length && !ring; h++) for (const [v, kk] of adjB[q[h]]) { if (kk === k || seen.has(v)) continue; if (v === b.b) { ring = true; break; } seen.add(v); q.push(v); }
          b.inRing = ring;
        }
      });
    })();
    bonds.forEach(b => {
      b.order = b.sy === '=' ? 2 : b.sy === '#' ? 3 : b.sy === ':' ? (b.inRing ? 'ar' : 1) : b.sy ? 1 : (atoms[b.a].arom && atoms[b.b].arom && b.inRing ? 'ar' : 1);
    });
    atoms.forEach((a, i) => { if (a.arom && !bonds.some(b => (b.a === i || b.b === i) && b.order === 'ar')) throw new Error('aromatic atom ' + a.el.toLowerCase() + ' is not in an aromatic ring'); });
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
    // neutral pentavalent N written as N(=O)=O / C[N](C)(C)=O: draw it charge-separated, N+ and O-, like PubChem does
    atoms.forEach((a, i) => {
      if (a.el !== 'N' || a.charge) return;
      const inc = bonds.filter(b => b.a === i || b.b === i);
      if (inc.reduce((t, b) => t + b.order, 0) !== 5) return;
      const ob = inc.find(b => b.order === 2 && (() => { const o = atoms[b.a === i ? b.b : b.a]; return o.el === 'O' && !o.charge && bonds.filter(x => x.a === (b.a === i ? b.b : b.a) || x.b === (b.a === i ? b.b : b.a)).length === 1; })());
      if (!ob) return;
      ob.order = 1; atoms[ob.a === i ? ob.b : ob.a].charge = -1; a.charge = 1;
    });
    // valence check: refuse structures the editor would silently draw differently (wrong hydrogen counts, impossible atoms)
    atoms.forEach((a, i) => {
      const bos = nbr[i].reduce((t, k) => t + bonds[k].order, 0), vals = C.valences(a.el, a.charge);
      const v = vals.find(x => x >= bos);
      const name = a.br ? '[' + (a.arom ? a.el.toLowerCase() : a.el) + (a.h ? 'H' + (a.h > 1 ? a.h : '') : '') + (a.charge ? (a.charge > 0 ? '+' : '-') + (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') : '') + ']' : a.el;
      if (v === undefined) throw new Error(name + ' has too many bonds (' + bos + ') for its valence');
      if (a.br && v - bos !== a.h) throw new Error(name + ' is written with ' + a.h + ' hydrogen' + (a.h === 1 ? '' : 's') + ' but the editor would give it ' + (v - bos) + ' (radicals and unusual charge states are not supported)');
    });
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
        enter(e.j, hh);
      }
    }
    /* ---- ring systems: every ring (single or fused) is laid out as a regular polygon, fused rings sharing an edge
       with their neighbour. The turn-by-turn walk used before cannot close a ring on atoms that are already placed
       (it mangled naphthalene, and any ring entered from a chain such as the benzene in phenetole). A ring system is
       placed as one block the moment the walk first reaches one of its atoms. */
    const rings = rg.rings, rc = rings.map(() => null);          // rc[r] = polygon centre once ring r is placed
    const nRing = rings.length, grp = rings.map((_, i) => i);
    const find = x => grp[x] === x ? x : (grp[x] = find(grp[x]));
    for (let a = 0; a < nRing; a++) for (let b = a + 1; b < nRing; b++)
      if (rings[a].atoms.filter(x => rings[b].atoms.indexOf(x) >= 0).length >= 2) grp[find(a)] = find(b);
    const clusters = {}, clusterOf = Array(n).fill(undefined), cdone = {};
    rings.forEach((r, i) => { const key = find(i); (clusters[key] = clusters[key] || []).push(i); r.atoms.forEach(at => { if (clusterOf[at] === undefined) clusterOf[at] = key; }); });
    const radius = k => BL / (2 * Math.sin(Math.PI / k));
    // walk ring order from placed atom `a` (which is at angle `ang` round `cen`) filling in the polygon
    function polygon(ring, a, cen, ang) {
      const k = ring.size, ia = ring.atoms.indexOf(a), R = radius(k);
      let sgn = 1;
      const nb = ring.atoms[(ia + 1) % k];
      if (pos[nb]) { const cr = (pos[a][0] - cen[0]) * (pos[nb][1] - cen[1]) - (pos[a][1] - cen[1]) * (pos[nb][0] - cen[0]); sgn = cr >= 0 ? 1 : -1; }
      else { const pv = ring.atoms[(ia + k - 1) % k]; if (pos[pv]) { const cr = (pos[a][0] - cen[0]) * (pos[pv][1] - cen[1]) - (pos[a][1] - cen[1]) * (pos[pv][0] - cen[0]); sgn = cr >= 0 ? -1 : 1; } }
      for (let t = 0; t < k; t++) {
        const at = ring.atoms[(ia + t) % k];
        if (pos[at]) continue;
        const th = ang + sgn * 2 * Math.PI * t / k;
        pos[at] = [cen[0] + R * Math.cos(th), cen[1] + R * Math.sin(th)];
      }
    }
    function placeCluster(key, a, cen) {
      cdone[key] = true;
      const cl = clusters[key], ri0 = cl.find(r => rings[r].atoms.indexOf(a) >= 0);
      rc[ri0] = cen;
      polygon(rings[ri0], a, cen, Math.atan2(pos[a][1] - cen[1], pos[a][0] - cen[0]));
      let progress = true;
      while (progress) {
        progress = false;
        for (const ri of cl) {
          if (rc[ri]) continue;
          const ring = rings[ri], k = ring.size;
          let pair = null;                                         // an edge of this ring whose two atoms are already placed
          for (let t = 0; t < k && !pair; t++) { const x = ring.atoms[t], y = ring.atoms[(t + 1) % k]; if (pos[x] && pos[y]) pair = [x, y]; }
          if (!pair) continue;
          let host = -1;                                           // placed ring that owns that edge
          cl.forEach(rj => { if (rc[rj] && rings[rj].atoms.indexOf(pair[0]) >= 0 && rings[rj].atoms.indexOf(pair[1]) >= 0) host = rj; });
          if (host < 0) continue;
          const mx = (pos[pair[0]][0] + pos[pair[1]][0]) / 2, my = (pos[pair[0]][1] + pos[pair[1]][1]) / 2;
          let nx = -(pos[pair[1]][1] - pos[pair[0]][1]), ny = pos[pair[1]][0] - pos[pair[0]][0];
          const nl = Math.hypot(nx, ny) || 1; nx /= nl; ny /= nl;
          if ((mx + nx - rc[host][0]) ** 2 + (my + ny - rc[host][1]) ** 2 < (mx - nx - rc[host][0]) ** 2 + (my - ny - rc[host][1]) ** 2) { nx = -nx; ny = -ny; }   // away from the host ring
          const ap = BL / (2 * Math.tan(Math.PI / k));
          rc[ri] = [mx + nx * ap, my + ny * ap];
          polygon(ring, pair[0], rc[ri], Math.atan2(pos[pair[0]][1] - rc[ri][1], pos[pair[0]][0] - rc[ri][0]));
          progress = true;
        }
      }
      // (a ring of the cluster that cannot be reached this way, e.g. in a bridged system, is left to the relaxer)
      const inCl = new Set(); cl.forEach(ri => rings[ri].atoms.forEach(x => inCl.add(x)));
      inCl.forEach(c => {
        if (!pos[c]) return;
        const free = adj[c].filter(e => !pos[e.j]);
        if (!free.length) return;
        let sx = 0, sy = 0;
        adj[c].forEach(e => { if (pos[e.j]) { const l = Math.hypot(pos[e.j][0] - pos[c][0], pos[e.j][1] - pos[c][1]) || 1; sx += (pos[e.j][0] - pos[c][0]) / l; sy += (pos[e.j][1] - pos[c][1]) / l; } });
        const h = Math.atan2(-sy, -sx);
        free.forEach((e, m) => {
          if (pos[e.j]) return;
          const hh = free.length === 1 ? h : h + (m === 0 ? 1 : m === 1 ? -1 : 0) * Math.PI / 3;
          ts[e.j] = 1;
          pos[e.j] = [pos[c][0] + BL * Math.cos(hh), pos[c][1] + BL * Math.sin(hh)];
          enter(e.j, hh);
        });
      });
    }
    // called whenever the walk has just put atom j at its position, arriving with heading hh
    function enter(j, hh) {
      const key = clusterOf[j];
      if (key !== undefined && !cdone[key]) {
        let placed = 0; clusters[key].forEach(ri => rings[ri].atoms.forEach(x => { if (pos[x]) placed++; }));
        if (placed === 1) {                                        // only j itself: the whole ring system can be dropped in, pointing away from the parent
          const R = radius(rings[clusters[key].find(r => rings[r].atoms.indexOf(j) >= 0)].size);
          placeCluster(key, j, [pos[j][0] + R * Math.cos(hh), pos[j][1] + R * Math.sin(hh)]);
          return;
        }
      }
      if (key !== undefined && !cdone[key]) {                      // spiro: the walk arrived through the one atom that two ring systems share
        const plc = []; clusters[key].forEach(ri => rings[ri].atoms.forEach(x => { if (pos[x] && plc.indexOf(x) < 0) plc.push(x); }));
        const sp = plc.length === 2 ? plc.find(x => x !== j) : undefined;
        if (sp !== undefined && adj[j].some(e => e.j === sp)) {
          const ri0 = clusters[key].find(r => rings[r].atoms.indexOf(j) >= 0), k = rings[ri0].size;
          const mx = (pos[j][0] + pos[sp][0]) / 2, my = (pos[j][1] + pos[sp][1]) / 2;
          let nx = -(pos[sp][1] - pos[j][1]), ny = pos[sp][0] - pos[j][0]; const nl = Math.hypot(nx, ny) || 1; nx /= nl; ny /= nl;
          let cx = 0, cy = 0, cnt = 0; adj[sp].forEach(e => { if (pos[e.j] && e.j !== j) { cx += pos[e.j][0]; cy += pos[e.j][1]; cnt++; } });
          if (cnt) { cx /= cnt; cy /= cnt; if ((mx + nx - cx) ** 2 + (my + ny - cy) ** 2 < (mx - nx - cx) ** 2 + (my - ny - cy) ** 2) { nx = -nx; ny = -ny; } }
          const ap = BL / (2 * Math.tan(Math.PI / k));
          placeCluster(key, j, [mx + nx * ap, my + ny * ap]);
          return;
        }
      }
      place(j, hh);
    }
    let ox = 0;
    for (let s0 = 0; s0 < n; s0++) {
      if (pos[s0]) continue;
      pos[s0] = [ox, 0]; enter(s0, -Math.PI / 3);
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

  /* ------------------------------------------- double-bond geometry from / and \ marks
     F/C=C/F is trans (E), F/C=C\\F is cis (Z). For every acyclic double bond that has a / or \ mark on both ends, the drawing is
     mirrored across the double bond (the side of the molecule beyond atom B) whenever the layout came out with the wrong cis / trans. */
  function applyDoubleBonds(p, lay) {
    const n = lay.atoms.length, adj = lay.atoms.map(() => []);
    lay.bonds.forEach((b, k) => { adj[b.a].push([b.b, k]); adj[b.b].push([b.a, k]); });
    const inCycle = (a, b, k) => {                                  // is bond k (a-b) part of a ring?
      const seen = new Set([a]), q = [a];
      for (let h = 0; h < q.length; h++) for (const [v, kk] of adj[q[h]]) { if (kk === k || seen.has(v)) continue; if (v === b) return true; seen.add(v); q.push(v); }
      return false;
    };
    const side = (X, k) => {                                         // +1 / -1 = the substituent on bond k lies "above" / "below" X
      const b = p.bonds[k], s = b.a === X ? b.b : b.a;
      return b.dfrom === X ? (b.dsym === '/' ? 1 : -1) : (b.dsym === '/' ? -1 : 1);
    };
    const mark = X => { for (const [v, k] of adj[X]) if (p.bonds[k].dsym && p.bonds[k].order === 1) return { v, k }; return null; };
    let applied = 0;
    lay.bonds.forEach((bd, kb) => {
      if (bd.order !== 2) return;
      const A = bd.a, B = bd.b;
      const mA = (() => { for (const [v, k] of adj[A]) if (k !== kb && p.bonds[k].dsym && p.bonds[k].order === 1) return { v, k }; return null; })();
      const mB = (() => { for (const [v, k] of adj[B]) if (k !== kb && p.bonds[k].dsym && p.bonds[k].order === 1) return { v, k }; return null; })();
      if (!mA || !mB || inCycle(A, B, kb)) return;
      const wantCis = side(A, mA.k) === side(B, mB.k);
      const ax = lay.atoms[B].x - lay.atoms[A].x, ay = lay.atoms[B].y - lay.atoms[A].y;
      const cr = (v, O) => ax * (lay.atoms[v].y - lay.atoms[O].y) - ay * (lay.atoms[v].x - lay.atoms[O].x);
      const isCis = (cr(mA.v, A) > 0) === (cr(mB.v, B) > 0);
      if (isCis === wantCis) return;
      const comp = new Set([B]), q = [B];                             // everything beyond B
      for (let h = 0; h < q.length; h++) for (const [v, k] of adj[q[h]]) { if (k === kb || comp.has(v)) continue; comp.add(v); q.push(v); }
      const L = Math.hypot(ax, ay) || 1, ux = ax / L, uy = ay / L, A0 = lay.atoms[A];
      comp.forEach(i => {
        const a = lay.atoms[i], dx = a.x - A0.x, dy = a.y - A0.y, along = dx * ux + dy * uy;
        a.x = A0.x + 2 * along * ux - dx; a.y = A0.y + 2 * along * uy - dy;      // mirror image across the line through A and B
      });
      applied++;
    });
    return applied;
  }

  function fromSMILES(str, BL) { const p = parse(str), lay = layout(p.atoms, p.bonds, BL || 40); applyDoubleBonds(p, lay); assignWedges(p, lay); return lay; }

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
