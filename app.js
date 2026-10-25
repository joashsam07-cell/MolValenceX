(() => {
'use strict';
const C = window.Chem, $ = s => document.querySelector(s), BL = 40;
const cv = $('#ed'), cx = cv.getContext('2d');
let mol = { atoms: [], bonds: [] }, nid = 1, ver = 0, tool = 'bond', vr = { o: 1, s: 0 }, el = 'C';
let sel = null, drag = null, an = null, hist = [], hi = -1, cache = null;

/* ---------- state / history ---------- */
const snap = () => JSON.stringify({ m: mol, n: nid });
function commit() { hist = hist.slice(0, hi + 1); hist.push(snap()); hi = hist.length - 1; changed(); }
function changed() { ver++; an = null; update(); draw(); }
function restore() { const s = JSON.parse(hist[hi]); mol = s.m; nid = s.n; sel = null; changed(); }
function undo() { if (hi > 0) { hi--; restore(); } }
function redo() { if (hi < hist.length - 1) { hi++; restore(); } }
function G() { if (!cache || cache.v !== ver) { const g = C.buildGraph(mol); cache = { v: ver, g, rg: C.findRings(g) }; } return cache; }
const atomById = id => mol.atoms.find(a => a.id === id);

/* ---------- canvas sizing ---------- */
function size() {
  const d = window.devicePixelRatio || 1, r = cv.getBoundingClientRect();
  if(!r.width || !r.height) return; // Prevent crashes on hidden canvas
  cv.width = r.width * d; cv.height = r.height * d; cx.setTransform(d, 0, 0, d, 0, 0); draw();
}
window.addEventListener('resize', () => { size(); size3(); });

/* ---------- hit tests ---------- */
const atomAt = p => mol.atoms.find(a => Math.hypot(a.x - p.x, a.y - p.y) < 14);
function bondAt(p) {
  let best = null, bd = 7;
  mol.bonds.forEach(b => {
    const A = atomById(b.a), B = atomById(b.b), dx = B.x - A.x, dy = B.y - A.y, l2 = dx * dx + dy * dy || 1;
    const t = Math.max(0, Math.min(1, ((p.x - A.x) * dx + (p.y - A.y) * dy) / l2));
    const d = Math.hypot(A.x + t * dx - p.x, A.y + t * dy - p.y);
    if (d < bd) { bd = d; best = b; }
  });
  return best;
}

/* ---------- editing ---------- */
function addAtom(e, x, y) { const a = { id: nid++, el: e, x, y, charge: 0 }; mol.atoms.push(a); ver++; an = null; return a; }
function setBond(ex, v, from) {
  if (v.s) {
    if (ex.stereo === v.s && ex.a === from.id) { [ex.a, ex.b] = [ex.b, ex.a]; }
    else { ex.stereo = v.s; ex.order = 1; if (ex.a !== from.id) [ex.a, ex.b] = [ex.b, ex.a]; }
  } else { ex.order = v.o; ex.stereo = 0; }
}
function addBond(a, b, v) {
  const ex = mol.bonds.find(k => (k.a === a.id && k.b === b.id) || (k.a === b.id && k.b === a.id));
  if (ex) setBond(ex, v, a);
  else mol.bonds.push({ id: nid++, a: a.id, b: b.id, order: v.s ? 1 : v.o, stereo: v.s });
}
function freeDir(s) {
  const nb = mol.bonds.filter(b => b.a === s.id || b.b === s.id).map(b => atomById(b.a === s.id ? b.b : b.a));
  const ang = nb.map(n => Math.atan2(n.y - s.y, n.x - s.x));
  const pos = t => ({ x: s.x + BL * Math.cos(t), y: s.y + BL * Math.sin(t) });
  const room = t => { const q = pos(t); return Math.min(1e9, ...mol.atoms.filter(a => a !== s).map(a => Math.hypot(a.x - q.x, a.y - q.y))); };
  if (!ang.length) return pos(-Math.PI / 6);
  if (ang.length === 1) { const c = [ang[0] + 2 * Math.PI / 3, ang[0] - 2 * Math.PI / 3]; return pos(room(c[0]) >= room(c[1]) ? c[0] : c[1]); }
  const sa = ang.slice().sort((a, b) => a - b); let bg = -1, bt = 0;
  sa.forEach((a, i) => { const nx = i + 1 < sa.length ? sa[i + 1] : sa[0] + 2 * Math.PI, g = nx - a; if (g > bg) { bg = g; bt = a + g / 2; } });
  return pos(bt);
}
const P = e => { const r = cv.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
function bondEnd(d, p) {
  const t = atomAt(p);
  if (t && t !== d.s) return { atom: t, x: t.x, y: t.y };
  const dx = p.x - d.s.x, dy = p.y - d.s.y;
  if (Math.hypot(dx, dy) < 10) return null;
  const a = Math.round(Math.atan2(dy, dx) / (Math.PI / 6)) * (Math.PI / 6);
  return { x: d.s.x + BL * Math.cos(a), y: d.s.y + BL * Math.sin(a) };
}
cv.onpointerdown = e => {
  cv.setPointerCapture(e.pointerId);
  const p = P(e), a = atomAt(p), b = a ? null : bondAt(p);
  if (tool === 'preset') { placePreset(curPre, p, a, b); return; }
  if (tool === 'erase') {
    if (a) { mol.atoms = mol.atoms.filter(x => x !== a); mol.bonds = mol.bonds.filter(k => k.a !== a.id && k.b !== a.id); commit(); }
    else if (b) { mol.bonds = mol.bonds.filter(k => k !== b); commit(); }
    return;
  }
  if (tool === 'chg+' || tool === 'chg-') { if (a) { a.charge = (a.charge || 0) + (tool === 'chg+' ? 1 : -1); commit(); } return; }
  if (tool === 'select') {
    sel = a ? a.id : null;
    drag = { m: a ? 'move' : 'pan', a, x: p.x, y: p.y, moved: false }; draw(); return;
  }
  if (!a && b) {
    const A = atomById(b.a), B = atomById(b.b);
    setBond(b, vr, Math.hypot(A.x - p.x, A.y - p.y) <= Math.hypot(B.x - p.x, B.y - p.y) ? A : B);
    commit(); return;
  }
  const fresh = !a, s = a || addAtom(el, p.x, p.y);
  drag = { m: 'bond', s, fresh, x: p.x, y: p.y, moved: false, end: null };
  draw();
};
cv.onpointermove = e => {
  if (!drag) return;
  const p = P(e), dx = p.x - drag.x, dy = p.y - drag.y;
  if (Math.hypot(dx, dy) > 4) drag.moved = true;
  if (drag.m === 'move') { drag.a.x += dx; drag.a.y += dy; drag.x = p.x; drag.y = p.y; ver++; an = null; }
  else if (drag.m === 'pan') { mol.atoms.forEach(a => { a.x += dx; a.y += dy; }); drag.x = p.x; drag.y = p.y; }
  else drag.end = bondEnd(drag, p);
  draw();
};
cv.onpointerup = e => {
  if (!drag) return;
  const d = drag; drag = null;
  if (d.m === 'move' || d.m === 'pan') { if (d.moved) commit(); else draw(); return; }
  if (d.end) {
    const t = d.end.atom || addAtom(el, d.end.x, d.end.y);
    addBond(d.s, t, vr); commit(); return;
  }
  if (d.fresh) { commit(); return; }
  if (tool === 'atom') { if (d.s.el !== el) { d.s.el = el; commit(); } else draw(); return; }
  const q = freeDir(d.s), n = addAtom(el, q.x, q.y); addBond(d.s, n, vr); commit();
};

/* ---------- drawing ---------- */
const hasLabel = (a, i, g) => a.el !== 'C' || !g.nbrs[i].length || a.charge;
const cen = ids => { let x = 0, y = 0; ids.forEach(i => { x += mol.atoms[i].x; y += mol.atoms[i].y; }); return [x / ids.length, y / ids.length]; };
function line(x1, y1, x2, y2, w) { cx.lineWidth = w || 2; cx.beginPath(); cx.moveTo(x1, y1); cx.lineTo(x2, y2); cx.stroke(); }
function draw() {
  const W = cv.clientWidth, H = cv.clientHeight; cx.clearRect(0, 0, W, H);
  const { g, rg } = G(), A = mol.atoms;
  if (an) an.rings.forEach(r => {
    const col = { aromatic: '34,197,94', antiaromatic: '239,68,68', nonaromatic: '148,163,184' }[r.type];
    cx.beginPath(); r.atoms.forEach((a, k) => k ? cx.lineTo(A[a].x, A[a].y) : cx.moveTo(A[a].x, A[a].y)); cx.closePath();
    cx.fillStyle = 'rgba(' + col + ',.2)'; cx.fill();
    const c = cen(r.atoms); cx.fillStyle = 'rgb(' + col + ')'; cx.font = 'bold 12px system-ui'; cx.textAlign = 'center'; cx.fillText(r.e + 'π', c[0], c[1] + 4);
  });
  cx.strokeStyle = '#222'; cx.fillStyle = '#222'; cx.lineCap = 'round';
  mol.bonds.forEach((b, bi) => {
    const i = g.idx.get(b.a), j = g.idx.get(b.b), a = A[i], c = A[j];
    let dx = c.x - a.x, dy = c.y - a.y; const l = Math.hypot(dx, dy) || 1; dx /= l; dy /= l;
    const g1 = hasLabel(a, i, g) ? 11 : 0, g2 = hasLabel(c, j, g) ? 11 : 0;
    const x1 = a.x + dx * g1, y1 = a.y + dy * g1, x2 = c.x - dx * g2, y2 = c.y - dy * g2, nx = -dy, ny = dx;
    if (b.stereo === 1) { cx.beginPath(); cx.moveTo(x1, y1); cx.lineTo(x2 + nx * 5, y2 + ny * 5); cx.lineTo(x2 - nx * 5, y2 - ny * 5); cx.closePath(); cx.fill(); }
    else if (b.stereo === 2) { for (let k = 1; k <= 7; k++) { const t = k / 7.5, w = 1 + 5 * t, px = x1 + (x2 - x1) * t, py = y1 + (y2 - y1) * t; line(px - nx * w, py - ny * w, px + nx * w, py + ny * w, 1.6); } }
    else if (b.order === 1) line(x1, y1, x2, y2);
    else if (b.order === 3) { line(x1, y1, x2, y2); line(x1 + nx * 4, y1 + ny * 4, x2 + nx * 4, y2 + ny * 4); line(x1 - nx * 4, y1 - ny * 4, x2 - nx * 4, y2 - ny * 4); }
    else {
      const r = rg.rings.filter(r => r.bonds.includes(bi)).sort((p, q) => p.size - q.size)[0];
      if (r) {
        const c0 = cen(r.atoms), s = Math.sign(nx * (c0[0] - (a.x + c.x) / 2) + ny * (c0[1] - (a.y + c.y) / 2)) || 1;
        line(x1, y1, x2, y2); line(x1 + nx * 5 * s + dx * 6, y1 + ny * 5 * s + dy * 6, x2 + nx * 5 * s - dx * 6, y2 + ny * 5 * s - dy * 6);
      } else { line(x1 + nx * 2.5, y1 + ny * 2.5, x2 + nx * 2.5, y2 + ny * 2.5); line(x1 - nx * 2.5, y1 - ny * 2.5, x2 - nx * 2.5, y2 - ny * 2.5); }
    }
  });
  A.forEach((a, i) => {
    if (a.id === sel) { cx.strokeStyle = '#6366f1'; cx.lineWidth = 2; cx.beginPath(); cx.arc(a.x, a.y, 15, 0, 7); cx.stroke(); }
    if (!g.ok[i]) { cx.strokeStyle = '#ef4444'; cx.setLineDash([3, 3]); cx.lineWidth = 2; cx.beginPath(); cx.arc(a.x, a.y, 16, 0, 7); cx.stroke(); cx.setLineDash([]); }
    if (!hasLabel(a, i, g)) return;
    const h = g.hc[i], col = C.EL[a.el] ? C.EL[a.el].c : '#111';
    cx.font = 'bold 16px system-ui'; cx.textAlign = 'center'; cx.textBaseline = 'middle';
    const sw = cx.measureText(a.el).width, hw = cx.measureText('H').width;
    let chx = a.x + sw / 2 + 7;
    cx.fillStyle = '#fff'; cx.beginPath(); cx.arc(a.x, a.y, 9, 0, 7); cx.fill();
    cx.fillStyle = col; cx.fillText(a.el, a.x, a.y + 1);
    if (h > 0) {
      const right = g.nbrs[i].reduce((s, e) => s + A[e.j].x - a.x, 0) <= 0;
      const hx = a.x + (right ? 1 : -1) * (sw / 2 + hw / 2 + 1);
      cx.fillText('H', hx, a.y + 1);
      if (right) chx = hx + hw / 2 + (h > 1 ? 14 : 6);
      if (h > 1) { cx.font = 'bold 11px system-ui'; cx.fillText(h, hx + (right ? 1 : -1) * (hw / 2 + 5), a.y + 6); }
    }
    if (a.charge) { cx.font = 'bold 11px system-ui'; cx.fillStyle = col; cx.fillText((Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') + (a.charge > 0 ? '+' : '−'), chx, a.y - 8); }
    cx.textBaseline = 'alphabetic';
  });
  if (drag && drag.m === 'bond' && drag.end) {
    cx.strokeStyle = '#6366f1'; cx.setLineDash([5, 4]); line(drag.s.x, drag.s.y, drag.end.x, drag.end.y, 2); cx.setLineDash([]);
  }
}

/* ---------- info panel ---------- */
function update() {
  const f = C.formula(mol);
  $('#formula').innerHTML = mol.atoms.length ? f.formulaHtml : '–';
  $('#mass').textContent = mol.atoms.length ? f.mass.toFixed(3) + ' g/mol' : '–';
  $('#dbe').textContent = f.dbe === null ? '–' : f.dbe;
  $('#cnt').textContent = mol.atoms.length ? f.heavy + ' / ' + f.hydrogens : '–';
  $('#result').innerHTML = '';
}
$('#arom').onclick = () => {
  const { g, rg } = G();
  if (!mol.atoms.length) { $('#result').innerHTML = '<div class="warn">Draw a molecule first.</div>'; return; }
  an = C.aromaticity(mol, g, rg);
  const bad = mol.atoms.map((a, i) => g.ok[i] ? '' : a.el + (i + 1)).filter(Boolean);
  const names = { aromatic: 'Aromatic', antiaromatic: 'Anti-aromatic', nonaromatic: 'Non-aromatic', mixed: 'Aromatic + anti-aromatic rings', acyclic: 'No ring' };
  let h = (bad.length ? '<div class="warn">Over-valent atom(s): ' + bad.join(', ') + ' — results may be unreliable.</div>' : '');
  h += '<span class="chip ' + an.overall + '">' + names[an.overall] + '</span><div>' + an.message + '</div>';
  an.rings.forEach((r, k) => {
    h += '<div class="ring ' + r.type + '"><b>Ring ' + (k + 1) + ' · ' + r.size + '-membered · ' + r.e + ' π e⁻ · ' + names[r.type] + '</b><br>' + r.reason +
      '<small>' + r.contrib.map(c => c.label + ':' + c.e).join('  ') + '</small></div>';
  });
  $('#result').innerHTML = h; draw();
};

/* ---------- toolbar ---------- */
document.querySelectorAll('#bondtools button').forEach(b => b.onclick = () => {
  document.querySelectorAll('#bondtools button, .pre').forEach(x => x.classList.remove('on')); b.classList.add('on');
  $('#status').textContent = 'New atoms: ' + el;
  tool = b.dataset.tool; if (tool === 'bond') vr = { o: +b.dataset.o, s: +b.dataset.s };
});
function pickTool(sel_) { const b = document.querySelector('#bondtools button' + sel_); if (b) b.click(); }
const ELS = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Se'];
function setEl(e) {
  el = e; $('#status').textContent = 'New atoms: ' + e;
  document.querySelectorAll('#elements button').forEach(x => x.classList.toggle('on', x.textContent === e));
  if (tool !== 'bond' && tool !== 'atom') { document.querySelectorAll('#bondtools button').forEach(x => x.classList.remove('on')); }
  tool = 'atom'; document.querySelectorAll('#bondtools button, .pre').forEach(x => x.classList.remove('on'));
}
ELS.forEach(e => { const b = document.createElement('button'); b.textContent = e; b.style.color = C.EL[e].c; b.onclick = () => setEl(e); $('#elements').appendChild(b); });
['Li', 'Na', 'K', 'Mg', 'Ca', 'Al'].forEach(e => { const o = document.createElement('option'); o.textContent = e; $('#more').appendChild(o); });
$('#more').onchange = e => { if (e.target.value) setEl(e.target.value); e.target.value = ''; };
$('#undo').onclick = undo; $('#redo').onclick = redo;
$('#clear').onclick = () => { mol = { atoms: [], bonds: [] }; sel = null; commit(); };
$('#clean').onclick = () => { const p = C.clean2D(mol, BL); mol.atoms.forEach((a, i) => { a.x = p[i][0]; a.y = p[i][1]; }); commit(); };
const dl = (name, text) => { const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text])); a.download = name; a.click(); };
$('#expmol').onclick = () => dl('molecule.mol', C.toMolBlock(mol, BL));
window.addEventListener('keydown', e => {
  if (e.ctrlKey && e.key === 'z') { e.preventDefault(); undo(); return; }
  if (e.ctrlKey && e.key === 'y') { e.preventDefault(); redo(); return; }
  if (e.key === 'Delete' || e.key === 'Backspace') {
    if (sel) { mol.atoms = mol.atoms.filter(a => a.id !== sel); mol.bonds = mol.bonds.filter(k => k.a !== sel && k.b !== sel); sel = null; commit(); } return;
  }
  if (e.key === 'Escape' && tool === 'preset') { pickTool('[data-o="1"][data-s="0"]'); return; }
  const m = { c: 'C', n: 'N', o: 'O', s: 'S', p: 'P', f: 'F', h: 'H', b: 'B', i: 'I', l: 'Cl', r: 'Br' }[e.key];
  if (m && !e.ctrlKey) setEl(m);
  if (e.key === '1') pickTool('[data-o="1"][data-s="0"]'); if (e.key === '2') pickTool('[data-o="2"]'); if (e.key === '3') pickTool('[data-o="3"]');
});

/* ---------- examples ---------- */
function B() { const m = { atoms: [], bonds: [] }; return { m, a(e, x, y, q) { m.atoms.push({ id: nid, el: e, x, y, charge: q || 0 }); return nid++; }, b(i, j, o, s) { if (!m.bonds.some(k => (k.a === i && k.b === j) || (k.a === j && k.b === i))) m.bonds.push({ id: nid++, a: i, b: j, order: o || 1, stereo: s || 0 }); }, pos(i) { return m.atoms.find(a => a.id === i); } }; }
function ring(b, x, y, n, els, ord) {
  const r = BL / (2 * Math.sin(Math.PI / n)), ids = [];
  for (let k = 0; k < n; k++) { const t = -Math.PI / 2 + 2 * Math.PI * k / n, e = els[k] || 'C', q = Array.isArray(e); ids.push(b.a(q ? e[0] : e, x + r * Math.cos(t), y + r * Math.sin(t), q ? e[1] : 0)); }
  ids.forEach((id, k) => b.b(id, ids[(k + 1) % n], ord[k])); return ids;
}
function fuse(b, ids0, i, j, n, els, ord, away) {
  const A = b.pos(i), Bp = b.pos(j), mx = (A.x + Bp.x) / 2, my = (A.y + Bp.y) / 2;
  let nx = -(Bp.y - A.y), ny = Bp.x - A.x; const l = Math.hypot(nx, ny); nx /= l; ny /= l;
  if (nx * (mx - away[0]) + ny * (my - away[1]) < 0) { nx = -nx; ny = -ny; }
  const ap = BL / (2 * Math.tan(Math.PI / n)), r = BL / (2 * Math.sin(Math.PI / n)), c = [mx + nx * ap, my + ny * ap];
  const ai = Math.atan2(A.y - c[1], A.x - c[0]); let d = Math.atan2(Bp.y - c[1], Bp.x - c[0]) - ai; while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI;
  const sg = d > 0 ? 1 : -1, ids = [i, j];
  for (let k = 2; k < n; k++) { const t = ai + sg * k * 2 * Math.PI / n; ids.push(b.a(els[k - 2] || 'C', c[0] + r * Math.cos(t), c[1] + r * Math.sin(t))); }
  ids.forEach((id, k) => b.b(id, ids[(k + 1) % n], ord[k])); return ids;
}
const K6 = [2, 1, 2, 1, 2, 1], K5 = [1, 2, 1, 2, 1], K7 = [1, 2, 1, 2, 1, 2, 1];
const EX = {
  'Benzene': b => ring(b, 300, 250, 6, [], K6),
  'Naphthalene': b => { const a = ring(b, 260, 250, 6, [], K6); fuse(b, a, a[1], a[2], 6, [], [1, 1, 2, 1, 2, 1], cen2(b, a)); },
  'Azulene': b => { const a = ring(b, 260, 250, 7, [], [1, 2, 1, 2, 1, 2, 1]); fuse(b, a, a[6], a[0], 5, [], [1, 2, 1, 2, 1], cen2(b, a)); },
  'Pyridine': b => ring(b, 300, 250, 6, ['N'], K6),
  'Pyrrole': b => ring(b, 300, 250, 5, ['N'], K5),
  'Furan': b => ring(b, 300, 250, 5, ['O'], K5),
  'Thiophene': b => ring(b, 300, 250, 5, ['S'], K5),
  'Cyclobutadiene': b => ring(b, 300, 250, 4, [], [2, 1, 2, 1]),
  'Cyclooctatetraene': b => ring(b, 300, 250, 8, [], [2, 1, 2, 1, 2, 1, 2, 1]),
  'Cyclopentadienyl anion': b => ring(b, 300, 250, 5, [['C', -1]], K5),
  'Tropylium cation': b => ring(b, 300, 250, 7, [['C', 1]], K7),
  'Cyclopentadiene': b => ring(b, 300, 250, 5, [], K5),
  'Cyclopentadienone': b => { const a = ring(b, 300, 250, 5, [], K5), p = b.pos(a[0]); b.b(a[0], b.a('O', p.x, p.y - BL), 2); },
  'Tropone': b => { const a = ring(b, 300, 250, 7, [], K7), p = b.pos(a[0]); b.b(a[0], b.a('O', p.x, p.y - BL), 2); },
  'Cyclohexane': b => ring(b, 300, 250, 6, [], [1, 1, 1, 1, 1, 1]),
  'CHFClBr (wedge F)': b => { const c = b.a('C', 300, 250), h = BL * 0.87; b.b(c, b.a('Br', 300 - h, 230)); b.b(c, b.a('Cl', 300 + h, 230)); b.b(c, b.a('F', 300, 250 + BL), 1, 1); },
  'trans-2-butene': b => { const c = [b.a('C', 220, 285), b.a('C', 240, 250), b.a('C', 280, 250), b.a('C', 300, 215)]; b.b(c[0], c[1]); b.b(c[1], c[2], 2); b.b(c[2], c[3]); },
  'cis-2-butene': b => { const c = [b.a('C', 220, 285), b.a('C', 240, 250), b.a('C', 280, 250), b.a('C', 300, 285)]; b.b(c[0], c[1]); b.b(c[1], c[2], 2); b.b(c[2], c[3]); }
};
function cen2(b, ids) { let x = 0, y = 0; ids.forEach(i => { x += b.pos(i).x; y += b.pos(i).y; }); return [x / ids.length, y / ids.length]; }
Object.keys(EX).forEach(k => { const o = document.createElement('option'); o.textContent = k; $('#examples').appendChild(o); });
$('#examples').onchange = e => {
  if (!EX[e.target.value]) return;
  nid = 1; const b = B(); EX[e.target.value](b); mol = b.m; sel = null; e.target.value = ''; commit();
};

/* ---------- presets: added onto the drawing, the drawing is never replaced ---------- */
let curPre = null;
const nOf = (id, lo, hi, d) => Math.max(lo, Math.min(hi, parseInt($(id).value, 10) || d));
const RING = {   // els = atoms around the ring, ord[k] = order of the bond joining atom k and k+1
  cyclo: () => { const n = nOf('#n-cyc', 3, 12, 6); return { els: Array(n).fill('C'), ord: Array(n).fill(1) }; },
  benzene: () => ({ els: Array(6).fill('C'), ord: K6 }),
  pyridine: () => ({ els: ['N', 'C', 'C', 'C', 'C', 'C'], ord: K6 }),
  pyrimidine: () => ({ els: ['N', 'C', 'N', 'C', 'C', 'C'], ord: K6 }),
  pyrrole: () => ({ els: ['N', 'C', 'C', 'C', 'C'], ord: K5 }),
  furan: () => ({ els: ['O', 'C', 'C', 'C', 'C'], ord: K5 }),
  thiophene: () => ({ els: ['S', 'C', 'C', 'C', 'C'], ord: K5 }),
  imidazole: () => ({ els: ['N', 'C', 'N', 'C', 'C'], ord: K5 })
};
const HX = BL * 0.866, VY = BL * 0.5;
const FRAG = {   // local frame: atom 0 is the attachment atom at the origin, the fragment grows along +y
  alkane: b => { const n = nOf('#n-alk', 1, 30, 4); let prev = null; for (let i = 0; i < n; i++) { const id = b.a('C', (i % 2) * VY, i * HX); if (prev !== null) b.b(prev, id, 1); prev = id; } },
  naphthalene: b => { const a = ring(b, 0, 0, 6, [], K6); fuse(b, a, a[1], a[2], 6, [], [1, 1, 2, 1, 2, 1], cen2(b, a)); },
  indole: b => { const a = ring(b, 0, 0, 6, [], K6); fuse(b, a, a[1], a[2], 5, ['N', 'C', 'C'], [1, 1, 1, 2, 1], cen2(b, a)); },
  cooh: b => { const c = b.a('C', 0, 0); b.b(c, b.a('O', -HX, VY), 2); b.b(c, b.a('O', HX, VY), 1); },
  cho: b => { const c = b.a('C', 0, 0); b.b(c, b.a('O', HX, VY), 2); },
  nitro: b => { const n = b.a('N', 0, 0, 1); b.b(n, b.a('O', -HX, VY), 2); b.b(n, b.a('O', HX, VY, -1), 1); },
  cn: b => { const c = b.a('C', 0, 0); b.b(c, b.a('N', 0, BL), 3); },
  vinyl: b => { const c = b.a('C', 0, 0); b.b(c, b.a('C', HX, VY), 2); },
  ethynyl: b => { const c = b.a('C', 0, 0); b.b(c, b.a('C', 0, BL), 3); },
  ipr: b => { const c = b.a('C', 0, 0); b.b(c, b.a('C', -HX, VY), 1); b.b(c, b.a('C', HX, VY), 1); },
  ome: b => { const o = b.a('O', 0, 0); b.b(o, b.a('C', HX, VY), 1); }
};
// fuse a new ring onto an existing bond, choosing the Kekulé pattern that keeps every valence legal
function fuseRing(name, bd) {
  const { g, rg } = G(), ia = g.idx.get(bd.a), ib = g.idx.get(bd.b), A1 = mol.atoms[ia], B1 = mol.atoms[ib];
  const mi = mol.bonds.indexOf(bd), gbi = g.bonds.findIndex(x => x.mi === mi);
  const sides = rg.rings.filter(r => r.bonds.includes(gbi)).sort((p, q) => p.size - q.size);
  if (sides.length >= 2) { $('#status').textContent = 'That bond is already shared by two rings.'; return true; }
  const base = RING[name](), n = base.els.length, maxV = a => Math.max(...C.valences(a.el, a.charge || 0));
  let best = null;
  for (const dir of [1, -1]) for (let s = 0; s < n; s++) {
    const pos = t => (((s + dir * t) % n) + n) % n, bi = t => dir === 1 ? pos(t) : pos(t + 1);
    if (base.els[pos(0)] !== 'C' || base.els[pos(1)] !== 'C') continue;
    const bad = (g.bos[ia] + base.ord[bi(n - 1)] > maxV(A1) ? 1 : 0) + (g.bos[ib] + base.ord[bi(1)] > maxV(B1) ? 1 : 0);
    const score = bad * 10 + (base.ord[bi(0)] === bd.order ? 0 : 1);
    if (!best || score < best.score) best = { score, pos, bi };
  }
  if (!best) { $('#status').textContent = 'This ring cannot be fused onto that bond.'; return true; }
  const els = [], ord = [];
  for (let t = 0; t < n; t++) { if (t >= 2) els.push(base.els[best.pos(t)]); ord.push(base.ord[best.bi(t)]); }
  let away;
  if (sides.length) away = cen(sides[0].atoms);
  else {
    const o = []; [ia, ib].forEach(x => g.nbrs[x].forEach(e => { if (e.j !== ia && e.j !== ib) o.push(mol.atoms[e.j]); }));
    away = o.length ? [o.reduce((s, a) => s + a.x, 0) / o.length, o.reduce((s, a) => s + a.y, 0) / o.length] : [A1.x - (B1.y - A1.y), A1.y + (B1.x - A1.x)];
  }
  const b = B();
  b.m.atoms.push({ id: A1.id, el: 'C', x: A1.x, y: A1.y, charge: 0 }, { id: B1.id, el: 'C', x: B1.x, y: B1.y, charge: 0 });
  fuse(b, [], A1.id, B1.id, n, els, ord, away);
  b.m.atoms = b.m.atoms.filter(a => a.id !== A1.id && a.id !== B1.id);
  b.m.bonds = b.m.bonds.filter(k => !((k.a === A1.id && k.b === B1.id) || (k.a === B1.id && k.b === A1.id)));
  mol.atoms.push(...b.m.atoms); mol.bonds.push(...b.m.bonds);
  commit(); return true;
}
function placePreset(name, pt, anchor, bond) {
  if (!RING[name] && !FRAG[name]) return;
  if (!anchor && bond && RING[name] && fuseRing(name, bond)) return;
  const b = B();
  if (RING[name]) { const r = RING[name](); ring(b, 0, 0, r.els.length, r.els, r.ord); } else FRAG[name](b);
  const A = b.m.atoms, p0 = { x: A[0].x, y: A[0].y };
  let rot = 0, tx, ty;
  const rp = (x, y) => ({ x: x * Math.cos(rot) - y * Math.sin(rot), y: x * Math.sin(rot) + y * Math.cos(rot) });
  if (anchor) {
    const q = freeDir(anchor), d = Math.hypot(q.x - anchor.x, q.y - anchor.y) || 1, ux = (q.x - anchor.x) / d, uy = (q.y - anchor.y) / d;
    rot = Math.atan2(uy, ux) - Math.PI / 2;
    const r0 = rp(p0.x, p0.y); tx = anchor.x + ux * BL - r0.x; ty = anchor.y + uy * BL - r0.y;
  } else {
    let sx = 0, sy = 0; A.forEach(a => { sx += a.x; sy += a.y; }); tx = pt.x - sx / A.length; ty = pt.y - sy / A.length;
  }
  A.forEach(a => { const r = rp(a.x, a.y); a.x = r.x + tx; a.y = r.y + ty; });
  mol.atoms.push(...A); mol.bonds.push(...b.m.bonds);
  if (anchor) mol.bonds.push({ id: nid++, a: anchor.id, b: A[0].id, order: 1, stereo: 0 });
  commit();
}
document.querySelectorAll('.pre').forEach(b => b.onclick = () => {
  document.querySelectorAll('#bondtools button, #elements button, .pre').forEach(x => x.classList.remove('on'));
  b.classList.add('on'); tool = 'preset'; curPre = b.dataset.p;
  $('#status').textContent = 'Placing ' + b.textContent + ' · click canvas / atom' + (RING[curPre] ? ' / bond' : '') + ' · Esc to stop';
});

/* ---------- resizable panels ---------- */
function splitter(handleSel, targetSel, side, min, key) {
  const h = $(handleSel), t = $(targetSel), other = () => $(side === 'l' ? '#panel' : '#tools').getBoundingClientRect().width;
  const maxW = () => Math.max(min, window.innerWidth - 14 - 300 - other());
  const apply = w => { t.style.width = w + 'px'; };
  try { const v = +localStorage.getItem(key); if (v && window.innerWidth > 900) apply(Math.max(min, Math.min(v, maxW()))); } catch (e) {}
  let st = null;
  h.onpointerdown = e => { h.setPointerCapture(e.pointerId); st = { x: e.clientX, w: t.getBoundingClientRect().width }; h.classList.add('drag'); };
  h.onpointermove = e => { if (!st) return; const dx = e.clientX - st.x; apply(Math.max(min, Math.min(maxW(), st.w + (side === 'l' ? dx : -dx)))); size(); size3(); };
  const end = () => { if (!st) return; st = null; h.classList.remove('drag'); try { localStorage.setItem(key, String(Math.round(t.getBoundingClientRect().width))); } catch (e) {} };
  h.onpointerup = end; h.onpointercancel = end;
  h.ondblclick = () => { t.style.width = ''; try { localStorage.removeItem(key); } catch (e) {} size(); size3(); };
}
splitter('#spl', '#tools', 'l', 110, 'ms-tools-w');
splitter('#spr', '#panel', 'r', 240, 'ms-panel-w');

/* ---------- 3D viewer ---------- */
const c3 = $('#v3'), x3 = c3.getContext('2d');
let m3 = null, dip3 = null, yaw = 0.6, pit = -0.35, zoom = 1, dirty = true, d3 = null;
function size3() { const d = window.devicePixelRatio || 1, r = c3.getBoundingClientRect(); if (!r.width || !r.height) return; c3.width = r.width * d; c3.height = r.height * d; x3.setTransform(d, 0, 0, d, 0, 0); dirty = true; }
function draw3() {
  const W = c3.clientWidth, H = c3.clientHeight; x3.clearRect(0, 0, W, H);
  if (!m3) { x3.fillStyle = '#94a3b8'; x3.font = '14px system-ui'; x3.textAlign = 'center'; x3.fillText('Draw a molecule, then press “Convert drawing to 3D”', W / 2, H / 2); return; }
  const cyw = Math.cos(yaw), syw = Math.sin(yaw), cp = Math.cos(pit), sp = Math.sin(pit), st = $('#style').value, showH = $('#showh').checked, lab = $('#lbl').checked;
  let R = 1; m3.atoms.forEach(a => { R = Math.max(R, Math.hypot(a.x, a.y, a.z) + 1); });
  const sc = 0.42 * Math.min(W, H) / R * zoom;
  const pts = m3.atoms.map(a => { const x = a.x * cyw + a.z * syw, z0 = -a.x * syw + a.z * cyw, y = a.y * cp - z0 * sp, z = a.y * sp + z0 * cp; return { X: W / 2 + x * sc, Y: H / 2 - y * sc, z, a, vis: showH || a.el !== 'H' }; });
  const rad = a => st === 'sf' ? C.EL[a.el].vdw * 0.85 : st === 'st' ? 0.13 : C.EL[a.el].vdw * 0.24;
  const items = [];
  pts.forEach((p, i) => { if (p.vis) items.push({ z: p.z, atom: p }); });
  if (st !== 'sf') m3.bonds.forEach(b => { const p = pts[b.a], q = pts[b.b]; if (p.vis && q.vis) items.push({ z: (p.z + q.z) / 2 - 0.01, bond: [p, q, b] }); });
  items.sort((u, v) => u.z - v.z);
  const lw = (st === 'st' ? 0.2 : 0.13) * sc;
  items.forEach(it => {
    if (it.bond) {
      const [p, q, b] = it.bond; let dx = q.X - p.X, dy = q.Y - p.Y; const l = Math.hypot(dx, dy) || 1; dx /= l; dy /= l;
      const offs = b.order === 3 ? [-1.6, 0, 1.6] : b.order === 2 ? [-0.9, 0.9] : [0], w = b.order > 1 ? lw * 0.6 : lw;
      offs.forEach(o => {
        const ox = -dy * o * lw, oy = dx * o * lw, mx = (p.X + q.X) / 2 + ox, my = (p.Y + q.Y) / 2 + oy;
        x3.lineCap = 'round'; x3.strokeStyle = '#1f2937'; x3.lineWidth = w + 2; x3.beginPath(); x3.moveTo(p.X + ox, p.Y + oy); x3.lineTo(q.X + ox, q.Y + oy); x3.stroke();
        [[p, mx, my], [q, mx, my]].forEach(([e, X, Y]) => { x3.strokeStyle = C.EL[e.a.el].c3; x3.lineWidth = w; x3.beginPath(); x3.moveTo(e.X + ox, e.Y + oy); x3.lineTo(X, Y); x3.stroke(); });
      });
    } else {
      const p = it.atom, r = rad(p.a) * sc, g = x3.createRadialGradient(p.X - r * .35, p.Y - r * .35, r * .1, p.X, p.Y, r);
      g.addColorStop(0, '#fff'); g.addColorStop(.25, C.EL[p.a.el].c3); g.addColorStop(1, '#111827');
      x3.fillStyle = g; x3.beginPath(); x3.arc(p.X, p.Y, r, 0, 7); x3.fill();
      if (lab) { x3.fillStyle = '#fff'; x3.strokeStyle = '#000'; x3.lineWidth = 3; x3.font = 'bold 12px system-ui'; x3.textAlign = 'center'; x3.textBaseline = 'middle'; const lt = p.a.el + (p.a.charge ? (Math.abs(p.a.charge) > 1 ? Math.abs(p.a.charge) : '') + (p.a.charge > 0 ? '+' : '−') : ''); x3.strokeText(lt, p.X, p.Y); x3.fillText(lt, p.X, p.Y); x3.textBaseline = 'alphabetic'; }
    }
  });
  if (dip3 && dip3.supported && dip3.mu > 0.05 && $('#showdip').checked) {
    const pr = (x, y, z) => { const X = x * cyw + z * syw, z0 = -x * syw + z * cyw, Y = y * cp - z0 * sp; return [W / 2 + X * sc, H / 2 - Y * sc]; };
    const u = dip3.vec.map(v => v / dip3.mu), Lw = Math.min(Math.max(0.9, 0.5 * dip3.mu), 0.9 * R);
    const t0 = pr(-u[0] * Lw / 2, -u[1] * Lw / 2, -u[2] * Lw / 2), t1 = pr(u[0] * Lw / 2, u[1] * Lw / 2, u[2] * Lw / 2);
    const ang = Math.atan2(t1[1] - t0[1], t1[0] - t0[0]), hd = 15;
    x3.lineCap = 'round'; x3.lineJoin = 'round';
    x3.strokeStyle = '#000'; x3.lineWidth = 9; x3.beginPath(); x3.moveTo(t0[0], t0[1]); x3.lineTo(t1[0] - Math.cos(ang) * hd * 0.6, t1[1] - Math.sin(ang) * hd * 0.6); x3.stroke();
    x3.strokeStyle = '#fbbf24'; x3.lineWidth = 5; x3.beginPath(); x3.moveTo(t0[0], t0[1]); x3.lineTo(t1[0] - Math.cos(ang) * hd * 0.6, t1[1] - Math.sin(ang) * hd * 0.6); x3.stroke();
    x3.beginPath(); x3.moveTo(t1[0], t1[1]); x3.lineTo(t1[0] - hd * Math.cos(ang - 0.45), t1[1] - hd * Math.sin(ang - 0.45)); x3.lineTo(t1[0] - hd * Math.cos(ang + 0.45), t1[1] - hd * Math.sin(ang + 0.45)); x3.closePath();
    x3.fillStyle = '#fbbf24'; x3.strokeStyle = '#000'; x3.lineWidth = 2; x3.fill(); x3.stroke();
    x3.fillStyle = '#fff'; x3.beginPath(); x3.arc(t0[0], t0[1], 7, 0, 7); x3.fill(); x3.strokeStyle = '#000'; x3.lineWidth = 1.5; x3.stroke();
    x3.fillStyle = '#000'; x3.font = 'bold 11px system-ui'; x3.textAlign = 'center'; x3.textBaseline = 'middle'; x3.fillText('+', t0[0], t0[1] + 0.5);
    x3.font = 'bold 13px system-ui'; x3.strokeStyle = '#000'; x3.lineWidth = 3; x3.fillStyle = '#fbbf24';
    const lx = t1[0] + Math.cos(ang) * 8, ly = t1[1] + Math.sin(ang) * 8 - 10, txt = 'μ ≈ ' + dip3.mu.toFixed(2) + ' D';
    x3.strokeText(txt, lx, ly); x3.fillText(txt, lx, ly); x3.textBaseline = 'alphabetic';
  }
}
function loop() { if ($('#spin').checked && m3 && !d3 && !$('#tab-d').hidden) { yaw += 0.01; dirty = true; } if (dirty && !$('#tab-d').hidden) { draw3(); dirty = false; } requestAnimationFrame(loop); }
c3.onpointerdown = e => { c3.setPointerCapture(e.pointerId); d3 = { x: e.clientX, y: e.clientY }; };
c3.onpointermove = e => { if (!d3) return; yaw += (e.clientX - d3.x) * 0.01; pit += (e.clientY - d3.y) * 0.01; d3 = { x: e.clientX, y: e.clientY }; dirty = true; };
c3.onpointerup = () => { d3 = null; };
c3.onwheel = e => { e.preventDefault(); zoom = Math.max(0.3, Math.min(4, zoom * (e.deltaY < 0 ? 1.1 : 0.91))); dirty = true; };
['style', 'showh', 'lbl', 'showdip'].forEach(id => $('#' + id).onchange = () => { dirty = true; });
$('#to3d').onclick = () => {
  if (!mol.atoms.length) { m3 = null; dirty = true; return; }
  const b = $('#to3d'); b.textContent = 'Computing…';
  setTimeout(() => { m3 = C.embed3D(mol); try { dip3 = C.dipole(mol, m3); } catch (er) { dip3 = null; } b.textContent = 'Convert drawing to 3D'; yaw = 0.6; pit = -0.35; dirty = true; }, 20);
};
$('#dip').onclick = () => {
  const out = $('#result');
  if (!mol.atoms.length) { out.innerHTML = '<div class="warn">Draw a molecule first.</div>'; return; }
  const { g } = G(), bad = mol.atoms.map((a, i) => g.ok[i] ? '' : a.el + (i + 1)).filter(Boolean);
  m3 = C.embed3D(mol); dip3 = C.dipole(mol, m3); dirty = true;
  const warn = bad.length ? '<div class="warn">Over-valent atom(s): ' + bad.join(', ') + ' — results may be unreliable.</div>' : '';
  if (!dip3.supported) { out.innerHTML = warn + '<div class="warn">' + dip3.reason + '</div>'; return; }
  const cls = { 'non-polar': 'dpn', 'weakly polar': 'dpw', 'polar': 'dpp', 'strongly polar': 'dps' }[dip3.cls];
  out.innerHTML = warn + '<div class="card"><div class="lbl">Dipole moment (estimate)</div><div class="big">' + (dip3.mu < 0.05 ? '≈ 0' : dip3.mu.toFixed(2)) + ' D</div>' +
    '<span class="chip ' + cls + '">' + dip3.cls + '</span>' +
    '<div class="hint">' + (dip3.mu < 0.1 ? 'The bond dipoles cancel (symmetric shape) or are negligible.' : 'The arrow (δ+ → δ−) is drawn in the 3D view tab; tick “dipole” there.') + '</div></div>' +
    '<div class="hint">Rough estimate: empirical bond moments added up on a quick 3D geometry. Typical error is about 0.3–0.5 D, larger for amides, acids and flexible molecules (conformation matters). Not defined for ions.</div>';
};
$('#expxyz').onclick = () => { if (m3) dl('molecule.xyz', C.toXYZ(m3)); };
document.querySelectorAll('.tab').forEach(t => t.onclick = () => {
  document.querySelectorAll('.tab').forEach(x => x.classList.toggle('on', x === t));
  $('#tab-a').hidden = t.dataset.tab !== 'a'; $('#tab-d').hidden = t.dataset.tab !== 'd';
  size3(); dirty = true;
});

/* ---------- start ---------- */
hist = [snap()]; hi = 0; update(); loop();

/* ---------- Layout Mode Toggle ---------- */
const modeBtn = $('#mode-toggle');
function setMode(isMobile) {
  if (isMobile) {
    document.body.classList.remove('laptop-mode');
    document.body.classList.add('mobile-mode');
    modeBtn.innerHTML = '💻 Switch to Laptop View';
  } else {
    document.body.classList.remove('mobile-mode');
    document.body.classList.add('laptop-mode');
    modeBtn.innerHTML = '📱 Switch to Mobile View';
  }
  
  // Force canvases to recalculate their bounding boxes after layout shift
  setTimeout(() => {
    window.dispatchEvent(new Event('resize')); 
  }, 100);
}

// Ensure the mode button functions correctly when clicked
if (modeBtn) {
  modeBtn.onclick = () => setMode(document.body.classList.contains('laptop-mode'));
}

// Boot up directly into mobile mode to ensure UI builds correctly
setMode(true);
})();
