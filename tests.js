/* tests.js — automated checks for MolValenceX.   Run it in the browser (tests.html) or with Node.
   Browser: open tests.html (same folder as chem.js and smiles.js) - no install needed.
   Node:    node tests.js  (optional, same folder).
   Add a new line to a table below every time you find a bug, so it can never come back silently. */
(function () {
const isNode = typeof require === 'function' && typeof module !== 'undefined';
if (isNode) { require('./chem.js'); require('./smiles.js'); }
const C = isNode ? require('./chem.js') : window.Chem;
if (!C || !C.fromSMILES) { throw new Error('chem.js and smiles.js must be loaded before tests.js'); }
const secs = []; let cur = null;
const out = s => { console.log(s); if (typeof window !== 'undefined' && window.__testOut) window.__testOut(s); };

let pass = 0, fail = 0, known = 0;
function check(name, ok, detail, knownBug) {
  if (!cur) { cur = { name: '(start)', pass: 0, fail: 0, known: 0 }; secs.push(cur); }
  if (ok) { pass++; cur.pass++; return; }
  if (knownBug) { known++; cur.known++; out('  KNOWN  ' + name + (detail ? '  -> ' + detail : '')); return; }
  fail++; cur.fail++; out('  FAIL   ' + name + (detail ? '  -> ' + detail : ''));
}
// SMILES -> drawing -> analysis objects
function load(smi) {
  const m = C.fromSMILES(smi, 40);
  const mol = { atoms: m.atoms.map((a, i) => ({ ...a, id: i })), bonds: m.bonds.map((b, k) => ({ ...b, id: k })) };
  const g = C.buildGraph(mol), rg = C.findRings(g);
  return { m, mol, g, rg };
}
const section = t => { cur = { name: t, pass: 0, fail: 0, known: 0 }; secs.push(cur); out('\n' + t); };

/* 1. formula, molar mass, degree of unsaturation — values from PubChem / textbooks */
section('Formula, mass, DBE');
[
  // name, SMILES, formula, mass (g/mol), DBE
  ['methane', 'C', 'CH4', 16.04, 0],
  ['ethanol', 'CCO', 'C2H6O', 46.07, 0],
  ['benzene', 'c1ccccc1', 'C6H6', 78.11, 4],
  ['aspirin', 'CC(=O)Oc1ccccc1C(=O)O', 'C9H8O4', 180.16, 6],
  ['caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'C8H10N4O2', 194.19, 6],
  ['paracetamol', 'CC(=O)Nc1ccc(O)cc1', 'C8H9NO2', 151.16, 5],
  ['ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'C13H18O2', 206.28, 5],
  ['glucose', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', 'C6H12O6', 180.16, 1],
  ['anthracene', 'c1ccc2cc3ccccc3cc2c1', 'C14H10', 178.23, 10],
  ['chloroform', 'ClC(Cl)Cl', 'CHCl3', 119.38, 0],
].forEach(([n, s, f, mass, dbe]) => {
  try {
    const { mol, g } = load(s), r = C.formula(mol, g);
    check(n + ' formula', r.formula === f, r.formula + ' != ' + f);
    check(n + ' mass', Math.abs(r.mass - mass) < 0.05, r.mass.toFixed(3) + ' != ' + mass);
    check(n + ' DBE', r.dbe === dbe, r.dbe + ' != ' + dbe);
  } catch (e) { check(n, false, 'threw ' + e.message); }
});

/* 2. aromaticity (Hückel) — expected type of every ring */
section('Aromaticity');
[
  ['benzene', 'c1ccccc1', ['aromatic']],
  ['pyridine', 'c1ccncc1', ['aromatic']],
  ['pyrrole', 'c1cc[nH]c1', ['aromatic']],
  ['furan', 'c1ccoc1', ['aromatic']],
  ['thiophene', 'c1ccsc1', ['aromatic']],
  ['imidazole', 'c1c[nH]cn1', ['aromatic']],
  ['naphthalene', 'c1ccc2ccccc2c1', ['aromatic', 'aromatic']],
  ['azulene', 'c1ccc2cccc2cc1', ['aromatic', 'aromatic']],
  ['indole', 'c1ccc2[nH]ccc2c1', ['aromatic', 'aromatic']],
  ['tropylium', 'C1=CC=CC=C[CH+]1', ['aromatic']],
  ['cyclopentadienide', '[CH-]1C=CC=C1', ['aromatic']],
  ['cyclobutadiene', 'C1=CC=C1', ['antiaromatic']],
  ['cyclohexane', 'C1CCCCC1', ['nonaromatic']],
  ['cyclohexene', 'C1=CCCCC1', ['nonaromatic']],
  ['cyclopentadiene', 'C1=CCC=C1', ['nonaromatic']],
].forEach(([n, s, exp]) => {
  try {
    const { mol, g, rg } = load(s), got = C.aromaticity(mol, g, rg).rings.map(r => r.type).sort();
    check(n, JSON.stringify(got) === JSON.stringify(exp.slice().sort()), got.join('/') + ' != ' + exp.join('/'));
  } catch (e) { check(n, false, 'threw ' + e.message); }
});

/* 3. R/S — sorted list of CIP labels (checked against PubChem) */
section('R/S (CIP)');
[
  ['L-alanine', 'C[C@H](N)C(=O)O', ['S']],
  ['D-alanine', 'C[C@@H](N)C(=O)O', ['R']],
  ['L-lactic acid', 'C[C@H](O)C(=O)O', ['S']],
  ['(R)-CHFClBr', 'F[C@](Cl)(Br)[H]', null],             // null = only check it does not crash
  ['D-glucose (pyranose)', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', null],
  ['achiral ethanol', 'CCO', []],
  ['achiral propan-2-ol', 'CC(C)O', []],
].forEach(([n, s, exp]) => {
  try {
    const { mol, g, rg } = load(s), got = C.rsStereo(mol, g, rg).centres.map(c => c.label).sort();
    if (exp) check(n, JSON.stringify(got) === JSON.stringify(exp), got.join(',') + ' != ' + exp.join(','));
    else check(n + ' (no crash)', true);
  } catch (e) { check(n, false, 'threw ' + e.message); }
});

/* 4. SMILES round trip: import -> export -> import must give the same formula, ring count and stereo labels */
section('SMILES round trip');
const RT = ['CCO', 'c1ccccc1', 'c1ccc2ccccc2c1', 'c1ccc2cc3ccccc3cc2c1', 'c1ccc2c(c1)ccc1ccccc12', 'CCOc1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O',
  'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'C[C@H](N)C(=O)O', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', 'C1CC2CCC1C2', 'C1CCC2(CC1)CCCC2',
  'c1ccc(cc1)-c1ccccc1', 'C[N+](C)(C)C', '[O-]C(=O)C', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'C#CC=C', 'FC(F)(F)c1ccccc1', 'C1CC1', 'C1CCCCCCC1'];
RT.forEach(s => {
  try {
    const a = load(s), out = C.toSMILES(a.mol), b = load(out);
    const fa = C.formula(a.mol, a.g).formula, fb = C.formula(b.mol, b.g).formula;
    const ra = a.rg.rings.length, rb = b.rg.rings.length;
    const sa = C.rsStereo(a.mol, a.g, a.rg).centres.map(c => c.label).sort().join(''), sb = C.rsStereo(b.mol, b.g, b.rg).centres.map(c => c.label).sort().join('');
    check(s, fa === fb && ra === rb && sa === sb, s + ' -> ' + out + '  (' + fa + '/' + fb + ', rings ' + ra + '/' + rb + ', stereo ' + sa + '/' + sb + ')');
  } catch (e) { check(s, false, 'threw ' + e.message); }
});

/* 5. PubChem-style input: PubChem sends upper-case Kekulé SMILES — must equal the aromatic form */
section('Kekulé == aromatic');
[['CCOC1=CC=CC=C1', 'CCOc1ccccc1'], ['C1=CC=CC=C1', 'c1ccccc1'], ['C1=CC=C2C=CC=CC2=C1', 'c1ccc2ccccc2c1'],
 ['C1=CC=C2C(=C1)C=CC3=CC=CC=C32', 'c1ccc2c(c1)ccc1ccccc12'], ['C1=CNC=C1', 'c1cc[nH]c1']].forEach(([k, a]) => {
  try {
    const A = load(k), B = load(a);
    check(k, C.formula(A.mol, A.g).formula === C.formula(B.mol, B.g).formula && A.rg.rings.length === B.rg.rings.length &&
      C.aromaticity(A.mol, A.g, A.rg).rings.map(r => r.type).join() === C.aromaticity(B.mol, B.g, B.rg).rings.map(r => r.type).join());
  } catch (e) { check(k, false, 'threw ' + e.message); }
});

/* 6. 2D layout quality: every bond about one bond-length long, no two atoms nearly on top of each other */
section('2D layout geometry');
['c1ccccc1', 'c1ccc2ccccc2c1', 'c1ccc2cc3ccccc3cc2c1', 'c1ccc2c(c1)ccc1ccccc12', 'CCOc1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', 'c1ccc(cc1)-c1ccccc1',
 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O', 'c1cc2ccc3cccc4ccc(c1)c2c34', 'C1CCCCC1', 'CCCCCCCCCC'].forEach(s => {
  try {
    const { m } = load(s);
    const bl = m.bonds.map(b => Math.hypot(m.atoms[b.a].x - m.atoms[b.b].x, m.atoms[b.a].y - m.atoms[b.b].y));
    let md = 1e9; for (let i = 0; i < m.atoms.length; i++) for (let j = i + 1; j < m.atoms.length; j++) md = Math.min(md, Math.hypot(m.atoms[i].x - m.atoms[j].x, m.atoms[i].y - m.atoms[j].y));
    check(s, Math.min(...bl) > 34 && Math.max(...bl) < 46 && md > 28, 'bond ' + Math.min(...bl).toFixed(0) + '-' + Math.max(...bl).toFixed(0) + ', closest atoms ' + md.toFixed(0));
  } catch (e) { check(s, false, 'threw ' + e.message); }
});
// bridged cages: known to be a little uneven before "Clean up"
check('norbornane layout', (() => { const { m } = load('C1CC2CCC1C2'); let md = 1e9; for (let i = 0; i < m.atoms.length; i++) for (let j = i + 1; j < m.atoms.length; j++) md = Math.min(md, Math.hypot(m.atoms[i].x - m.atoms[j].x, m.atoms[i].y - m.atoms[j].y)); return md > 36; })(), 'atoms closer than 36 px', true);

/* 7. bad input must be rejected with an error, never crash or draw garbage */
section('Invalid SMILES rejected');
['', 'C(', 'C)', 'C1CC', 'c1ccc', '[Xx]', 'C[', 'CC==C', '123', 'C%'].forEach(s => {
  let threw = false; try { C.fromSMILES(s, 40); } catch (e) { threw = true; }
  check('"' + s + '"', threw, 'was accepted');
});

/* 8. E/Z from SMILES / and \ marks */
section('E/Z from SMILES');
[['F/C=C/F', 'E'], ['F/C=C\\F', 'Z'], ['C/C=C/C', 'E'], ['C/C=C\\C', 'Z'], ['C(/F)=C/F', 'Z'], ['C(\\F)=C/F', 'E'], ['F/C=C/C=C/F', 'E,E'], ['F/C=C\\C=C/F', 'Z,Z'], ['C/C=C/C(=O)O', 'E'], ['OC(=O)/C=C\\C(=O)O', 'Z'], ['OC(=O)/C=C/C(=O)O', 'E'], ['Cl/C(F)=C(/Br)I', 'Z']].forEach(([s, exp]) => {
  try { const { mol, g, rg } = load(s), got = C.ezStereo(mol, g, rg).bonds.map(b => b.label).join();
    check(s + ' -> ' + exp, got === exp, 'got "' + got + '"'); } catch (e) { check(s, false, 'threw ' + e.message); }
});


/* 9. importer: more valid SMILES with the formula they must give (checked by hand / PubChem) */
section('Importer: formulas');
[
  ['biphenyl, implicit bond', 'c1ccccc1c1ccccc1', 'C12H10'], ['biphenyl, explicit', 'c1ccccc1-c1ccccc1', 'C12H10'],
  ['explicit [H]', '[H]C([H])([H])[H]', 'CH4'], ['explicit H chiral', '[H][C@](F)(Cl)Br', 'CHBrClF'], ['deuterium', '[2H]C([2H])([2H])[2H]', 'CH4'],
  ['ammonium', '[NH4+]', 'H4N⁺'], ['pyridinium', 'C[n+]1ccccc1', 'C6H8N⁺'], ['pyridine N-oxide', '[O-][n+]1ccccc1', 'C5H5NO'],
  ['nitromethane', 'C[N+](=O)[O-]', 'CH3NO2'], ['nitro, pentavalent', 'CN(=O)=O', 'CH3NO2'], ['DMSO', 'CS(C)=O', 'C2H6OS'], ['sulfate', 'OS(=O)(=O)O', 'H2O4S'],
  ['sodium chloride', '[Na+].[Cl-]', null], ['acetate', 'CC(=O)[O-]', 'C2H3O2⁻'], ['pyrrole', 'c1cc[nH]c1', 'C4H5N'], ['indole', 'c1ccc2[nH]ccc2c1', 'C8H7N'],
  ['2-pyridone', 'O=c1cccc[nH]1', 'C5H5NO'], ['uracil', 'O=c1cc[nH]c(=O)[nH]1', 'C4H4N2O2'], ['purine', 'c1ncc2[nH]cnc2n1', 'C5H4N4'],
  ['imidazole 1-Me', 'Cn1ccnc1', 'C4H6N2'], ['thiophene', 'c1ccsc1', 'C4H4S'], ['furan', 'c1ccoc1', 'C4H4O'], ['pyrene', 'c1cc2ccc3cccc4ccc(c1)c2c34', 'C16H10'],
  ['ring label %10', 'C%10CCCCC%10', 'C6H12'], ['label reuse', 'C1CC1C1CC1', 'C6H10'], ['dot, two molecules', 'CCO.CC(=O)O', null], ['trailing name', 'CCO ethanol', 'C2H6O'],
  ['triple bond', 'C#CC#N', 'C3HN'], ['cumulene', 'C=C=C', 'C3H4'], ['quinone', 'O=C1C=CC(=O)C=C1', 'C6H4O2'], ['azulene', 'c1ccc2cccc2cc1', 'C10H8'],
  ['paracetamol', 'CC(=O)Nc1ccc(O)cc1', 'C8H9NO2'], ['penicillin G', 'CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O', 'C16H18N2O4S'],
  ['morphine', 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O', 'C17H19NO3'],
  ['cholesterol', 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C', 'C27H46O'],
  ['ring closure bond =', 'C=1CCCCC=1', 'C6H10'], ['bracket B', '[BH3]', 'BH3'], ['silane', '[SiH4]', 'H4Si'],
].forEach(([n, s, f]) => {
  try { const { mol, g } = load(s), got = C.formula(mol, g).formula; if (f) check(n, got === f, s + ' -> ' + got + ' != ' + f); else check(n + ' (no crash)', true); }
  catch (e) { check(n, false, 'threw ' + e.message); }
});
// the three stereocentres of these must survive import
section('Importer: stereo kept');
[['penicillin G', 'CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O', 'RRS'], ['L-alanine', 'C[C@H](N)C(=O)O', 'S'], ['(S)-CHBrClF explicit H', '[H][C@](F)(Cl)Br', null]].forEach(([n, s, exp]) => {
  try { const { mol, g, rg } = load(s), got = C.rsStereo(mol, g, rg).centres.map(c => c.label).sort().join('');
    if (exp) check(n, got === exp.split('').sort().join(''), got + ' != ' + exp); else check(n + ' has 1 centre', got.length === 1, got); } catch (e) { check(n, false, 'threw ' + e.message); }
});
section('Importer: must be rejected');
['CC==C', 'C(C)(C)(C)(C)C', 'F(C)C', '[CH2]', '[Fe]', '*', 'n1cccc1', 'C1CC1CC1', 'C1C1', 'C1CC', 'CC(', '[NH5]', 'C=', '=C', 'C..C', 'c', '[C@@H', 'CC#', '1C', 'C%1', 'C-'].forEach(s => {
  let threw = false; try { C.fromSMILES(s, 40); } catch (e) { threw = true; }
  check('"' + s + '"', threw, 'was accepted');
});

out('\n' + pass + ' passed, ' + fail + ' failed, ' + known + ' known issues');
if (typeof window !== 'undefined' && window.__testDone) window.__testDone(secs, pass, fail, known);
if (isNode) process.exit(fail ? 1 : 0);
})();
