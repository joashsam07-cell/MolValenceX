<div align="center">

# ⬡ MolValenceX

**A free, instant, no-install molecule drawer and chemistry analyzer that runs entirely in your browser.**

Draw a structure → get its formula, aromaticity, stereochemistry, dipole, solubility parameters and a 3D model. No sign-up, no server, no tracking.

[![Live demo](https://img.shields.io/badge/Live%20demo-open%20app-6366f1?style=for-the-badge)](https://joashsam07-cell.github.io/MolValenceX/)
[![License: MIT](https://img.shields.io/badge/License-MIT-16a34a?style=for-the-badge)](LICENSE)
![Dependencies](https://img.shields.io/badge/dependencies-0-0ea5e9?style=for-the-badge)
![Made for](https://img.shields.io/badge/made%20for-students%20%26%20educators-d97706?style=for-the-badge)

**[▶ Launch MolValenceX](https://joashsam07-cell.github.io/MolValenceX/)**

<!-- Add a screenshot or GIF here:  ![MolValenceX screenshot](docs/screenshot.png) -->

</div>

---

## Why MolValenceX?

Most free viewers only *draw*. MolValenceX also *explains*: every result comes with the reasoning behind it (for example, which ring has how many π electrons and why that makes it aromatic), so it works as a study aid, not just a calculator.

Built by a B.Pharm student for organic chemistry, medicinal chemistry and drug-discovery learning.

## Features

### ✏️ Drawing
- Single, double and triple bonds, plus **wedge / dash** stereo bonds
- Elements C, H, N, O, S, P, F, Cl, Br, I, B, Si, Se, plus Li, Na, K, Mg, Ca, Al
- Formal charges (+ / −), eraser, select / move / pan
- One-click **presets**: n-alkanes, cycloalkanes (3–12), benzene, naphthalene, pyridine, pyrimidine, pyrrole, furan, thiophene, imidazole, indole
- **Functional groups**: COOH, CHO, NO₂, C≡N, vinyl, ethynyl, i-Pr, OMe (fuse rings onto existing bonds, attach groups to atoms)
- Undo / redo, **Clean up** (2D relaxation), ready-made **example molecules**
- Resizable panels, **laptop and mobile layouts**
- Built-in **guided tutorial** of every feature

### 🧪 Analysis
| Tool | What you get |
|---|---|
| **Formula & mass** | Molecular formula, molar mass, degree of unsaturation, heavy-atom / H count, live as you draw |
| **Aromaticity** | Ring-by-ring **Hückel (4n+2)** verdict: aromatic, antiaromatic or non-aromatic, with π-electron counts and a written reason. Handles fused rings and heterocycles (pyrrole, furan, azulene, tropylium, cyclopentadienide…) |
| **Dipole moment** | Estimated μ in Debye with a polarity class, drawn as an arrow in the 3D view |
| **E/Z, cis / trans** | CIP-priority E/Z assignment for double bonds |
| **R/S & optical activity** | CIP R/S for every stereocentre, from wedge/dash or Fischer drawings; chirality and meso detection |
| **Solubility (Hansen)** | δD, δP, δH estimates (Stefanis–Panayiotou, with a Hoftyzer–Van Krevelen cross-check) and a ranked list of good / partial / poor solvents |

### 🧊 3D viewer
- One-click **2D → 3D** conversion
- Ball-and-stick, sticks and space-filling styles
- Toggle hydrogens, atom labels, auto-spin and the dipole arrow
- Drag to rotate, scroll to zoom

### 🔁 Import & export
- **Import SMILES**, including stereo (`@` / `@@`), E/Z (`/` `\`), charges, isotopes, explicit `[H]`, aromatic and Kekulé forms
- **Export SMILES** (with tetrahedral stereo), **Export MOL**, **Export XYZ**
- **Search PubChem by name** (aspirin, caffeine, glucose…): only the name you type is sent to PubChem; nothing you draw ever leaves your browser

### ⌨️ Shortcuts
| Key | Action |
|---|---|
| `C` `N` `O` `S` `P` `F` `H` `B` `I` | Choose element (`L` = Cl, `R` = Br) |
| `1` `2` `3` | Single / double / triple bond |
| `Ctrl+Z` / `Ctrl+Y` | Undo / redo |
| `Del` | Erase selection |
| `Esc` | Cancel the current preset |

## Quick start

**Use it online:** open **<https://joashsam07-cell.github.io/MolValenceX/>**

**Run it locally:** there is nothing to build.

```bash
git clone https://github.com/joashsam07-cell/MolValenceX.git
cd MolValenceX
# open index.html in any modern browser
```

## Project structure

```
MolValenceX/
├── index.html     UI layout and SEO metadata
├── style.css      Glassmorphism theme, laptop and mobile modes
├── app.js         Canvas editor, tools, 3D viewer, PubChem search, tutorial
├── chem.js        Chemistry engine (graph, rings, aromaticity, CIP, dipole, Hansen, 2D/3D)
├── smiles.js      SMILES import / export and 2D layout
├── tests.js       Automated tests for chem.js and smiles.js
├── robots.txt
└── sitemap.xml
```

`chem.js` has no DOM code, so the same engine runs in the browser and in Node.js.

## Testing

The chemistry and SMILES code is covered by an automated test suite (formulas, aromaticity, R/S, E/Z, SMILES round trips, 2D layout quality, and rejection of invalid input).

```bash
node tests.js
```

Requires [Node.js](https://nodejs.org) (LTS). Every bug fix should come with a new test line, so it can never silently come back.

## Accuracy and limitations

MolValenceX is **an educational tool, not an analytical reference**.

- Formula, mass, degree of unsaturation, Hückel aromaticity and CIP R/S/E/Z are rule-based and tested against known molecules.
- **Dipole moment, Hansen solubility parameters and 3D geometry are estimates.** Use them to compare and learn, not to report values.
- Hansen verdicts follow "like dissolves like"; they ignore crystal lattice energy, ionisation, temperature and concentration, and predict water poorly.
- Axial and planar chirality are not assessed. Radicals and unusual charge states are not supported.
- Large bridged polycyclic structures may need **Clean up** after import.

## Privacy

Everything runs locally in your browser. The only network request is the optional PubChem name search, and it sends **only the text you typed** in the search box. The app stores just your panel widths and a "tutorial seen" flag in local storage.

## Roadmap

- [ ] Drug-likeness: Lipinski rule of five, logP, TPSA, H-bond donors / acceptors
- [ ] E/Z written out in exported SMILES
- [ ] Support for 2+ charged metal ions
- [ ] PNG / SVG export of the drawing
- [ ] More example molecules and a drug-structure library

## Contributing

Bug reports and ideas are welcome. Open an issue with the molecule (a SMILES string is perfect) and what you expected to see. If you send a pull request, please run `node tests.js` first and add a test for what you changed.

## Author

**Joash Sam** · B.Pharm student · built with Claude and prompt engineering

## License

Released under the [MIT License](LICENSE). © 2026 Joash Sam.

> **Disclaimer:** For educational purposes only and not a reference for analytics. Predicted values may differ from reality.
