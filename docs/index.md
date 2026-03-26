---
template: home.html
hide:
  - navigation
  - toc
---

# MolPy

<p class="lead" markdown>Composable molecular modeling in Python.</p>

<div class="badges" markdown>
  [![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)
  [![License: BSD](https://img.shields.io/badge/license-BSD-green.svg)](https://github.com/MolCrafts/molpy/blob/master/LICENSE)
  [![Docs](https://img.shields.io/badge/docs-mkdocs-blue.svg)](https://molcrafts.github.io/molpy)
  [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
  [![Type checked: ty](https://img.shields.io/badge/type%20checked-ty-blue.svg)](https://github.com/astral-sh/ty)
</div>

<div class="button-group" markdown>
  [Get Started](getting-started/index.md){ .md-button .md-button--primary }
  [Guides](user-guide/index.md){ .md-button }
  [API Reference](api/index.md){ .md-button }
</div>

---

## See it work

<div class="hero-code" markdown>

```python
import molpy as mp

# editable chemistry
mol = mp.parser.parse_molecule("CCO")
mol = mp.tool.generate_3d(mol)
mol.get_topo(gen_angle=True, gen_dihe=True)

# inspectable force field (OPLS-AA ships with MolPy)
ff = mp.io.read_xml_forcefield("oplsaa.xml")
typed = mp.typifier.OplsTypifier(ff).typify(mol)

# portable export
mp.io.write_lammps_system("output/", typed.to_frame(), ff)
```

</div>

---

## Why MolPy

<div class="feature-rows" markdown>

<a class="feature-row" href="tutorials/index.md">
<span class="feature-icon">:material-swap-horizontal:</span>
<span class="feature-body">
<strong>Explicit representations</strong>
Chemistry editing, system snapshots, and force field data are separate objects with visible transitions.
</span>
<span class="feature-arrow">:octicons-arrow-right-24:</span>
</a>

<a class="feature-row" href="tutorials/04_force_field.md">
<span class="feature-icon">:material-database-search:</span>
<span class="feature-body">
<strong>Force fields as data</strong>
Parameters stay inspectable until you convert them to potentials. Catch errors before export.
</span>
<span class="feature-arrow">:octicons-arrow-right-24:</span>
</a>

<a class="feature-row" href="api/index.md">
<span class="feature-icon">:material-puzzle:</span>
<span class="feature-body">
<strong>Composable modules</strong>
Parser, builder, typifier, packer, and I/O are independent. Use one alone or wire several together.
</span>
<span class="feature-arrow">:octicons-arrow-right-24:</span>
</a>

</div>

---

## Learn MolPy

<div class="feature-rows" markdown>

<a class="feature-row" href="getting-started/index.md">
<span class="feature-icon">:material-rocket-launch:</span>
<span class="feature-body">
<strong>Getting Started</strong>
Install, run one example, learn the core objects.
</span>
<span class="feature-arrow">:octicons-arrow-right-24:</span>
</a>

<a class="feature-row" href="tutorials/index.md">
<span class="feature-icon">:material-book-open-variant:</span>
<span class="feature-body">
<strong>Concepts</strong>
Atomistic, Block, Frame, Box, Trajectory, ForceField.
</span>
<span class="feature-arrow">:octicons-arrow-right-24:</span>
</a>

<a class="feature-row" href="user-guide/index.md">
<span class="feature-icon">:material-hammer-wrench:</span>
<span class="feature-body">
<strong>Guides</strong>
Parse, build polymers, typify, export to LAMMPS.
</span>
<span class="feature-arrow">:octicons-arrow-right-24:</span>
</a>

<a class="feature-row" href="developer/index.md">
<span class="feature-icon">:material-code-braces:</span>
<span class="feature-body">
<strong>Developer Guide</strong>
Extend the tool layer or add new core types.
</span>
<span class="feature-arrow">:octicons-arrow-right-24:</span>
</a>

</div>
