# Installation

MolPy is designed to work across multiple platforms and Python versions. This guide covers all installation methods and requirements.

## Requirements

- **Python**: 3.10 or higher
- **Operating Systems**: Windows, macOS, Linux
- **Architecture**: x86_64, ARM64 (Apple Silicon)

## Quick Install

The simplest way to install MolPy:

```bash
pip install molcrafts-molpy
```

## Installation with Extras

MolPy provides optional dependency groups for different use cases:

### Development Dependencies
```bash
pip install molcrafts-molpy[dev]
```
Includes: pytest, black, isort, pytest-cov

### Documentation Dependencies
```bash
pip install molcrafts-molpy[doc]
```
Includes: mkdocs, mkdocs-material, mkdocstrings

### All Extras
```bash
pip install molcrafts-molpy[all]
```

## Core Dependencies

MolPy automatically installs these required packages:

- **numpy** (≥2.0) - Numerical computing
- **molcrafts-molq** - Core quantum chemistry backend
- **python-igraph** (≥0.9.0) - Graph analysis
- **lark** (≥1.1.0) - SMILES/SMARTS parsing
- **pint** - Unit handling
- **freud-analysis** - Analysis algorithms

### Getting Help

- **GitHub Issues**: [Report installation problems](https://github.com/MolCrafts/molpy/issues)
- **Discussions**: [Community support](https://github.com/MolCrafts/molpy/discussions)
- **Documentation**: Check the [FAQ](faq.md) for common issues

## Next Steps

Once installed, proceed to:
- **[Quickstart Guide](quickstart.md)** - Your first MolPy simulation
- **[Quickstart](quickstart.md)** - Getting started with MolPy
- **[Tutorials](../tutorials/block-frame.md)** - Learning core concepts
