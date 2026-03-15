# pcp-thermodynamics

[![CI](https://github.com/enrik03/pcp_thermo/actions/workflows/ci.yml/badge.svg)](https://github.com/enrik03/pcp_thermo/actions/workflows/ci.yml)

Numerical thermodynamic analysis and publication-quality figures for the
Post-Carroll-Pauli (PCP) sector, including pressure scaling, background-field
corrections, and PCP–UR comparisons.

---

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Generated Figures](#generated-figures)
- [Numerical Results](#numerical-results)
- [Running Tests](#running-tests)
- [Project Structure](#project-structure)
- [License](#license)

---

## Overview

This package reproduces and documents the thermodynamic pressure analysis of
the post-Carroll-Pauli (PCP) sector.  The main deliverable is a self-contained
Python script that:

- Computes the key PCP coefficients **C₄**, **C₇**, **T_c^(Φ)**, **T_x**
  from one-loop thermal integrals.
- Generates three publication-quality figures (PNG + PDF).
- Supports command-line configuration of output directory, DPI, and formats.

---

## Installation

### Using pip (recommended)

```bash
pip install -r requirements.txt
```

### Using conda / mamba

```bash
conda env create -f environment.yml
conda activate pcp-thermo
```

### Requirements

- Python ≥ 3.9
- numpy, scipy, matplotlib, pytest

---

## Quick Start

```bash
# Generate all figures (PNG + PDF) into the figures/ directory
python pcp_thermo_figures_corrected.py

# Custom output directory and resolution
python pcp_thermo_figures_corrected.py --outdir out/ --dpi 300

# PNG only
python pcp_thermo_figures_corrected.py --format png

# Use the Makefile
make figures
```

---

## Generated Figures

| File stem | Description |
|-----------|-------------|
| `pcp_pressure_comparison_main` | Normalised PCP vs UR pressure P(T)/T⁴ on a log scale (main text) |
| `pcp_pressure_crossover_appendix` | Ratio P_PCP / P_UR showing the crossover at T_x (appendix) |
| `pcp_pressure_decomposition_appendix` | Decomposition of P_PCP into C₄·T⁴ and C₇·T⁷/M³ terms (appendix) |

---

## Numerical Results

Key values for M = 1 (dimensionless units):

| Quantity | Value |
|----------|-------|
| C_UR (Stefan–Boltzmann) | ≈ 4.058 |
| C₄ | ≈ 2.906 |
| C₇ | ≈ 0.0399 |
| T_c^(Φ) | ≈ 6.93 M |
| T_x (crossover) | ≈ 3.07 M |

See `pcp_thermo_package_notes.md` for full derivations and LaTeX fragments.

---

## Running Tests

```bash
# Run all unit tests
pytest tests/

# Or via Makefile
make test
```

The tests verify:
- Correct sign and order-of-magnitude of C₄ and C₇.
- T_c^(Φ) within expected range.
- All three figure files are created for each requested format.

---

## Project Structure

```
pcp_thermo/
├── pcp_thermo_figures_corrected.py   # Main analysis script
├── pcp_thermo_package_notes.md       # Technical notes, captions, LaTeX
├── README.md
├── requirements.txt
├── environment.yml
├── LICENSE
├── .gitignore
├── tests/
│   └── test_coefficients.py          # Unit tests
├── tasks/
│   └── Makefile                      # make figures / make test
└── .github/
    └── workflows/
        └── ci.yml                    # GitHub Actions CI
```

---

## License

MIT – see [LICENSE](LICENSE).

