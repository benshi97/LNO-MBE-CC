# LNO-MBE-CCSD(T) - Efficient Modelling of Molecular Crystals

`LNO-MBE-CCSD(T)` is a computational framework for performing efficient and accurate predictions of molecular crystals:

- `LNO-MBE-CCSD(T)` is highly flexible, combining the many-body expansion with [MRCC](https://mrcc.hu/), together with periodic solid-state calculations in [VASP](https://www.vasp.at/)

- `LNO-MBE-CCSD(T)` is powered by [QuAcc](https://github.com/Quantum-Accelerators/quacc) that can be efficiently dispatched locally or on HPCs.

## Installation

`LNO-MBE-CCSD(T)` requires python >= 3.11. We recommend starting in a new python environment using [miniconda](https://docs.anaconda.com/miniconda/):

```
conda create --name lnombecc python=3.11
conda activate lnombecc
```

For local development of the code:

1. Clone the repository

```
git clone https://github.com/benshi97/LNO-MBE-CC.git
```

2. Then install the package in editable mode

```
pip install -e .
```

where this command is run in the root directory. All dependences (i.e., QuAcc) will be automatically installed. By using the `-e` flag, the package will be installed in editable mode, meaning that changes to the code will be reflected in the installed package. Installation should only take a few minutes.

Note: You also will need to have [pMBE](https://github.com/kmherman/pMBE) installed to run the code. This can be installed after `LNO-MBE-CCSD(T)` is installed (another simple `pip install`).

## Requirements

### Python packages
- `numpy`
- `ase`
- `quacc`
- MRCC support via `quacc.calculators.mrcc`
- VASP support via `ase.calculators.vasp` (only needed if running periodic HF)

### External executables
- `pmbe` (pMBE command-line tool) — required for `create_fragments`
- `mrcc` executable — required for fragment calculations
- `vasp_[gam,std]` executable — required for periodic HF calculations (optional)

## Environment setup

This workflow relies on `quacc` to locate and run MRCC and VASP.  
Before running any calculations, configure your environment variables.

```bash
# --- MRCC ---
export QUACC_MRCC_CMD="/path/to/mrcc/bin/dmrcc"
export PATH="/path/to/mrcc/bin:$PATH"

# --- VASP ---
export QUACC_VASP_PARALLEL_CMD="mpirun"        # or srun, etc.
export QUACC_VASP_PP_PATH="/path/to/vasp/pseudopotentials"
```
---

## How to run
