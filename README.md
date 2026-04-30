# LNO-MBE-CCSD(T) - Efficient Modeling of Molecular Crystals

`LNO-MBE-CCSD(T)` is a computational framework for performing efficient and accurate predictions of molecular crystals:

- `LNO-MBE-CCSD(T)` is highly flexible, combining the many-body expansion with [MRCC](https://mrcc.hu/), together with periodic solid-state calculations in [VASP](https://www.vasp.at/)

- `LNO-MBE-CCSD(T)` is powered by [QuAcc](https://github.com/Quantum-Accelerators/quacc), allowing it to be efficiently dispatched locally or on HPCs.

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

where this command is run in the root directory. All dependences (i.e., QuAcc) will be automatically installed. By using the `-e` flag, the package will be installed in editable mode, meaning that changes to the code will be reflected in the installed package. The installation time should only be a few minutes.

Note: You also will need to have [pMBE](https://github.com/kmherman/pMBE) installed to run the code. This can be installed after `LNO-MBE-CCSD(T)` is installed (another simple `pip install`).

## Requirements

`LNO-MBE-CCSD(T)` requires only a standard computer with e.g., 16 GB of RAM. The CCSD(T) calculations in MRCC can require more RAM depending on the size of system studied.

This package is supported for both macOS and Linux. The package has been tested on the following systems:
+ Apple: macOS 15 Sequoia
+ Linux: Ubuntu 22.04.5

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

To perform CCSD(T) calculations, [MRCC](https://mrcc.hu/) (>= 2025) must also be installed; free for academics. Calculations with the periodic Hartree-Fock currently supports [VASP](https://vasp.at/), which requires the purchase of a license.

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
export VASP_PP_PATH="/path/to/vasp/pseudopotentials"
```
---


## Quickstart: Run the full workflow

The package is designed to be modular, with each stage of the workflow callable independently.  
A typical end-to-end calculation proceeds through the following steps:

1. **Generate fragments (pMBE)**  
   Use `create_fragments` to decompose the periodic structure into 1-body, 2-body, and 3-body fragments and build the corresponding ASE databases.

2. **Set up LNO-MBE-CCSD(T) inputs**  
   Use `setup_lnombecc_inputs` to:
   - Attach MRCC calculators to each fragment,
   - Write MRCC input files (MINP),
   - Save calculator objects for later execution.

3. **Run fragment calculations (MRCC via QuAcc)**  
   Use `run_lnombecc` to execute all 1B/2B/3B fragment calculations.  
   *Optional:* You may instead run the generated input files manually using your own job submission workflow.

4. **Run periodic HF calculations (optional)**  
   Use `run_periodic_hf` to compute the periodic Hartree–Fock correction using VASP via ASE/QuAcc.

5. **Analyze results**  
   - `analyze_lnombecc_outputs` computes the 1B, 2B, and 3B contributions to the lattice energy.
   - `analyze_periodic_hf_outputs` extracts the periodic HF correction.
   - Combine these terms to obtain the final lattice energy.

Because each stage is independent, you can:
- Generate inputs on one machine and run calculations on another,
- Skip periodic HF entirely,
- Replace the execution step with your own job scheduler or workflow system.


### Example script

Create a file `run_workflow.py`:

```python
from pathlib import Path

from lnombecc.lnombecc import (
    create_fragments,
    setup_lnombecc_inputs,
    run_lnombecc,
    run_periodic_hf,
    analyze_lnombecc_outputs,
    analyze_periodic_hf_outputs,
)

# Inputs
poscar = Path("POSCAR")          # periodic crystal structure. Must be VASP POSCAR format
gas = Path("geometry.in")        # gas-phase geometry (optional). Any ASE format accepted

# Output / run directory
run_dir = Path("LNOMBECC_calcs")

# 1) Create fragments (pMBE) -> system_1b.db, system_2b.db, system_3b.db
create_fragments(
    poscar_filepath=poscar,
    gas_filepath=gas,            # or omit to use the first fragment as gas
    twob_cutoff=12.0,            # cutoff distance for the 2B terms
    threeb_cutoff=300.0,         # cutoff distance for the 3B terms
    full_uc=False,               # use if Z'=1 asymmetric unit cell
    include_periodic_hf=True,    # also create system_periodic.db
)

# 2) Setup inputs + write MINP files + save calculators_*.npy
setup_lnombecc_inputs(
    one_body_db_filepath="system_1b.db",
    two_body_db_filepath="system_2b.db",
    three_body_db_filepath="system_3b.db",
    write_inputs_dir=run_dir,
    calculation_type="lattice",  # or "relative" which uses smaller basis sets and looser LNO thresholds
    memory="46GB",               # memory for the MRCC calculations
    include_periodic_hf=True,    # writes periodic_HF inputs + saves calculators_periodic.npy
)

# 3) Run fragment calculations (MRCC)
run_lnombecc(
    run_directory=run_dir,
    calculators_1b_filepath="calculators_1b.npy",
    calculators_2b_filepath="calculators_2b.npy",
    calculators_3b_filepath="calculators_3b.npy",
)

# 4) Run periodic HF (VASP), optional
run_periodic_hf(
    run_directory=run_dir,
    calculators_periodic_filepath="calculators_periodic.npy",
)

# 5) Analyze fragment outputs
frag_contribs = analyze_lnombecc_outputs(
    run_directory=run_dir,
    calculators_1b_filepath="calculators_1b.npy",
    calculators_2b_filepath="calculators_2b.npy",
    calculators_3b_filepath="calculators_3b.npy",
    calculation_type="lattice",  # must match setup choice
)

# 6) Analyze periodic HF outputs (optional)
hf_elatt = analyze_periodic_hf_outputs(
    run_directory=run_dir,
    calculators_periodic_filepath="calculators_periodic.npy",
)

print("Fragment contributions (eV):")
print(frag_contribs)

print("Periodic HF correction:")
print(hf_elatt)

print("Total lattice energy (fragments + HF correction):")
print(frag_contribs["1B"] + frag_contribs["2B"] + frag_contribs["3B"] + hf_elatt)
```

## Demo

We provide a demonstration for running LNO-MBE-CCSD(T) in [`example/run_workflow.py`](example/run_workflow.py). It features pre-calculated CCSD(T) and (periodic) HF data for the ammonia crystal, found in `example/LNOMBECC_calcs`. As this workflow is restarting from completed calculations, it does not perform any quantum chemistry calculations. If the user would like to perform these calculations, please move or delete the example/LNOMBECC_calcs folder and follow the guidance within the demo to initialise these calculations. The expected output from this demonstration should be:
```
Fragment contributions (eV):
{'1B': -0.010775053988799854, '2B': -0.29135135812023716, '3B': 0.014853126064387823}

Periodic HF correction:
-0.093841577500001

Total lattice energy (fragments + HF correction):
-0.3811148635446502
```

## Citation and Reproducing Data

If you use `LNO-MBE-CCSD(T)` in your work, please cite it as follows:

- Efficient first-principles modeling of complex molecular crystals at sub-chemical accuracy, [arXiv:2603.02180](https://arxiv.org/abs/2603.02180)

In the companion repository found at [benshi97/Data_LNOMBECC](https://github.com/benshi97/Data_LNOMBECC), we have compiled the data and outputs for **all** of the calculations in [arXiv:2603.02180](https://arxiv.org/abs/2603.02180) with detailed explanation/codes for analysing and reproducing **all** the outcomes. This can be viewed online on [Colab](https://colab.research.google.com/github/benshi97/Data_LNOMBECC/blob/master/analyse.ipynb).


## License

`LNO-MBE-CCSD(T)` is released under a [BSD 3-Clause license](https://github.com/quantum-accelerators/quacc/blob/main/LICENSE.md).
