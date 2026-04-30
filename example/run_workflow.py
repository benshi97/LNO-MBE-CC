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
    twob_cutoff=8.0,            # cutoff distance for the 2B terms
    threeb_cutoff=150.0,         # cutoff distance for the 3B terms
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
    memory="200GB",               # memory for the MRCC calculations
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
