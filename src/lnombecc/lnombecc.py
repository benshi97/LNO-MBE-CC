from __future__ import annotations


import json
import tempfile
import shutil
import subprocess
from pathlib import Path
from ase.db import connect
import ase.io



def create_fragments(
        poscar_filepath: str | Path,
        gas_filepath: str | Path | None = None,
        twob_cutoff: float = 12.0,
        threeb_cutoff: float = 300,
        full_uc: bool = False,
):
    """
    Create fragments using pMBE and return the ASE database.
    """

    poscar_filepath = Path(poscar_filepath)

    # Remove temporary directory if it already exists
    if Path('./tmp_pMBE').exists():
        shutil.rmtree('./tmp_pMBE')

    # Make a temporary directory to store the pMBE input and output files
    Path('./tmp_pMBE/system').mkdir(exist_ok=True,parents=True)

    # Copy the POSCAR file to the temporary directory
    shutil.copy(poscar_filepath, './tmp_pMBE/system/POSCAR')

    # Generate the .json file used in pMBE
    config = {
        "file_path": str(Path('./tmp_pMBE')),
        "crys_names": ["system"],
        "nb_term": 2,
        "cutoff": twob_cutoff,
        "mult_cutoff": threeb_cutoff,
        "full_unit_cell": full_uc,
        "out_path": str(Path('./tmp_pMBE'))
    }
    
    # Write the 2B JSON file
    json_2b_filepath = Path('./input.json')
    with open(json_2b_filepath, "w") as f:
        json.dump(config, f, indent=4)

    subprocess.run(["pmbe", str(json_2b_filepath)], check=True)
    
    # Write the 3B JSON file
    config["nb_term"] = 3
    json_3b_filepath = Path('./input.json')
    with open(json_3b_filepath, "w") as f:
        json.dump(config, f, indent=4)

    # Execute the pMBE code (assuming it's available as a command-line tool)
    subprocess.run(["pmbe", str(json_3b_filepath)], check=True)

    # Create a database to store the 1B fragments
    fragment_structure_list = ase.io.read("./tmp_pMBE/system_fragments.xyz", index=":")
    if gas_filepath is not None:
        gas_filepath = Path(gas_filepath)
        # Read the gas phase structure using ASE
        gas_structure = ase.io.read(gas_filepath)
        # Make sure gas_structure has pbc = False
        gas_structure.set_pbc(False)
        gas_structure.translate(-gas_structure.get_center_of_mass())
        gas_structure.set_cell([0, 0, 0])
    else:
        # If no gas phase structure is provided, use the first fragment as the gas phase structure
        gas_structure = fragment_structure_list[0].copy()


    with connect('system_1b.db',append=False) as db:
        db.write(gas_structure, name='gas_phase', key_value_pairs={'fragment': 'gas_phase'})
        for i, fragment_structure in enumerate(fragment_structure_list):
            db. write(fragment_structure, name=f'fragment_{i}', key_value_pairs={'fragment': f'fragment_{i}'})
            if i == 0 and not full_uc:
                break


    
    # Copy ./tmp_pMBE/*system_2b_*_filtered.db to ./system_2b_filtered.db
    src_dir = Path("./tmp_pMBE")
    pattern = "*system_2b_*_filtered.db"

    matches = list(src_dir.glob(pattern))

    src_file = matches[0]
    dst_file = Path("./system_2b.db")

    shutil.copy(src_file, dst_file)

    # Copy ./tmp_pMBE/*system_3b_*_filtered.db to ./system_3b_filtered.db
    src_dir = Path("./tmp_pMBE")
    pattern = "*system_3b_*_filtered.db"

    matches = list(src_dir.glob(pattern))

    src_file = matches[0]
    dst_file = Path("./system_3b.db")

    shutil.copy(src_file, dst_file)

    # Remove the temporary directory
    shutil.rmtree('./tmp_pMBE')

# def setup_lnombecc_inputs(
#         one_body_db_filepath: str | Path,
#         two_body_db_filepath: str | Path,
#         three_body_db_filepath: str | Path,
#         output_dir: str | Path,
# ):
#     """
#     Set up the input files for the LNO-MBE-CCSD(T) calculations.
#     """
#     pass



# def create_lnombecc_inputs


# def run_periodic_hf


# def run_lnombecc


# def parse_lnombecc_outputs