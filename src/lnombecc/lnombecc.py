from __future__ import annotations


import json
import tempfile
import shutil
import subprocess
from pathlib import Path
from ase.db import connect
import ase.io
import numpy as np

from typing import Literal

from lnombecc.data import calculation_defaults
from quacc import get_settings, change_settings

from quacc.calculators.mrcc.mrcc import MrccProfile, MRCC
from quacc.calculators.mrcc.io import write_mrcc, read_mrcc_outputs
from quacc.recipes.mrcc.core import static_job




def create_fragments(
        poscar_filepath: str | Path,
        gas_filepath: str | Path | None = None,
        twob_cutoff: float = 12.0,
        threeb_cutoff: float = 300,
        full_uc: bool = False,
):
    """
    Create fragments using pMBE to create the ASE database for the 1B, 2B and 3B fragments.
    
    Parameters
    ----------
    poscar_filepath : str | Path
        The file path to the POSCAR file of the periodic system.
    gas_filepath : str | Path | None, optional
        The file path to the gas phase structure (if None, the first fragment will be used as the gas phase structure), by default None.
    twob_cutoff : float, optional
        The cutoff distance for the 2B fragments in Angstrom, by default 12.0.
    threeb_cutoff : float, optional
        The cutoff distance for the 3B fragments in Angstrom, by default 300.
    full_uc : bool, optional
        Whether to use the full unit cell for the fragments (if False, only the first fragment will be used), by default False.
    
    Returns
    -------
    None
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

def setup_lnombecc_inputs(
    one_body_db_filepath: str | Path = "./system_1b.db",
    two_body_db_filepath: str | Path = "./system_2b.db",
    three_body_db_filepath: str | Path = "./system_3b.db",
    write_inputs_dir: str | Path = "./LNOMBECC_calcs",
    calculation_type: Literal["lattice", "relative"] = "lattice",
    memory: str = "46GB",
):
    """
    Set up input files for LNO-MBE-CCSD(T) calculations (Composite scheme).
    
    Parameters
    ----------
    one_body_db_filepath : str | Path, optional
        The file path to the ASE database containing the 1-body fragments, by default "./system_1b.db".
    two_body_db_filepath : str | Path, optional
        The file path to the ASE database containing the 2-body fragments, by default "./system_2b.db".
    three_body_db_filepath : str | Path, optional
        The file path to the ASE database containing the 3-body fragments, by default "./system_3b.db".
    write_inputs_dir : str | Path | None, optional
        The directory to write the input files for the calculations (if None, the input files will not be written), by default None.
    calculation_type : Literal["lattice", "relative"], optional
        The type of calculation to set up the inputs for, by default "lattice". This determines which set of calculation defaults to use from the `calculation_defaults` dictionary in `data.py`.
    memory : str, optional
        The amount of memory to use for the calculations (e.g., "46GB"), by default "46GB".

    Returns
    -------
    None
    """

    calc_type_defaults = calculation_defaults[calculation_type]

    one_body_calculators = {}
    with connect(one_body_db_filepath) as db:
        # Determine number of leading zeros for file naming based on number of rows in the database
        num_rows = db.count()
        width = max(1, len(str(max(num_rows - 1, 0)))) + 1

        for row_idx, row in enumerate(db.select()):
            folder = f"{row_idx:0{width}d}"
            one_body_calculators[row_idx] = {}
            for calc_key, calc_settings in calc_type_defaults["1B"].items():
                # let per-step mem override the function default if present
                inputs = {**calc_settings, "mem": calc_settings.get("mem", memory)}

                struct = row.toatoms().copy()
                struct.calc = MRCC(
                    profile=MrccProfile(command=get_settings().MRCC_CMD),
                    **inputs,
                )
                struct.info["folder"] = folder
                one_body_calculators[row_idx][calc_key] = struct
                if write_inputs_dir is not None:
                    Path(write_inputs_dir,"1B_calcs",folder, str(calc_key)).mkdir(exist_ok=True, parents=True)
                    write_mrcc(Path(write_inputs_dir, "1B_calcs", folder, str(calc_key)) / f"MINP", struct,inputs)


    
    two_body_calculators = {}
    with connect(two_body_db_filepath) as db:
        num_rows = db.count()
        width = max(1, len(str(max(num_rows - 1, 0)))) + 1
        for row_idx, row in enumerate(db.select()):
            folder = f"{row_idx:0{width}d}"
            two_body_calculators[row_idx] = {0: {}, 1: {}, 2: {}}
            dimer_length = len(row.toatoms())
            monomer_length = int(dimer_length/2)
            ghost_atoms_input = {
                0: "serialno\n\n",
                1: f"serialno\n1-{monomer_length}\n",
                2: f"serialno\n{monomer_length+1}-{dimer_length}\n",
            }
            for sub_calc in [0, 1, 2]:
                for calc_key, calc_settings in calc_type_defaults["2B"].items():
                    inputs = {
                        **calc_settings,
                        "mem": memory,
                        "ghost": ghost_atoms_input[sub_calc]
                    }
                    struct = row.toatoms().copy()
                    struct.calc = MRCC(
                        profile=MrccProfile(command=get_settings().MRCC_CMD),
                        **inputs,
                    )
                    struct.info["folder"] = folder

                    two_body_calculators[row_idx][sub_calc][calc_key] = struct
                    if write_inputs_dir is not None:
                        Path(write_inputs_dir,"2B_calcs",folder,str(sub_calc),str(calc_key)).mkdir(exist_ok=True, parents=True)
                        write_mrcc(Path(write_inputs_dir, "2B_calcs", folder,str(sub_calc),str(calc_key)) / f"MINP", struct,inputs)

    three_body_calculators = {}
    with connect(three_body_db_filepath) as db:
        num_rows = db.count()
        width = max(1, len(str(max(num_rows - 1, 0)))) + 1
        for row_idx, row in enumerate(db.select()):
            folder = f"{row_idx:0{width}d}"
            three_body_calculators[row_idx] = {0: {}, 1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
            trimer_length = len(row.toatoms())
            monomer_length = int(trimer_length/3)
            dimer_length = int(2*trimer_length/3)
            ghost_atoms_input = {
                0: "serialno\n\n",
                1: f"serialno\n1-{monomer_length}\n",
                2: f"serialno\n{monomer_length+1}-{dimer_length}\n",
                3: f"serialno\n{dimer_length+1}-{trimer_length}\n",
                4: f"serialno\n1-{dimer_length}\n",
                5: f"serialno\n{monomer_length+1}-{trimer_length}\n",
                6: f"serialno\n1-{monomer_length},{dimer_length+1}-{trimer_length}\n",
            }
            for sub_calc in [0, 1, 2, 3, 4, 5, 6]:
                for calc_key, calc_settings in calc_type_defaults["3B"].items():
                    inputs = {
                        **calc_settings,
                        "mem": memory,
                        "ghost": ghost_atoms_input[sub_calc]
                    }
                    struct = row.toatoms().copy()
                    struct.calc = MRCC(
                        profile=MrccProfile(command=get_settings().MRCC_CMD),
                        **inputs,
                    )
                    struct.info["folder"] = folder

                    three_body_calculators[row_idx][sub_calc][calc_key] = struct
                    if write_inputs_dir is not None:
                        Path(write_inputs_dir,"3B_calcs",folder,str(sub_calc),str(calc_key)).mkdir(exist_ok=True, parents=True)
                        write_mrcc(Path(write_inputs_dir, "3B_calcs", folder,str(sub_calc),str(calc_key)) / f"MINP", struct,inputs)
    
    # Save the calculators to an npy file
    np.save("calculators_1b.npy", one_body_calculators, allow_pickle=True)
    np.save("calculators_2b.npy", two_body_calculators, allow_pickle=True)
    np.save("calculators_3b.npy", three_body_calculators, allow_pickle=True)

def run_lnombecc(
    run_directory: str | Path = "./LNOMBECC_calcs",
    calculators_1b_filepath: str | Path = "calculators_1b.npy",
    calculators_2b_filepath: str | Path = "calculators_2b.npy",
    calculators_3b_filepath: str | Path = "calculators_3b.npy",
    ):
    """
    Run the LNO-MBE-CCSD(T) calculations for the 1-body, 2-body and 3-body fragments.
    
    Parameters
    ----------
    calculators_1b_filepath : str | Path, optional
        The file path to the npy file containing the calculators for the 1-body fragments, by default "calculators_1b.npy".
    calculators_2b_filepath : str | Path, optional
        The file path to the npy file containing the calculators for the 2-body fragments, by default "calculators_2b.npy".
    calculators_3b_filepath : str | Path, optional
        The file path to the npy file containing the calculators for the 3-body fragments, by default "calculators_3b.npy".
        
    Returns
    -------
    None
    """
    calculators_1b = np.load(calculators_1b_filepath, allow_pickle=True).item()
    calculators_2b = np.load(calculators_2b_filepath, allow_pickle=True).item()
    calculators_3b = np.load(calculators_3b_filepath, allow_pickle=True).item()

    for row_idx, calc_dict in calculators_1b.items():
        for calc_key, struct in calc_dict.items():
            calc_folder = Path(run_directory, "1B_calcs", struct.info["folder"], str(calc_key))
            calc_parameters = struct.calc.parameters
            with change_settings(
                {
                    "RESULTS_DIR": calc_folder,
                    "CREATE_UNIQUE_DIR": False,
                    "GZIP_FILES": False,
                }                
                ):
                static_job(calculators_1b[row_idx][calc_key], **calc_parameters)

    # for row_idx, sub_calc_dict in calculators_2b.items():
    #     for sub_calc, calc_dict in sub_calc_dict.items():
    #         for calc_key, struct in calc_dict.items():
    #             calc_folder = Path(run_directory, "2B_calcs", struct.info["folder"], str(sub_calc))
    #             calc_parameters = struct.calc.parameters
    #             with change_settings(
    #                 {
    #                     "RESULTS_DIR": calc_folder,
    #                     "CREATE_UNIQUE_DIR": False,
    #                     "GZIP_FILES": False,
    #                 }                
    #                 ):
    #                 static_job(calculators_2b[row_idx][sub_calc][calc_key], **calc_parameters)

    # for row_idx, sub_calc_dict in calculators_3b.items():
    #     for sub_calc, calc_dict in sub_calc_dict.items():
    #         # If calculations are finished, skip
    #         calc_finished = True
    #         for calc_key, struct in calc_dict.items():
    #             calc_folder = Path(run_directory, "3B_calcs", struct.info["folder"], str(sub_calc))
    #             # Check the "Normal termination" string in the output file to determine if the calculation is finished
    #             output_file = Path(calc_folder, f"mrcc_{calc_key}.out")
    #             if output_file.exists():
    #                 with open(output_file, "r") as f:
    #                     # Read last 5 lines of the output file to check for "Normal termination of mrcc."
    #                     last_lines = f.readlines()[-5:]
    #                     if any("Normal termination of mrcc." in line for line in last_lines):
    #                         print(f"Calculation in {calc_folder} is already finished. Skipping.")
    #                         continue
    #             else:
    #                 calc_finished = False
    #                 break




    #         for calc_key, struct in calc_dict.items():
    #             calc_folder = Path(run_directory, "3B_calcs", struct.info["folder"], str(sub_calc))
    #             calc_parameters = struct.calc.parameters
    #             with change_settings(
    #                 {
    #                     "RESULTS_DIR": calc_folder,
    #                     "CREATE_UNIQUE_DIR": False,
    #                     "GZIP_FILES": False,
    #                 }                
    #                 ):
    #                 static_job(calculators_3b[row_idx][sub_calc][calc_key], **calc_parameters)



# def run_periodic_hf

# def parse_lnombecc_outputs