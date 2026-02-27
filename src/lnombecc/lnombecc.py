from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from ase.db import connect
import ase.io
import numpy as np

from typing import Literal
from logging import getLogger

from importlib import resources
from contextlib import contextmanager

from lnombecc.data import calculation_defaults, vasp_calculation_defaults, vasp_calculation_setup
from quacc import get_settings, change_settings, flow

from quacc.calculators.mrcc.mrcc import MrccProfile, MRCC
from quacc.calculators.mrcc.io import write_mrcc, read_mrcc_outputs
from quacc.recipes.mrcc.core import static_job

from ase.calculators.vasp import Vasp as ASE_Vasp
from quacc.recipes.vasp.core import static_job

LOGGER = getLogger(__name__)

def create_fragments(
        poscar_filepath: str | Path,
        gas_filepath: str | Path | None = None,
        twob_cutoff: float = 12.0,
        threeb_cutoff: float = 300,
        full_uc: bool = False,
        include_periodic_hf: bool = False,
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

    if include_periodic_hf:
        # Set up the periodic HF database
        with connect('system_periodic.db',append=False) as db:
            # Add vacuum to the gas structure to avoid issues with the periodic HF calculation
            vac_gas_structure = gas_structure.copy()
            vac_gas_structure.set_pbc(True)
            vac_gas_structure.center(vacuum=7.5)
            db.write(vac_gas_structure, name='gas_phase')
            db.write(ase.io.read("./tmp_pMBE/system/POSCAR"), name='periodic_system')


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
    include_periodic_hf: bool = False
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
    include_periodic_hf : bool, optional
        Whether to include the periodic HF calculation in the input setup (if False, only the fragment calculations will be set up), by default False.

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
                    struct.info["dist"] = row.data["dist"]
                    struct.info["count"] = row.data["count"]

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
                    struct.info["dist"] = row.data["dist"]
                    struct.info["count"] = row.data["count"]
                    three_body_calculators[row_idx][sub_calc][calc_key] = struct
                    if write_inputs_dir is not None:
                        Path(write_inputs_dir,"3B_calcs",folder,str(sub_calc),str(calc_key)).mkdir(exist_ok=True, parents=True)
                        write_mrcc(Path(write_inputs_dir, "3B_calcs", folder,str(sub_calc),str(calc_key)) / f"MINP", struct,inputs)

    if include_periodic_hf:
        periodic_calculators = {}
        # Set up the calculators for the periodic HF calculation
        with connect("system_periodic.db") as db:
            periodic_structure = db.get_atoms(id=2) # get the periodic structure (assuming it's the second entry in the database)
            periodic_structure.calc = ASE_Vasp(**vasp_calculation_defaults, setups=vasp_calculation_setup)
            periodic_structure.info["folder"] = "periodic_HF/crystal"

            gas_structure = db.get_atoms(id=1) # get the gas phase structure (assuming it's the first entry in the database)
            gas_structure.calc = ASE_Vasp(**vasp_calculation_defaults, setups=vasp_calculation_setup)
            gas_structure.info["folder"] = "periodic_HF/molecule"

            if write_inputs_dir is not None:
                for struct in [periodic_structure, gas_structure]:
                    struct_folder = Path(write_inputs_dir, struct.info["folder"])
                    struct_folder.mkdir(exist_ok=True, parents=True)
                    struct.calc.directory = struct_folder
                    ase.io.write(struct_folder / "POSCAR", struct)
                    struct.calc.write_input(struct)
        periodic_calculators["crystal"] = periodic_structure
        periodic_calculators["molecule"] = gas_structure
        np.save("calculators_periodic.npy", periodic_calculators, allow_pickle=True)
        
        
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
            run_mrcc_calculation(calculators_1b[row_idx][calc_key], calc_parameters, calc_folder)


    for row_idx, sub_calc_dict in calculators_2b.items():
        for sub_calc, calc_dict in sub_calc_dict.items():
            for calc_key, struct in calc_dict.items():
                calc_folder = Path(run_directory, "2B_calcs", struct.info["folder"], str(sub_calc), str(calc_key))
                calc_parameters = struct.calc.parameters
                run_mrcc_calculation(calculators_2b[row_idx][sub_calc][calc_key], calc_parameters, calc_folder)


    for row_idx, sub_calc_dict in calculators_3b.items():
        for sub_calc, calc_dict in sub_calc_dict.items():
            for calc_key, struct in calc_dict.items():
                calc_folder = Path(run_directory, "3B_calcs", struct.info["folder"], str(sub_calc), str(calc_key))
                calc_parameters = struct.calc.parameters
                run_mrcc_calculation(calculators_3b[row_idx][sub_calc][calc_key], calc_parameters, calc_folder)

def run_periodic_hf(
    run_directory: str | Path = "./LNOMBECC_calcs",
    calculators_periodic_filepath: str | Path = "calculators_periodic.npy",
):
    """
    Run the periodic HF calculations for the crystal and gas phase structures.
    
    Parameters
    ----------
    run_directory : str | Path, optional
        The directory to run the calculations in, by default "./LNOMBECC_calcs".
    calculators_periodic_filepath : str | Path, optional
        The file path to the npy file containing the calculators for the periodic HF calculations, by default "calculators_periodic.npy".

    Returns
    -------
    None
    """

    periodic_calculators = np.load(calculators_periodic_filepath, allow_pickle=True).item()
    for calc_key, struct in periodic_calculators.items():
        calc_folder = Path(run_directory, "periodic_HF", calc_key)
        calc_parameters = struct.calc.parameters

        # Change kspacing to kpts for molecule
        if calc_key == "molecule":
            calc_parameters.pop("kspacing", None)
            calc_parameters["kpts"] = (1, 1, 1)

        # Skip the calculation if "aborting loop because EDIFF is reached"
        if (calc_folder / "OUTCAR").exists():
            with open(calc_folder / "OUTCAR", "r") as f:
                if any("aborting loop because EDIFF is reached" in line for line in f):
                    LOGGER.info(f"Skipping calculation for {calc_key}.")
                    continue

        with change_settings(
            {
                "RESULTS_DIR": calc_folder,
                "CREATE_UNIQUE_DIR": False,
                "GZIP_FILES": False,
            }        
            ):
            with lnombecc_preset_path() as preset_path:
                static_job(
                    struct,
                    preset=str(preset_path),
                    **calc_parameters
                )

def analyze_lnombecc_outputs(
    run_directory: str | Path = "./LNOMBECC_calcs",
    calculators_1b_filepath: str | Path = "calculators_1b.npy",
    calculators_2b_filepath: str | Path = "calculators_2b.npy",
    calculators_3b_filepath: str | Path = "calculators_3b.npy",
    calculation_type: Literal["lattice", "relative"] = "lattice",
) -> dict[str, float] :
    """
    Analyze the outputs of the LNO-MBE-CCSD(T) calculations and compile the results into a JSON file.
    
    Parameters
    ----------
    run_directory : str | Path, optional
        The directory where the calculation outputs are stored, by default "./LNOMBECC_calcs".
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

    elatt_contributions = {
        "1B": 0,
        "2B": 0,
        "3B": 0,
    }

    if calculation_type == "lattice":
        basis_family = 'acc'
    elif calculation_type == "relative":
        basis_family = 'mixcc'

    # Analyse the 1B outputs
    calculators_1b = np.load(calculators_1b_filepath, allow_pickle=True).item()

    energies_1b = []
    for row_idx, calc_dict in calculators_1b.items():
        method_calc_energies = {calc_key: {} for calc_key in calc_dict.keys()}
        for calc_key, struct in calc_dict.items():
            calc_folder = Path(run_directory, "1B_calcs", struct.info["folder"], str(calc_key))
            output_file = Path(calc_folder, f"mrcc.out")
            if not output_file.exists():
                raise FileNotFoundError(f"Output file {output_file} does not exist.")

            with open(output_file, "r") as f:
                last_lines = f.readlines()[-5:]
                if not any("Normal termination of mrcc." in line for line in last_lines):
                    raise ValueError(f"Calculation in {calc_folder} did not finish normally. Please check the output file.")
            
            method_calc_energies[calc_key] = read_mrcc_outputs(output_file)

        mp2_cbs_energy = get_cbs_extrapolation(
            hf_X=0.0,
            corr_X=method_calc_energies[3]["mp2_corr_energy"],
            hf_Y=0.0,
            corr_Y=method_calc_energies[4]["mp2_corr_energy"],
            X_size = "QZ",
            Y_size = "5Z",
            family = basis_family,
        )[-1]
        if calculation_type == "lattice":
            calc_1_delta = method_calc_energies[1]["ccsdt_corr_energy"] - method_calc_energies[1]["mp2_corr_energy"]
            calc_2_delta = method_calc_energies[2]["ccsdt_corr_energy"] - method_calc_energies[2]["mp2_corr_energy"]
            delta_cc_mp2_energy = calc_2_delta + 0.5*(calc_2_delta - calc_1_delta)
        elif calculation_type == "relative":
            delta_cc_mp2_energy = method_calc_energies[1]["ccsdt_corr_energy"] - method_calc_energies[1]["mp2_corr_energy"]
        
        energies_1b += [mp2_cbs_energy + delta_cc_mp2_energy]
    elatt_contributions["1B"] = np.mean(energies_1b[1:]) - energies_1b[0]

    energies_2b = []
    calculators_2b = np.load(calculators_2b_filepath, allow_pickle=True).item()
    for row_idx, sub_calc_dict in calculators_2b.items():
        # print(sub_calc_dict[0][1])
        count = sub_calc_dict[0][1].info["count"]
        sub_calc_energies = {sub_calc: 0 for sub_calc in sub_calc_dict.keys()}
        for sub_calc, calc_dict in sub_calc_dict.items(): # sub_calc is the index for the different fragment calculations (e.g., dimer with ghost atoms on monomer 1, dimer with ghost atoms on monomer 2, etc.)
            method_calc_energies = {calc_key: {} for calc_key in calc_dict.keys()}
            for calc_key, struct in calc_dict.items():
                calc_folder = Path(run_directory, "2B_calcs", struct.info["folder"], str(sub_calc), str(calc_key))
                output_file = Path(calc_folder, f"mrcc.out")
                if not output_file.exists():
                    raise FileNotFoundError(f"Output file {output_file} does not exist.")

                with open(output_file, "r") as f:
                    last_lines = f.readlines()[-5:]
                    if not any("Normal termination of mrcc." in line for line in last_lines):
                        raise ValueError(f"Calculation in {calc_folder} did not finish normally. Please check the output file.")
                
                method_calc_energies[calc_key] = read_mrcc_outputs(output_file)
            mp2_cbs_energy = get_cbs_extrapolation(
                hf_X=0.0,
                corr_X=method_calc_energies[3]["mp2_corr_energy"],
                hf_Y=0.0,
                corr_Y=method_calc_energies[4]["mp2_corr_energy"],
                X_size = "TZ",
                Y_size = "QZ",
                family = basis_family,
            )[-1]
            if calculation_type == "lattice":
                calc_1_delta = method_calc_energies[1]["ccsdt_corr_energy"] - method_calc_energies[1]["mp2_corr_energy"]
                calc_2_delta = method_calc_energies[2]["ccsdt_corr_energy"] - method_calc_energies[2]["mp2_corr_energy"]
                delta_cc_mp2_energy = calc_2_delta + 0.5*(calc_2_delta - calc_1_delta)
            elif calculation_type == "relative":
                delta_cc_mp2_energy = method_calc_energies[1]["ccsdt_corr_energy"] - method_calc_energies[1]["mp2_corr_energy"]

            sub_calc_energies[sub_calc] = (mp2_cbs_energy + delta_cc_mp2_energy)
        energies_2b += [(sub_calc_energies[0] - sub_calc_energies[1] - sub_calc_energies[2])*count]
    elatt_contributions["2B"] = np.sum(energies_2b)

    energies_3b = []
    calculators_3b = np.load(calculators_3b_filepath, allow_pickle=True).item()
    for row_idx, sub_calc_dict in calculators_3b.items():
        count = sub_calc_dict[0][1].info["count"]
        sub_calc_energies = {sub_calc: 0 for sub_calc in sub_calc_dict.keys()}
        for sub_calc, calc_dict in sub_calc_dict.items(): # sub_calc is the index for the different fragment calculations (e.g., trimer with ghost atoms on monomer 1, trimer with ghost atoms on monomer 2, etc.)
            method_calc_energies = {calc_key: {} for calc_key in calc_dict.keys()}
            for calc_key, struct in calc_dict.items():
                calc_folder = Path(run_directory, "3B_calcs", struct.info["folder"], str(sub_calc), str(calc_key))
                output_file = Path(calc_folder, f"mrcc.out")
                if not output_file.exists():
                    raise FileNotFoundError(f"Output file {output_file} does not exist.")

                with open(output_file, "r") as f:
                    last_lines = f.readlines()[-5:]
                    if not any("Normal termination of mrcc." in line for line in last_lines):
                        raise ValueError(f"Calculation in {calc_folder} did not finish normally. Please check the output file.")
                
                method_calc_energies[calc_key] = read_mrcc_outputs(output_file)
            mp2_cbs_energy = get_cbs_extrapolation(
                hf_X=0.0,
                corr_X=method_calc_energies[3]["mp2_corr_energy"],
                hf_Y=0.0,
                corr_Y=method_calc_energies[4]["mp2_corr_energy"],
                X_size = "DZ",
                Y_size = "TZ",
                family = basis_family,
            )[-1]
            if calculation_type == "lattice":
                calc_1_delta = method_calc_energies[1]["ccsdt_corr_energy"] - method_calc_energies[1]["mp2_corr_energy"]
                calc_2_delta = method_calc_energies[2]["ccsdt_corr_energy"] - method_calc_energies[2]["mp2_corr_energy"]
                delta_cc_mp2_energy = calc_2_delta + 0.5*(calc_2_delta - calc_1_delta)
            elif calculation_type == "relative":
                delta_cc_mp2_energy = method_calc_energies[1]["ccsdt_corr_energy"] - method_calc_energies[1]["mp2_corr_energy"]
            sub_calc_energies[sub_calc] = (mp2_cbs_energy + delta_cc_mp2_energy)
        energies_3b += [(sub_calc_energies[0] + sub_calc_energies[4] + sub_calc_energies[5] + sub_calc_energies[6] - sub_calc_energies[1] - sub_calc_energies[2] - sub_calc_energies[3])*count]
    elatt_contributions["3B"] = np.sum(energies_3b)
    return elatt_contributions

def analyze_periodic_hf_outputs(
    run_directory: str | Path = "./LNOMBECC_calcs",
    calculators_periodic_filepath: str | Path = "calculators_periodic.npy",
) -> dict[str, float]:
    """
    Analyze the outputs of the periodic HF calculations and extract the energies.
    
    Parameters
    ----------
    run_directory : str | Path, optional
        The directory where the calculation outputs are stored, by default "./LNOMBECC_calcs".
    calculators_periodic_filepath : str | Path, optional
        The file path to the npy file containing the calculators for the periodic HF calculations, by default "calculators_periodic.npy".

    Returns
    -------
    dict[str, float]
        A dictionary containing the periodic HF energies for the crystal and gas phase structures.
    """

    periodic_calculators = np.load(calculators_periodic_filepath, allow_pickle=True).item()
    energies = {}
    for calc_key, struct in periodic_calculators.items():
        calc_folder = Path(run_directory, "periodic_HF", calc_key)
        output_file = Path(calc_folder, f"OUTCAR")
        if not output_file.exists():
            raise FileNotFoundError(f"Output file {output_file} does not exist.")

        with open(output_file, "r") as f:
            if not any("aborting loop because EDIFF is reached" in line for line in f.readlines()):
                raise ValueError(f"Calculation for {calc_key} did not finish normally. Please check the output file.")
        
        # Extract the final energy from the OUTCAR file (look for "free  energy   TOTEN" in the OUTCAR)
        with open(output_file, "r") as f:
            # search in reverse
            for line in reversed(list(f)):
                if "energy  without entropy=" in line:
                    energy = float(line.split()[-1])
                    energies[calc_key] = energy
                    break
    
    num_monomers = float(len(periodic_calculators["crystal"].get_chemical_symbols())) / float(len(periodic_calculators["molecule"].get_chemical_symbols()))
    periodic_hf_elatt = energies["crystal"]/num_monomers - energies["molecule"]
    return periodic_hf_elatt


def get_cbs_extrapolation(
    hf_X: float,
    corr_X: float,
    hf_Y: float,
    corr_Y: float,
    X_size: Literal["DZ", "TZ", "QZ"] = "DZ",
    Y_size: Literal["TZ", "QZ", "5Z"] = "TZ",
    family: Literal["def2", "cc", "acc", "mixcc"] = "mixcc",
) -> tuple[float, float, float]:
    """
    Function to perform basis set extrapolation of HF and correlation energies for both the cc-pVXZ and def2-XZVP basis sets

    Parameters
    ----------
    hf_X : float
        HF energy in X basis set
    corr_X : float
        Correlation energy in X basis set
    hf_Y : float
        HF energy in Y basis set where Y = X+1 cardinal zeta number
    corr_Y : float
        Correlation energy in Y basis set
    X_size : str
        Cardinal zeta number of X basis set
    Y_size : str
        Cardinal zeta number of Y basis set
    family : str
        Basis set family. Options are `cc`, `def2`, `acc`, and `mixcc`. Where cc is for non-augmented correlation consistent basis sets, def2 is for def2 basis sets, acc is for augmented correlation consistent basis sets while mixcc is for mixed augmented + non-augmented correlation consistent basis sets

    Returns
    -------
    hf_cbs : float
        HF CBS energy
    corr_cbs : float
        Correlation CBS energy
    tot_cbs : float
        Total CBS energy
    """

    # Dictionary of alpha parameters followed by beta parameters in CBS extrapoation. Refer to: Neese, F.; Valeev, E. F. Revisiting the Atomic Natural Orbital Approach for Basis Sets: Robust Systematic Basis Sets for Explicitly Correlated and Conventional Correlated Ab Initio Methods. J. Chem. Theory Comput. 2011, 7 (1), 33-43. https://doi.org/10.1021/ct100396y.
    alpha_dict = {
        "def2_2_3": 10.39,
        "def2_3_4": 7.88,
        "cc_2_3": 4.42,
        "cc_3_4": 5.46,
        "cc_4_5": 5.46,
        "acc_2_3": 4.30,
        "acc_3_4": 5.79,
        "acc_4_5": 5.79,
        "mixcc_2_3": 4.36,
        "mixcc_3_4": 5.625,
        "mixcc_4_5": 5.625,
    }

    beta_dict = {
        "def2_2_3": 2.40,
        "def2_3_4": 2.97,
        "cc_2_3": 2.46,
        "cc_3_4": 3.05,
        "cc_4_5": 3.05,
        "acc_2_3": 2.51,
        "acc_3_4": 3.05,
        "acc_4_5": 3.05,
        "mixcc_2_3": 2.485,
        "mixcc_3_4": 3.05,
        "mixcc_4_5": 3.05,
    }

    size_to_num = {"DZ": 2, "TZ": 3, "QZ": 4, "5Z": 5}

    X = size_to_num[X_size]
    Y = size_to_num[Y_size]

    # Check if X and Y are consecutive cardinal zeta numbers
    if Y != X + 1:
        raise ValueError("The cardinal number of Y does not equal X+1")

    # Get the corresponding alpha and beta parameters depending on the basis set family
    alpha = alpha_dict[f"{family}_{X}_{Y}"]
    beta = beta_dict[f"{family}_{X}_{Y}"]

    # Perform CBS extrapolation for HF and correlation components
    hf_cbs = hf_X - np.exp(-alpha * np.sqrt(X)) * (hf_Y - hf_X) / (
        np.exp(-alpha * np.sqrt(Y)) - np.exp(-alpha * np.sqrt(X))
    )
    corr_cbs = (X ** (beta) * corr_X - Y ** (beta) * corr_Y) / (
        X ** (beta) - Y ** (beta)
    )

    return hf_cbs, corr_cbs, (hf_cbs + corr_cbs)

def cleanup_folder(folder: str | Path = "."):
    folder = Path(folder)

    patterns = [
        "55*", "56*", "fort.*", "DF*", "MOCOEF*", "MOLDEN*",
        "OCCUP", "OEINT", "OLDMO", "OSVFILE", "PRINT",
        "S12MAT*", "SCFDEN*", "SCHOL", "SROOT", "SYMTRA",
        "TEDAT", "VARS", "localcc.restart", "KEYWD",
        "iface", "FOCK", "EXIT", "DAO", "COORD.xyz",
        "pt.rst", "TEINT", "ccsd.rst",
    ]

    for pattern in patterns:
        for path in folder.glob(pattern):
            if path.is_file():
                path.unlink()

@flow
def run_mrcc_calculation(struct, calc_parameters, calc_folder):
    output_file = Path(calc_folder, f"mrcc.out")
    if output_file.exists():
        with open(output_file, "r") as f:
            last_lines = f.readlines()[-5:]
            if any("Normal termination of mrcc." in line for line in last_lines):
                LOGGER.info(f"Calculation in {calc_folder} is already finished. Skipping.")
                return
    with change_settings(
        {
            "RESULTS_DIR": calc_folder,
            "CREATE_UNIQUE_DIR": False,
            "GZIP_FILES": False,
        }        
        ):
        static_job(struct, **calc_parameters)
    cleanup_folder(calc_folder)


@contextmanager
def lnombecc_preset_path():
    ref = resources.files("lnombecc.presets").joinpath("lnombecc.yaml")
    with resources.as_file(ref) as path:
        yield Path(path)