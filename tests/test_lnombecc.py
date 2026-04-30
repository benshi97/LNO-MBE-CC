from __future__ import annotations

import re
import shutil
from pathlib import Path

import ase.io
import numpy as np
import pytest
from ase.db import connect
from numpy.testing import assert_allclose

import lnombecc.lnombecc as mod
from lnombecc.data import calculation_defaults
from lnombecc.lnombecc import (
    analyze_lnombecc_outputs,
    analyze_periodic_hf_outputs,
    create_fragments,
    get_cbs_extrapolation,
    run_lnombecc,
    run_periodic_hf,
    setup_lnombecc_inputs,
)

FILE_DIR = Path(__file__).parent

def test_calculation_defaults():

    # Top-level keys
    assert "lattice" in calculation_defaults
    assert "relative" in calculation_defaults

    print(calculation_defaults)

    # Lattice has all 4 steps
    lattice_1b = calculation_defaults["lattice"]["1B"]
    assert set(lattice_1b.keys()) == {1, 2, 3, 4}
    assert lattice_1b[1]["iesttol"] == 0.001

    # Relative skips step 2
    relative_1b = calculation_defaults["relative"]["1B"]
    assert 2 not in relative_1b

    # Relative overrides iesttol
    for step, settings in relative_1b.items():
        if settings["calc"] == "LNO-CCSD(T)":
            assert settings["iesttol"] == 0.0





def test_create_and_write(tmp_path, monkeypatch):
    poscar = tmp_path / "POSCAR"
    poscar.write_text((FILE_DIR / "data" / "POSCAR").read_text())

    gas = tmp_path / "geometry.in"
    gas.write_text((FILE_DIR / "data" / "geometry.in").read_text())

    pmbe_calls = []

    def fake_run_pmbe(command, check):
        assert command[0] == "pmbe"
        assert check is True
        pmbe_calls.append(command)

        pmbe_dir = tmp_path / "tmp_pMBE"
        with connect(FILE_DIR / "data" / "system_1b.db") as db:
            ase.io.write(pmbe_dir / "system_fragments.xyz", [db.get_atoms(id=2)])

        shutil.copy(
            FILE_DIR / "data" / "system_2b.db",
            pmbe_dir / "test_system_2b_0_filtered.db",
        )
        shutil.copy(
            FILE_DIR / "data" / "system_3b.db",
            pmbe_dir / "test_system_3b_0_filtered.db",
        )

    monkeypatch.setattr(mod.subprocess, "run", fake_run_pmbe)

    # Change working directory to tmp_path
    monkeypatch.chdir(tmp_path)

    create_fragments(
        poscar_filepath=poscar,
        gas_filepath=gas,
    )
    assert len(pmbe_calls) == 2

    # --- Check databases exist and sizes ---
    with connect(tmp_path / "system_1b.db") as db1:
        n1b = len(db1)

    with connect(tmp_path / "system_2b.db") as db2:
        n2b = len(db2)

    with connect(tmp_path / "system_3b.db") as db3:
        n3b = len(db3)

    assert n1b == 2
    assert n2b == 9
    assert n3b == 27

    # -------------------------------------------------------------------------
    # Patch write_mrcc so we can reliably inspect the written MINP contents
    # (and not depend on MRCC/quacc formatting details in unit tests).
    # -------------------------------------------------------------------------
    def fake_write_mrcc(minp_path: Path, struct, inputs: dict):
        # Minimal deterministic content for assertions
        lines = [
            f"calc={inputs.get('calc')}",
            f"basis={inputs.get('basis')}",
            f"mem={inputs.get('mem')}",
        ]
        # include ghost if present (2B/3B)
        if "ghost" in inputs:
            lines.append("ghost=" + inputs["ghost"].replace("\n", "\\n"))
        minp_path.write_text("\n".join(lines) + "\n")

    monkeypatch.setattr(mod, "write_mrcc", fake_write_mrcc)

    # Now write inputs
    outdir = tmp_path / "LNOMBECC_calcs"
    setup_lnombecc_inputs(write_inputs_dir=outdir)

    # -------------------------------------------------------------------------
    # Read ONE MINP file and assert contents
    # We'll pick: 1B_calcs / first folder / calc_key=1 / MINP
    # Folder zero-padding depends on DB size; your code uses 0..N-1 with width.
    # With n1b==2, width becomes 2 and folders are "00", "01".
    # -------------------------------------------------------------------------
    minp = outdir / "1B_calcs" / "00" / "1" / "MINP"
    assert minp.exists(), f"Expected MINP to exist at {minp}"

    text = minp.read_text()
    # Basic sanity checks on what was written
    assert "calc=" in text
    assert "basis=" in text
    assert "mem=" in text
    # Optional: make these stricter if you know the expected defaults for step 1
    # e.g. assert "calc=LNO-CCSD(T)" in text

    # Check that no folder called periodic_HF
    # Assert periodic_HF folder does NOT exist
    assert not (outdir / "periodic_HF").exists()

def test_run_lnombecc_periodic_hf():

    run_periodic_hf(
        run_directory=FILE_DIR / "data/LNOMBECC_calcs",
        calculators_periodic_filepath=FILE_DIR / "data/calculators_periodic.npy"
    )

    run_lnombecc(
        run_directory=FILE_DIR / "data/LNOMBECC_calcs",
        calculators_1b_filepath=FILE_DIR / "data/calculators_1b.npy",
        calculators_2b_filepath=FILE_DIR / "data/calculators_2b.npy",
        calculators_3b_filepath=FILE_DIR / "data/calculators_3b.npy"
    )

def test_analyze_lnombecc_periodic_hf():

    periodic_hf_elatt = analyze_periodic_hf_outputs(
        run_directory=FILE_DIR / "data/LNOMBECC_calcs",
        calculators_periodic_filepath=FILE_DIR / "data/calculators_periodic.npy")

    lnombecc_elatt = analyze_lnombecc_outputs(
        run_directory=FILE_DIR / "data/LNOMBECC_calcs",
        calculators_1b_filepath=FILE_DIR / "data/calculators_1b.npy",
        calculators_2b_filepath=FILE_DIR / "data/calculators_2b.npy",
        calculators_3b_filepath=FILE_DIR / "data/calculators_3b.npy"
    )

    assert np.isclose(periodic_hf_elatt, -0.093841577500001, rtol=1e-5, atol=1e-7)
    assert np.isclose(lnombecc_elatt["1B"], -0.010775053988799854, rtol=1e-5, atol=1e-7)
    assert np.isclose(lnombecc_elatt["2B"], -0.29135135812023716, rtol=1e-5, atol=1e-7)
    assert np.isclose(lnombecc_elatt["3B"], 0.014853126064387823, rtol=1e-5, atol=1e-7)



def test_get_cbs_extrapolation():
    cbs_extrapolated_ene = get_cbs_extrapolation(
        0.05710643952079408,
        -0.0899534679109486,
        0.027170816450961865,
        -0.12719774638527426,
        X_size="DZ",
        Y_size="TZ",
        family="mixcc",
    )

    assert_allclose(
        (0.017185301080196724, -0.1486152208466108, -0.13142991976641408),
        cbs_extrapolated_ene,
        rtol=1e-05,
        atol=1e-07,
    )

    with pytest.raises(
        ValueError, match=re.escape("The cardinal number of Y does not equal X+1")
    ):
        get_cbs_extrapolation(
            0.05710643952079408,
            -0.0899534679109486,
            0.027170816450961865,
            -0.12719774638527426,
            family="mixcc",
            X_size="TZ",
            Y_size="DZ",
        )
