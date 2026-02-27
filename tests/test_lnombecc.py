from __future__ import annotations

from pathlib import Path

from ase.db import connect

import pytest
import re
import lnombecc.lnombecc as mod
from lnombecc.lnombecc import create_fragments, setup_lnombecc_inputs, get_cbs_extrapolation
from lnombecc.data import calculation_defaults
from numpy.testing import assert_allclose
from pathlib import Path




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
    # Work inside temporary directory
    poscar_src = Path("data/POSCAR")  # <-- provide a small test POSCAR
    poscar = tmp_path / "POSCAR"
    poscar.write_text(poscar_src.read_text())

    # Optional gas file
    gas_src = Path("data/geometry.in")  # or remove gas argument below
    gas = tmp_path / "geometry.in"
    gas.write_text(gas_src.read_text())

    # Change working directory to tmp_path
    monkeypatch.chdir(tmp_path)

    # Run fragment creation (assumes pmbe is available in your test env)
    create_fragments(
        poscar_filepath="./POSCAR",
        gas_filepath="./geometry.in",
    )

    # --- Check databases exist and sizes ---
    with connect(tmp_path / "system_1b.db") as db1:
        n1b = len(db1)

    with connect(tmp_path / "system_2b.db") as db2:
        n2b = len(db2)

    with connect(tmp_path / "system_3b.db") as db3:
        n3b = len(db3)

    assert n1b == 2
    assert n2b == 27
    assert n3b == 136

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

def test_run_periodic_hf_mrcc(tmp_path):
    



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
