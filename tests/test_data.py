from __future__ import annotations

from pathlib import Path

from lnombecc.data import calculation_defaults

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
