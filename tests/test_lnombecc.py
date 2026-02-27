from __future__ import annotations

from pathlib import Path

from ase.db import connect

from lnombecc.lnombecc import create_fragments

FILE_DIR = Path(__file__).parent

def test_fragment_db_lengths(tmp_path,monkeypatch):
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

    # Run function
    create_fragments(
        poscar_filepath="./POSCAR",
        gas_filepath="./geometry.in",
    )

    # --- Check 1B database ---
    with connect(tmp_path / "system_1b.db") as db1:
        n1b = len(db1)

    # --- Check 2B database ---
    with connect(tmp_path / "system_2b.db") as db2:
        n2b = len(db2)

    # --- Check 3B database ---
    with connect(tmp_path / "system_3b.db") as db3:
        n3b = len(db3)

    # ---- Assertions ----
    # Adjust these expected values to your known test system
    assert n1b == 2
    assert n2b == 27
    assert n3b == 136
