from __future__ import annotations

from copy import deepcopy

COMMON = {
    "verbosity": 3,
    "unit": "angs",
    "scfalg": "locfit1",
    "scfmaxit": 200,
    "symm": "off",
    "scftol": 8,
    "basis_sm": "def2-SVP",
}

LNO_COMMON = {
    **COMMON,
    "calc": "LNO-CCSD(T)",
    "usedisk": 0,
    "aocd": "extra",
    "chtol": "1d-6",
    "iesttol": 0.001,
}

MP2_COMMON = {**COMMON, "calc": "DF-MP2"}

STEP_SETTINGS = {
    1: {"scfiguess": "small", "lcorthr": "tight"},
    2: {"scfiguess": "small", "lcorthr": "vtight"},
    3: {"scfiguess": "small"},
    4: {"scfiguess": "small"},
}

def with_basis(base: dict, basis: str) -> dict:
    d = deepcopy(base)
    d["basis"] = basis
    return d

def build_body(
    basis_map: dict[int, str | None],
    *,
    renumber: bool = False,
    extra_overrides: dict | None = None,
) -> dict[int, dict]:

    body: dict[int, dict] = {}
    items = [(step, basis) for step, basis in basis_map.items() if basis is not None]

    for new_step, (orig_step, basis) in enumerate(items, start=1):
        step_key = new_step if renumber else orig_step

        base = LNO_COMMON if orig_step in (1, 2) else MP2_COMMON
        d = with_basis(base, basis)
        d.update(STEP_SETTINGS[orig_step])

        if extra_overrides:
            d.update(extra_overrides)

        body[step_key] = d

    return body

elatt_calculation_defaults = {
    "1B": build_body({1: "aug-cc-pVQZ", 2: "aug-cc-pVQZ", 3: "aug-cc-pVQZ", 4: "aug-cc-pV5Z"}),
    "2B": build_body({1: "aug-cc-pVTZ", 2: "aug-cc-pVTZ", 3: "aug-cc-pVTZ", 4: "aug-cc-pVQZ"}),
    "3B": build_body({1: "aug-cc-pVDZ", 2: "aug-cc-pVDZ", 3: "aug-cc-pVDZ", 4: "aug-cc-pVTZ"}),
}

calculation_defaults = {
    "lattice": {
    "1B": build_body({1: "aug-cc-pVQZ", 2: "aug-cc-pVQZ", 3: "aug-cc-pVQZ", 4: "aug-cc-pV5Z"}),
    "2B": build_body({1: "aug-cc-pVTZ", 2: "aug-cc-pVTZ", 3: "aug-cc-pVTZ", 4: "aug-cc-pVQZ"}),
    "3B": build_body({1: "aug-cc-pVDZ", 2: "aug-cc-pVDZ", 3: "aug-cc-pVDZ", 4: "aug-cc-pVTZ"}),
},
    "relative": {
    "1B": build_body(
        {1: "aug'-cc-pVQZ", 2: None, 3: "aug'-cc-pVQZ", 4: "aug'-cc-pV5Z"},
        extra_overrides={"iesttol": 0.0},
    ),
    "2B": build_body(
        {1: "aug'-cc-pVTZ", 2: None, 3: "aug'-cc-pVTZ", 4: "aug'-cc-pVQZ"},
        extra_overrides={"iesttol": 0.0},
    ),
    "3B": build_body(
        {1: "aug'-cc-pVDZ", 2: None, 3: "aug'-cc-pVDZ", 4: "aug'-cc-pVTZ"},
        extra_overrides={"iesttol": 0.0},
    ),
}
}

vasp_calculation_defaults = {
    "ismear": 0,
    "sigma": 0.05,
    "nelm": 1000,
    "algo": "Normal",
    "prec": "Accurate",
    "precfock": "Normal",
    "ediff": 1e-5,
    "encut": 1000,
    "ibrion": -1,
    "hfrcut": -1,
    "kgamma": True,
    "kspacing": 0.25,
    "lhfcalc": True,
    "aexx": 1.0,
    "lwave": False,
    "lcharg": False,
}

vasp_calculation_setup = {
    "H": "_h",
    "B": "_h",
    "C": "_h",
    "N": "_h",
    "O": "_h",
    "F": "_h",
    "P": "_h",
    "S": "_h",
    "Cl": "_h",
    "Br": "_sv_GW"
}
