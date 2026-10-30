import pytest

from ionff.geometry import SOLVENTS, load_monomer


EXPECTED = {
    "water": ["OW", "HW", "HW"],
    "methanol": {"OH": 1, "HO": 1, "CA": 1, "H": 3},
    "dme": {"OE": 1, "CA": 2, "H": 6},
    "ec": {"OC": 1, "OE": 2, "CC": 1, "CA": 2, "H": 4},
    "dmc": {"OC": 1, "OE": 2, "CC": 1, "CA": 2, "H": 6},
    "acetone": {"OC": 1, "CC": 1, "CA": 2, "H": 6},
}


@pytest.mark.parametrize("name", SOLVENTS)
def test_types(name):
    _, _, types = load_monomer(name)
    want = EXPECTED[name]
    if isinstance(want, list):
        assert types == want
    else:
        assert {t: types.count(t) for t in set(types)} == want
