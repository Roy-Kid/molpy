import numpy as np
import xarray as xr
import molpy as mp


def make_frame(ids):
    ds = xr.Dataset({
        "id": ("index", ids),
        "x": ("index", np.arange(len(ids), dtype=float)),
        "y": ("index", np.arange(len(ids), dtype=float)),
        "z": ("index", np.arange(len(ids), dtype=float)),
    })
    return mp.Frame({"atoms": ds})


def test_basic_init():
    f = make_frame([1, 2])
    assert "atoms" in f
    assert isinstance(f["atoms"], xr.Dataset)
    assert f["atoms"].sizes["index"] == 2


def test_concat_split():
    f1 = make_frame([1, 2])
    f2 = make_frame([3, 4])
    joined = mp.Frame.from_frames([f1, f2])
    assert joined["atoms"].sizes["index"] == 4
    parts = joined.split([0, 0, 1, 1])
    assert len(parts) == 2
    assert parts[0]["atoms"].sizes["index"] == 2
    assert parts[1]["atoms"].sizes["index"] == 2
