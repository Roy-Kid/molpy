"""OpenMMEmitter: XML + PDB + a starter script with the requested knobs."""

from molpy.io.emit import OpenMMEmitter


def test_script_carries_the_run_parameters(tmp_path, water, tip3p):
    paths = OpenMMEmitter().emit(
        water,
        tip3p,
        tmp_path,
        prefix="w",
        temperature_K=310.0,
        timestep_fs=2.0,
        steps=42,
    )
    assert [p.name for p in paths] == ["w.xml", "w.pdb", "w.py"]
    assert all(p.exists() for p in paths)
    script = paths[2].read_text()
    assert "w.pdb" in script and "w.xml" in script
    assert "310.0" in script and "2.0" in script and "42" in script
