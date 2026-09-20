"""LammpsEmitter: data + settings from the force field, init and run scripts."""

from molpy.io.emit import LammpsEmitter


def test_writes_the_four_files_and_derives_styles(tmp_path, water, tip3p):
    paths = LammpsEmitter().emit(water, tip3p, tmp_path, prefix="w", units="real")
    assert [p.name for p in paths] == ["w.data", "w.in.settings", "w.in.init", "w.in"]
    assert all(p.exists() for p in paths)
    init = paths[2].read_text()
    assert "units real" in init
    assert "atom_style full" in init
    assert "bond_style harmonic" in init
    assert "pair_style lj/cut" in init
    run = paths[3].read_text()
    assert "read_data w.data" in run and "include w.in.settings" in run
    assert "3 atoms" in paths[0].read_text()
