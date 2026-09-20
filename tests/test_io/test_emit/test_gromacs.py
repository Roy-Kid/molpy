"""GromacsEmitter: gro + top from the force field, plus the two mdp files."""

from molpy.io.emit import GromacsEmitter


def test_writes_coordinates_topology_and_mdps(tmp_path, water, tip3p):
    paths = GromacsEmitter().emit(
        water, tip3p, tmp_path, prefix="w", temperature_K=280.0
    )
    assert [p.name for p in paths] == ["w.gro", "w.top", "em.mdp", "nvt.mdp"]
    assert all(p.exists() for p in paths)
    assert paths[0].read_text().splitlines()[1].strip() == "3"
    assert "[ defaults ]" in paths[1].read_text()
    assert "integrator      = steep" in paths[2].read_text()
    assert "ref_t           = 280.0" in paths[3].read_text()
