"""XMLEmitter: the force field as XML plus the coordinates as PDB."""

from molpy.io.emit import XMLEmitter


def test_writes_the_xml_and_pdb_pair(tmp_path, water, tip3p):
    paths = XMLEmitter().emit(water, tip3p, tmp_path, prefix="w")
    assert [p.name for p in paths] == ["w.xml", "w.pdb"]
    assert all(p.exists() for p in paths)
    assert "<ForceField" in paths[0].read_text()
    assert (
        sum(
            line.startswith(("ATOM", "HETATM"))
            for line in paths[1].read_text().splitlines()
        )
        == 3
    )
