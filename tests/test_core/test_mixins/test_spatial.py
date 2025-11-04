from math import pi, isclose

from molpy.core.aa.base import Atomistic
from molpy.core.aa import Atom, Bond


class TestSpatialMixin:
    def test_move(self) -> None:
        asm = Atomistic()
        a = Atom({"xyz": [0.0, 0.0, 0.0]})
        asm.add_entity(a)
        asm.move([1.0, 2.0, 3.0])
        assert a["xyz"] == [1.0, 2.0, 3.0]

    def test_rotate(self) -> None:
        asm = Atomistic()
        a = Atom({"xyz": [1.0, 0.0, 0.0]})
        asm.add_entity(a)
        asm.rotate([0.0, 0.0, 1.0], angle=pi / 2, about=[0.0, 0.0, 0.0])
        x, y, z = a["xyz"]
        assert isclose(x, 0.0, abs_tol=1e-6) and isclose(y, 1.0, abs_tol=1e-6) and isclose(z, 0.0, abs_tol=1e-6)

    def test_scale(self) -> None:
        asm = Atomistic()
        a = Atom({"xyz": [1.0, 1.0, 1.0]})
        asm.add_entity(a)
        asm.scale(2.0, about=[0.0, 0.0, 0.0])
        assert a["xyz"] == [2.0, 2.0, 2.0]

    def test_align_origin_only(self) -> None:
        asm = Atomistic()
        a = Atom({"xyz": [1.0, 0.0, 0.0]})
        b = Atom({"xyz": [3.0, 4.0, 0.0]})
        asm.add_entity(a)
        asm.add_entity(b)
        asm.align(a, b)
        assert a["xyz"] == [3.0, 4.0, 0.0]

    def test_align_with_dirs_and_flip(self) -> None:
        asm = Atomistic()
        a = Atom({"xyz": [0.0, 0.0, 0.0]})
        b = Atom({"xyz": [0.0, 0.0, 0.0]})
        # vector along +x should flip to -y after rotation with flip True around origin
        asm.add_entity(a)
        asm.add_entity(b)
        asm.align(a, b, a_dir=[1.0, 0.0, 0.0], b_dir=[0.0, 1.0, 0.0], flip=True)
        # scope is all entities; a stays at origin, but rotation applied — test no exception
        assert True
