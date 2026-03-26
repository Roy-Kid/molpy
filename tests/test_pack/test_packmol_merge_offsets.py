"""Unit tests for Packmol frame merge/index offset logic."""

import numpy as np

from molpy.core import Block, Frame
from molpy.pack import InsideBoxConstraint, Target
from molpy.pack.packer.packmol import Packmol


def _box() -> InsideBoxConstraint:
    return InsideBoxConstraint(
        length=np.array([20.0, 20.0, 20.0]),
        origin=np.array([0.0, 0.0, 0.0]),
    )


def _optimized_frame(n_atoms: int) -> Frame:
    coords = np.arange(n_atoms * 3, dtype=float).reshape(n_atoms, 3)
    return Frame(
        {
            "atoms": Block(
                {
                    "x": coords[:, 0],
                    "y": coords[:, 1],
                    "z": coords[:, 2],
                }
            )
        }
    )


def test_build_final_frame_offsets_cumulative_with_atomi_keys() -> None:
    """Offsets must use cumulative atoms of all previous packed instances."""
    frame_a = Frame(
        {
            "atoms": Block(
                {
                    "x": np.zeros(3),
                    "y": np.zeros(3),
                    "z": np.zeros(3),
                    "id": np.array([1, 2, 3], dtype=int),
                }
            ),
            "bonds": Block(
                {
                    "id": np.array([1], dtype=int),
                    "i": np.array([0], dtype=int),
                    "j": np.array([1], dtype=int),
                }
            ),
            "angles": Block(
                {
                    "id": np.array([1], dtype=int),
                    "i": np.array([0], dtype=int),
                    "j": np.array([1], dtype=int),
                    "k": np.array([2], dtype=int),
                }
            ),
        }
    )
    frame_b = Frame(
        {
            "atoms": Block(
                {
                    "x": np.zeros(2),
                    "y": np.zeros(2),
                    "z": np.zeros(2),
                    "id": np.array([1, 2], dtype=int),
                }
            ),
            "bonds": Block(
                {
                    "id": np.array([1], dtype=int),
                    "i": np.array([0], dtype=int),
                    "j": np.array([1], dtype=int),
                }
            ),
        }
    )

    # total atoms = 2 * 3 + 1 * 2 = 8
    optimized = _optimized_frame(8)
    packer = Packmol()
    result = packer._build_final_frame(
        targets=[
            Target(frame=frame_a, number=2, constraint=_box(), name="A"),
            Target(frame=frame_b, number=1, constraint=_box(), name="B"),
        ],
        optimized_frame=optimized,
    )

    np.testing.assert_array_equal(result["atoms"]["id"], np.arange(1, 9, dtype=int))
    np.testing.assert_array_equal(
        result["atoms"]["mol"], np.array([1, 1, 1, 2, 2, 2, 3, 3], dtype=int)
    )

    # bonds: A(0-1), A(3-4), B(6-7)
    np.testing.assert_array_equal(result["bonds"]["i"], np.array([0, 3, 6], dtype=int))
    np.testing.assert_array_equal(result["bonds"]["j"], np.array([1, 4, 7], dtype=int))

    # angles only from A instances: (0,1,2), (3,4,5)
    np.testing.assert_array_equal(result["angles"]["i"], np.array([0, 3], dtype=int))
    np.testing.assert_array_equal(result["angles"]["j"], np.array([1, 4], dtype=int))
    np.testing.assert_array_equal(result["angles"]["k"], np.array([2, 5], dtype=int))


def test_build_final_frame_offsets_impropers_with_unified_keys() -> None:
    frame = Frame(
        {
            "atoms": Block(
                {
                    "x": np.zeros(4),
                    "y": np.zeros(4),
                    "z": np.zeros(4),
                    "id": np.array([1, 2, 3, 4], dtype=int),
                }
            ),
            "impropers": Block(
                {
                    "id": np.array([1], dtype=int),
                    "i": np.array([0], dtype=int),
                    "j": np.array([1], dtype=int),
                    "k": np.array([2], dtype=int),
                    "l": np.array([3], dtype=int),
                }
            ),
        }
    )

    optimized = _optimized_frame(8)
    result = Packmol()._build_final_frame(
        targets=[Target(frame=frame, number=2, constraint=_box(), name="X")],
        optimized_frame=optimized,
    )

    np.testing.assert_array_equal(result["impropers"]["i"], np.array([0, 4], dtype=int))
    np.testing.assert_array_equal(result["impropers"]["l"], np.array([3, 7], dtype=int))


def test_build_final_frame_rejects_legacy_topology_keys() -> None:
    """Legacy atomi/atomj keys are no longer accepted."""
    legacy = Frame(
        {
            "atoms": Block(
                {
                    "x": np.zeros(2),
                    "y": np.zeros(2),
                    "z": np.zeros(2),
                    "id": np.array([1, 2], dtype=int),
                }
            ),
            "bonds": Block(
                {
                    "id": np.array([1], dtype=int),
                    "atomi": np.array([0], dtype=int),
                    "atomj": np.array([1], dtype=int),
                }
            ),
        }
    )

    with np.testing.assert_raises(KeyError):
        Packmol()._build_final_frame(
            targets=[Target(frame=legacy, number=1, constraint=_box(), name="legacy")],
            optimized_frame=_optimized_frame(2),
        )
