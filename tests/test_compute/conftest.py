"""Shared builders for the compute tests."""

from __future__ import annotations

import numpy as np
import pytest

import molrs

import molpy as mp


@pytest.fixture
def random_periodic_frame():
    """Factory: ``n`` uniformly random points in a cubic periodic box."""

    def build(n: int = 200, box_len: float = 12.0, seed: int = 0) -> molrs.Frame:
        rng = np.random.default_rng(seed)
        xyz = rng.uniform(0.0, box_len, size=(n, 3))
        frame = molrs.Frame()
        frame["atoms"] = {"x": xyz[:, 0], "y": xyz[:, 1], "z": xyz[:, 2]}
        frame.box = mp.Box.cubic(box_len)
        return frame

    return build


@pytest.fixture
def frame_coords_snapshot():
    """Factory: an owned (n, 3) copy of a frame's coordinates."""

    def snapshot(frame: molrs.Frame) -> np.ndarray:
        block = frame["atoms"]
        return np.column_stack([block["x"], block["y"], block["z"]]).copy()

    return snapshot
