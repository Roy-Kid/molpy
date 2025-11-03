from __future__ import annotations

from math import cos, sin
from typing import Iterable

from ..assembly import Assembly
from ..entity import Entity


def _iter_scope(all_ents: Iterable[Entity], scope: set[Entity] | None) -> Iterable[Entity]:
    if scope is None:
        return all_ents
    return (e for e in all_ents if e in scope)


def _vec_add(a: list[float], b: list[float]) -> list[float]:
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def _vec_sub(a: list[float], b: list[float]) -> list[float]:
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def _vec_scale(a: list[float], s: float) -> list[float]:
    return [a[0] * s, a[1] * s, a[2] * s]


def _dot(a: list[float], b: list[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: list[float], b: list[float]) -> list[float]:
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _norm(v: list[float]) -> float:
    return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5


def _unit(v: list[float]) -> list[float]:
    n = _norm(v)
    if n == 0:
        return [0.0, 0.0, 0.0]
    return [v[0] / n, v[1] / n, v[2] / n]


def _rodrigues_rotate(p: list[float], k: list[float], angle: float, about: list[float]) -> list[float]:
    # p' = about + R (p - about)
    v = _vec_sub(p, about)
    kx, ky, kz = k
    c = cos(angle)
    s = sin(angle)
    # rotation matrix components
    r00 = c + kx * kx * (1 - c)
    r01 = kx * ky * (1 - c) - kz * s
    r02 = kx * kz * (1 - c) + ky * s

    r10 = ky * kx * (1 - c) + kz * s
    r11 = c + ky * ky * (1 - c)
    r12 = ky * kz * (1 - c) - kx * s

    r20 = kz * kx * (1 - c) - ky * s
    r21 = kz * ky * (1 - c) + kx * s
    r22 = c + kz * kz * (1 - c)

    x = r00 * v[0] + r01 * v[1] + r02 * v[2]
    y = r10 * v[0] + r11 * v[1] + r12 * v[2]
    z = r20 * v[0] + r21 * v[1] + r22 * v[2]
    return _vec_add([x, y, z], about)


class SpatialMixin:
    """Geometry operations on entities with a "pos" key only."""

    def _all_entities(self: Assembly) -> list[Entity]:
        ents: list[Entity] = []
        for ecls in self.entities.classes():
            ents.extend(self.entities.bucket(ecls).items)
        return ents

    def move(self: Assembly, delta: list[float], *, scope: set[Entity] | None = None) -> None:
        for e in _iter_scope(self._all_entities(), scope):
            pos = e.get("pos")
            if isinstance(pos, list) and len(pos) == 3:
                e["pos"] = _vec_add(pos, delta)

    def rotate(
        self: Assembly,
        axis: list[float],
        angle: float,
        about: list[float] | None = None,
        *,
        scope: set[Entity] | None = None,
    ) -> None:
        k = _unit(axis)
        o = [0.0, 0.0, 0.0] if about is None else about
        for e in _iter_scope(self._all_entities(), scope):
            pos = e.get("pos")
            if isinstance(pos, list) and len(pos) == 3:
                e["pos"] = _rodrigues_rotate(pos, k, angle, o)

    def scale(
        self: Assembly,
        factor: float,
        about: list[float] | None = None,
        *,
        scope: set[Entity] | None = None,
    ) -> None:
        o = [0.0, 0.0, 0.0] if about is None else about
        for e in _iter_scope(self._all_entities(), scope):
            pos = e.get("pos")
            if isinstance(pos, list) and len(pos) == 3:
                v = _vec_sub(pos, o)
                e["pos"] = _vec_add(o, _vec_scale(v, factor))

    def align(
        self: Assembly,
        a: Entity,
        b: Entity,
        *,
        a_dir: list[float] | None = None,
        b_dir: list[float] | None = None,
        flip: bool = False,
        scope: set[Entity] | None = None,
    ) -> None:
        pa = a.get("pos")
        pb = b.get("pos")
        if not (isinstance(pa, list) and isinstance(pb, list) and len(pa) == 3 and len(pb) == 3):
            return  # silently skip if missing positions

        ents = set(self._all_entities()) if scope is None else scope

        # rotate if directions provided
        if a_dir is not None and b_dir is not None:
            va = _unit(a_dir)
            vb = _unit(b_dir)
            if flip:
                vb = _vec_scale(vb, -1.0)
            # axis = va x vb; angle = atan2(|axis|, dot)
            axis = _cross(va, vb)
            na = _norm(axis)
            if na > 0:
                # angle via sin/cos components
                from math import atan2

                angle = atan2(na, _dot(va, vb))
                for e in _iter_scope(ents, None):
                    pos = e.get("pos")
                    if isinstance(pos, list) and len(pos) == 3:
                        e["pos"] = _rodrigues_rotate(pos, _vec_scale(axis, 1.0 / na), angle, pa)
        # translate so that a -> b
        new_pa = a.get("pos")
        if isinstance(new_pa, list) and len(new_pa) == 3:
            delta = _vec_sub(pb, new_pa)
            self.move(delta, scope=ents)
