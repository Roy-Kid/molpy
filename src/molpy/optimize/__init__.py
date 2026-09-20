"""Geometry optimization: the molrs L-BFGS minimizer.

``LBFGS(potentials, *, fmax=0.05, max_steps=500, max_step=0.2, memory=8)``
takes the :class:`molrs.ff.Potentials` compiled from a force field for the
frame under study (``forcefield.to_potentials(frame)``); ``run(frame)``
returns ``(frame, OptReport)``. Composition — typify, compile, relax — is the
caller's, exactly as with any other molrs primitive.
"""

from molrs.optimize import LBFGS, OptReport

__all__ = ["LBFGS", "OptReport"]
