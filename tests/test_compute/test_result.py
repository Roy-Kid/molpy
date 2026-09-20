"""Dielectric result containers — Debye fit of an analytic susceptibility."""

from __future__ import annotations

import numpy as np
import pytest

from molpy.compute.result import DielectricResult


def test_dielectric_result_fit_debye():
    tau0, delta, eps_inf = 5.0, 40.0, 1.0
    n, dt = 513, 0.5
    n_pad = 2 * (n - 1)
    freq = 2.0 * np.pi * np.fft.rfftfreq(n_pad, d=dt)
    x = freq * tau0
    denom = 1.0 + x * x
    er = eps_inf + delta / denom
    ei = delta * x / denom
    er[0] = eps_inf + delta
    ei[0] = 0.0
    res = DielectricResult(
        frequency=freq,
        epsilon_real=er,
        epsilon_imag=ei,
        epsilon_static=eps_inf + delta,
        epsilon_inf=eps_inf,
        route="einstein-helfand",
        component="full",
    )
    fit = res.fit_debye()
    assert fit.tau == pytest.approx(tau0, rel=0.2)
