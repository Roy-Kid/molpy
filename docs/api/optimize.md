# Optimization

Geometry optimization with the molrs L-BFGS minimizer.

## Quick reference

| Symbol | Summary | Preferred for |
|--------|---------|---------------|
| `LBFGS` | Limited-memory BFGS over the `molrs.ff.Potentials` compiled for a frame | Geometry relaxation of small/medium structures |
| `OptReport` | Outcome record: `converged`, `final_energy`, `final_fmax`, `n_steps` | Inspecting why a run stopped |

Both are molrs types re-exported unchanged; see the
[user guide](../user-guide/08_geometry_optimization.md) for the composition
(typify → `to_potentials` → `run`).

## Related

- [Potential](potential.md) -- energy/force implementations the optimizer drives

---

## Full API

::: molpy.optimize
