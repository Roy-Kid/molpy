# Migrating to 0.14

molpy 0.14 pairs with `molcrafts-molrs` 0.14 (`>=0.14.0,<0.15`). The two
release together and the native core ships first; a mismatched **minor** fails
at import, a different patch does not.

Everything below is spelled `molpy`. The native types you meet (`Frame`,
`Block`, `Box`, `Atomistic`, `Trajectory`, `Element`, the force-field model)
are re-exported by identity — `mp.Frame` *is* the native class — so there is
never a reason for application code to import the core directly.

## Frame metadata

`frame.metadata` is gone. `frame.meta` is a mapping with the full `dict`
protocol:

```python
frame.meta["title"] = "water box"
if "dt" in frame.meta:
    dt = frame.meta["dt"]
for key, value in frame.meta.items():
    ...
```

There is no helper for reading or writing it; use it like a dict.

## Field names

Column and component names come from `molpy.core.fields` as constants
(`fields.X`, `fields.CHARGE`, `fields.RES_ID`, …) plus the grouped tuples
`fields.COORDS`, `fields.VELOCITIES`, `fields.DIPOLE`, `fields.QUAT`,
`fields.ENDPOINTS`. Spell a canonical field through the constant, not as a
string literal you typed yourself.

## Columns are typed and never zero-filled

`graph.column(key)` returns the column in the component's own dtype: a
float column is a zero-copy view that writes through to the graph, integer,
boolean and string columns are copies. A component that is missing on *any*
entity raises `KeyError` naming how many are missing — it is never read back
as `0`. Use `graph.validity(key)` to find the holes, `graph.get(handle, key)`
for entity-wise reads, and `graph.columns()` to list what is registered.

`atomistic.symbols` and `atomistic.xyz` are built on these reads, so a
structure whose atoms only partly carry coordinates now raises instead of
placing the rest at the origin.

## Neighbor lists

`NeighborQuery.free(points, cutoff).query_self()` (or `.query(other_points)`)
returns a `Neighbors` table. Its columns are **methods**:
`query_point_indices()`, `point_indices()`, `dist_sq()`, `disp()`
(`disp = r_j − r_i`); `n_pairs` is an attribute.

## Typifiers

`molpy.typifier` exports `Typifier`, `Match`, `OPLSAATypifier`,
`MMFFTypifier`, `ClpTypifier`, `AmberToolsTypifier`, `SmartsTypifier`,
`LocalTypifier`, `TypeScope`, `ForceFieldParams` and the region helpers.
Spellings that appeared in older docs (`OplsTypifier`, `PairTypifier`, a
Python `UFFTypifier`) do not exist. `typify()` returns a new graph.

## MD, optimization and potentials

- `molpy.md` is the native MD surface re-exported verbatim.
- `molpy.optimize.LBFGS` (also `molpy.LBFGS`) is the native L-BFGS and
  returns `(frame, OptReport)`. The Python `optimize.lbfgs` module, the
  `SoftPotential` and the `Optimizer.run(inplace=...)` wrapper are gone;
  compose typify → `to_potentials` → `LBFGS(...).run(frame)` yourself.
- `molpy.potential` re-exports the kernels; nothing is defined in molpy.

## Compute

- Shells are plain classes: construct with the measurement parameters, then
  call the `compute` verb with the data —
  `RDF(n_bins=100, r_max=10.0).compute(frames, neighbors)`. Calling the object
  itself no longer works.
- `molpy.compute.Compute` is the molrs `Protocol`
  (`molpy.compute.Compute is molrs.compute.Compute`): a class conforms by
  defining `compute(...)`, never by subclassing. The molpy base class is gone,
  and with it `dump()` and the `**config` catch-all that fed it.
- `RadicalVoronoi` and `VoronoiIntegration` are now the molrs classes
  themselves, which keep their own verbs: `RadicalVoronoi().build(...)` and
  `VoronoiIntegration().integrate(...)`.
- `Workflow` is removed: a DAG of analyses is a script, not a library object.
- `ACFAnalyzer` and `SpectralAnalyzer` are removed; `IonicConductivity` and
  `DielectricSusceptibility` stream frames and no longer print progress (the
  `progress_every` key is gone).
- `molpy.compute.spectra` classes are the native ones re-exported.

## I/O

- `read_lammps_log` returns the native `LammpsLog`: `log.runs` holds one
  `LammpsRun` per `run`, whose `thermo` is a `LammpsThermo` with
  `columns()`, `rows()`, `["Temp"]`, `"Step" in thermo`, `len(thermo)` and
  `to_dict()`. The molpy dataclasses are gone.
- `emit_all(...)` is removed; loop over `emit(name, ...)`.
- `XMLForceFieldReader` / `OPLSAAForceFieldReader` shells are removed; use
  `read_xml_forcefield` / `read_opls_xml`.
- `from molpy.io import *` exports only names that exist.

## Builders

- `build_crystal(...)` → `Lattice(...).build(...)`.
- `create_polydisperse_from_ir(ir)` → `ir.build()` on a `DistributionIR`.
- `get_forcefield_path(name)` lives in `molpy.data` and returns a `Path`; the
  `molpy.data.forcefield` copy that returned a string is gone.
- `GraphAssembler` raises on an unknown component map number instead of
  silently skipping it.

## Packing

`InsideBoxConstraint`, `OutsideBoxConstraint`, `InsideSphereConstraint` and
`OutsideSphereConstraint` now score the squared distance to the boundary and
`dpenalty` is its gradient (before: a count of violators and a unit push
direction). `MinDistanceConstraint` uses a neighbor query and also returns
the gradient.

## RDKit adapter

`RDKitAdapter` joins the two representations by the `mp_id` atom component
and RDKit atom property only. It no longer assigns `id` or `atomic_num`. An
RDKit atom with a negative `mp_id` becomes a new atom on the next
`sync_to_internal()`; an untagged RDKit atom is an error.

## Removed without replacement

`molpy.reacter`, `io.data.amber_prep`, the moltemplate `emit_all`, the
compute `Workflow`, the test-data download step (fixtures are committed under
`tests/tests-data/`).
