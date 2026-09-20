# Performance Standards

MolPy-specific NumPy / algorithm performance rules. Consumed by `/mol:review --axis=perf`
and the `optimizer` agent. Migrated from the former local `molpy-perf` skill and
`molpy-optimizer` agent (2026-06-10).

## Hot paths (profile these first)

- **Compute operators** (RDF, MSD, dielectric, …): called per-frame, must be fast.
- **Pairwise distance**: the most common inner loop across compute/ and builder/.
- **Builder placement**: coordinate transforms in tight loops.
- **Parser / SMARTS**: SMILES/SMARTS run in molrs — avoid re-compiling the same
  SMARTS pattern in a hot loop; reuse `SmartsPattern` instances.

## Vectorization rules

- No Python for-loops over atoms, bonds, or frames — use NumPy vectorization.
- Per-atom distances: `np.linalg.norm(r, axis=1)`, not a loop.
- Per-pair operations: broadcasting or `scipy.spatial.distance`.
- `np.einsum` for complex tensor contractions/reductions.
- `np.add.at` for scatter operations.
- `np.empty` + fill instead of `np.zeros` when values are immediately overwritten.

## Memory rules

- No unnecessary `.copy()` of large coordinate arrays; prefer views (boolean masks)
  over fancy indexing where possible.
- Never materialize full pairwise distance matrices for large systems.
- Stream large trajectories — never load-all-into-memory.
- `del` large intermediates in multi-step calculations.
- `np.float32` for coordinates is acceptable when float64 precision is not needed.
- Contiguous memory layout for iteration-heavy arrays.

## Algorithm complexity

- Neighbor search: O(N) cell lists or KD-tree, never O(N²) all-pairs for large N.
- Topology/graph algorithms: use the molrs graph kernels (`Atomistic` topology,
  `NeighborQuery`), do not hand-roll traversals; cache repeated traversals.
- Avoid repeated sorts; maintain sorted invariants or cache sort results.
- Document complexity in docstrings: O(N), O(N²), etc.

## I/O

- Text parsing for large files is a bottleneck — prefer binary formats (HDF5).
- No repeated file open/close in loops; buffer large writes.

## Discipline

- Never sacrifice correctness for speed. There is no benchmark harness in this
  repo (the benchmark/regression system is being redesigned); a performance
  change is verified for correctness by the unit suite and its speed claim is
  recorded here as owed, not asserted.

## Profiling commands

```bash
python -m cProfile -o profile.out script.py
kernprof -l -v script.py                      # pip install line_profiler
python -m memory_profiler script.py           # pip install memory_profiler
```

## Owed (2026-09-20)

Hot-path findings from the 0.14 cleanup, recorded rather than fixed:

- `builder/assembly/_proximity.py` — without a cutoff the site pairing is an
  O(N_a×N_b) Python double loop, `xyz` is read with three PyO3 calls per atom,
  and adjacency / connected components are hand-rolled although molrs has them.
- `builder/assembly/_placer.py` — every `place()` rescans `RES_ID` over the
  whole world and reads three coordinates per atom through PyO3.
- `builder/assembly/_assembler.py` — `_total_charge` is a per-atom Python scan
  twice per `assemble()`. Blocked on molrs: `Atomistic.column(name)` fills a
  missing component with 0 instead of failing, so a column read cannot tell
  "no charge column" from "all zero".
- `core/atomistic.py` — `symbols` reads per atom; `def_*s` call Rust once per
  element (needs a molrs batch entry); `select` / `get_neighbors` scan O(E).
- `adapter/rdkit.py` — ≥ 12 full passes over `atomistic.atoms`, two of them
  nested (optional dependency; not exercised by the gate).
- `typifier/clp.py` — `clp.xml` (364 KB) is parsed twice at first use (once by
  ElementTree to strip bonded sections for the molrs SMARTS typifier, once by
  molrs for the force field); one-time cost behind `lru_cache`.
- `pack/constraint.py` — `InsideBoxConstraint` / sphere `dpenalty` return a
  push direction, not the gradient the base class documents;
  `MinDistanceConstraint` now returns the gradient.
- `compute/dielectric.py` — `DielectricSusceptibility` prints progress from
  library code.
