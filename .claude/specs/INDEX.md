# Specs

> Spec files are deleted on completion (CLAUDE.md: specs are alive, not archived);
> the rows below are the history. The 12 files listed were removed 2026-09-20.

## release-0-14 (closed 2026-09-20)

The molrs-side chain (08–12) is closed; see molrs `.claude/specs/INDEX.md` and
`.claude/notes/release.md` § v0.14.0 for what shipped in-tree and the manual
release steps left. In molpy: `molpy.md` re-exports by identity, the shared
formats are molrs-backed, docs speak molpy and carry
`docs/getting-started/migration-0-14.md`. Open as its own decision: the
callable `compute.base.Compute` shells versus the molrs `Compute` Protocol
(verb unification).

## Active

| Slug | Status |
|---|---|

**api-verb-unification chain — closed 2026-09-20.** All four sub-specs landed on `ci/precommit-uv-parity` (commits 128dc7e, 1bbad0c, ce3981b and the 04 commit): verb table + `__call__` policy in architecture.md § 铁律 4; compute shells expose `compute()` with `Compute` the molrs Protocol; `molpy.pack` removed; `GraphAssembler.apply`. Routed follow-ups live in the verb table's debt list (Selector name collision, free `emit`/`build_*` functions, `_polymer.py` boundary signal, dielectric recipe classes + `from_dipole_series`, 26 forwarding shells, `/mol:map` and `/mol:bootstrap` re-runs).

| Slug | Status |
|---|---|
| `lammps-ff-p0-single-boundary` | **done** (implemented) — single boundary + Frame sovereignty + `map_type` + data_coeffs R/W |

## release-0-12-molpy

| Slug | Status |
|---|---|
| `release-0-12-molpy-compute-sink` | **done** |
| `release-0-12-molpy-01-api-cleanup` | **done** |
| `release-0-12-molpy-05-docs-harness` | **done** (ac-005 doc-blocks: run when env has molrs) |

**Superseded / do not implement as written:**

| Slug | Why |
|---|---|
| `release-0-12-molpy-02-units-fs` | SI units live in molrs; do not re-fix prefactors in molpy |
| `release-0-12-molpy-03-charge-transport` | no molpy charge-weighting — use molrs `GreenKuboConductivity` / `EinsteinConductivity` |
| `release-0-12-molpy-04-compute-tests` | rewrite as identity asserts against molrs types |

> Closed historically: graph-assembler-01..04 + molrs-core-cutover (2026-07-22).
