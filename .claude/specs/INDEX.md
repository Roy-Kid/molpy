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
| [lammps-ff-p0-single-boundary](lammps-ff-p0-single-boundary.md) | **done** (implemented) — single boundary + Frame sovereignty + `map_type` + data_coeffs R/W |

## release-0-12-molpy

| Slug | Status |
|---|---|
| [release-0-12-molpy-compute-sink](release-0-12-molpy-compute-sink.md) | **done** |
| [release-0-12-molpy-01-api-cleanup](release-0-12-molpy-01-api-cleanup.md) | **done** |
| [release-0-12-molpy-05-docs-harness](release-0-12-molpy-05-docs-harness.md) | **done** (ac-005 doc-blocks: run when env has molrs) |

**Superseded / do not implement as written:**

| Slug | Why |
|---|---|
| `release-0-12-molpy-02-units-fs` | SI units live in molrs; do not re-fix prefactors in molpy |
| `release-0-12-molpy-03-charge-transport` | no molpy charge-weighting — use molrs `GreenKuboConductivity` / `EinsteinConductivity` |
| `release-0-12-molpy-04-compute-tests` | rewrite as identity asserts against molrs types |

> Closed historically: graph-assembler-01..04 + molrs-core-cutover (2026-07-22).
