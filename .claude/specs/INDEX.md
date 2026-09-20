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
| [api-verb-unification-02-compute](api-verb-unification-02-compute.md) | **approved** — 删 compute/base.py，31 壳 `__call__`→`compute`，`Compute` = molrs Protocol 身份再导出，voronoi 两壳改 identity re-export；depends_on 01 |
| [api-verb-unification-03-pack](api-verb-unification-03-pack.md) | **approved** — 整包删除 molpy.pack（零消费者，打包归 molpack）；depends_on 01 |
| [api-verb-unification-04-assembler](api-verb-unification-04-assembler.md) | **approved** — `GraphAssembler.assemble`→`apply` + 迁移页动词表；depends_on 01 |

**api-verb-unification chain — five rulings approved by the maintainer 2026-09-20; 01 landed (verb table in architecture.md § 铁律 4), 02–04 in progress:**
(1) `assemble → apply`（不是 `build`）；(2) 整包删除 `molpy.pack`（含 `constraint.py`/`target.py`）；(3) compute 壳不留 `__call__` 糖、`dump()`/`**config` 随基类消失（breaking）；(4) voronoi 两壳改 identity re-export 而非改名；(5) `tests/test_compute/test_dielectric.py` 保留并路由（不删）。全部裁定即 `/mol:impl-all api-verb-unification`；任何一项改判则先 `/mol:spec` supersede 对应 sub-spec。

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
