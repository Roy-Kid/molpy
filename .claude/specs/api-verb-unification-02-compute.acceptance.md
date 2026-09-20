---
spec: api-verb-unification-02-compute
created: 2026-09-20
criteria:
  - id: ac-001
    summary: molpy no longer declares a Compute base of its own
    type: code
    pass_when: |
      src/molpy/compute/base.py does not exist, and no file under src/ or
      tests/ contains "compute.base" or "from .base import Compute".
    status: pending
  - id: ac-002
    summary: every molpy-owned compute shell exposes the verb compute()
    type: code
    pass_when: |
      No line matching "def __call__" exists under src/molpy/compute/, and the
      31 shell classes in the 17 modules cluster (3), decomposition (2),
      density (2), dielectric (2), diffraction (1), distribution (4),
      environment (1), hbond (1), msd (1), neighborlist (1), order (4),
      pmft (1), rdf (1), reorientation (1), shape (4), spatial (1),
      van_hove (1) each define "def compute(" with the same parameter list
      their previous __call__ had. No class under src/molpy/compute/ declares
      Compute as a base class and no "__call__ = compute" alias exists.
    status: pending
  - id: ac-003
    summary: the _config / dump() / **config_kwargs surface is gone
    type: code
    pass_when: |
      No occurrence of "super().__init__(" (any argument form — the shells use
      keyword and bare forms, never the literal `**config`), "_config",
      "def dump(" or "config_kwargs" remains under src/molpy/compute/; no file in
      the repo calls .dump() on a compute object; and the package docstring at
      src/molpy/compute/__init__.py:3-4 no longer says "call it".
    status: pending
  - id: ac-004
    summary: Compute is the molrs Protocol, reachable from molpy
    type: runtime
    pass_when: |
      `uv run --extra dev python -m pytest tests/test_compute/test_init.py`
      passes, asserting molpy.compute.Compute is molrs.compute.Compute, and
      "Compute" is present in src/molpy/compute/__init__.py __all__.
    status: pending
  - id: ac-005
    summary: voronoi shells are gone, replaced by identity re-exports
    type: code
    pass_when: |
      src/molpy/compute/voronoi.py declares no class and no def; it binds
      RadicalVoronoi and VoronoiIntegration directly to the corresponding
      molrs.compute.voronoi attributes, so
      molpy.compute.RadicalVoronoi is molrs.compute.voronoi.RadicalVoronoi
      and molpy.RadicalVoronoi is molpy.compute.RadicalVoronoi (the facade
      split at src/molpy/__init__.py:274,:276 is closed), and likewise for
      VoronoiIntegration.
    status: pending
  - id: ac-006
    summary: per-module unit tests green on the compute verb, numerics unchanged
    type: runtime
    pass_when: |
      `uv run --extra dev python -m pytest tests/test_compute/` passes; every
      shell exercised there is called as shell.compute(...); every numeric
      expectation carried over from the previous test bodies (RDF bins and the
      free-box ValueError, Nematic order > 0.9, PMFTXY (20, 20) counts shape,
      SDF orientation presence, DielectricResult.fit_debye tau ~ 5.0 rel 0.2,
      the three TypeError rejections) is byte-identical to its pre-change
      value, and no test file imports Compute except test_init.py.
    status: pending
  - id: ac-007
    summary: no test file spans two source modules; molrs re-export tests are gone
    type: code
    pass_when: |
      tests/test_compute/test_orientations.py, test_conductivity.py,
      test_onsager.py and test_persist.py no longer exist; test_order.py,
      test_spatial.py, test_pmft.py, test_result.py and test_init.py exist and
      no test_pmsd.py was added; each file under tests/test_compute/ imports
      symbols owned by exactly one src/molpy/compute/ module; the shared
      orientations-block builder lives in tests/test_compute/conftest.py as a
      fixture and no test file imports another test file.
    status: pending
  - id: ac-008
    summary: docs teach the owning verb and nothing points at compute.base
    type: docs
    pass_when: |
      No page under docs/ calls a molpy.compute object as obj(...) (the
      MaskPredicate examples in docs/tutorials/06_selector.md excepted);
      docs/compute/voronoi.md uses .build(...) / .integrate(...) and no longer
      carries the "you call the object itself; there is no .compute method"
      note; "::: molpy.compute.base" is gone from docs/api/compute.md and the
      JACF / PMSDCompute alias sentence and table are deleted while the
      raw-Compute → Fit → scale rule survives and the page states that
      IonicConductivity / DielectricSusceptibility currently exist, violate it,
      are routed to /mol:refactor and must not be used in new code; docs/developer/extending-compute.md
      describes subclass-free structural conformance with a compute() method
      and mentions neither super().__init__(**config) nor dump();
      docs/compute/index.md no longer says "configure, then call"; the Compute
      bullet in docs/getting-started/migration-0-14.md no longer calls shells
      "plain callables"; docs/developer/molrs-backend.md does not describe the
      operators as argument-forwarding shells without qualification.
    status: pending
  - id: ac-009
    summary: full lint + test gate green
    type: runtime
    pass_when: |
      `uv run --no-project --with 'tox>=4.23' --with ruff==0.16.1 --with
      ty==0.0.65 tox -e lint` and
      `uv run --extra dev python -m pytest tests/ -n auto` both exit 0.
    status: pending
  - id: ac-010
    summary: docs data scripts call the owning verb
    type: runtime
    pass_when: |
      The six modules under scripts/docs_data/ (structure, aggregate, order,
      dynamics, transport, angles) contain no `)(` call on a molpy.compute
      object — every former obj(...) site reads .compute(...) and
      dynamics.py:54 reads RadicalVoronoi().build(...) — and
      `PYTHONPATH=scripts python -c "import docs_data.structure, docs_data.aggregate, docs_data.order, docs_data.dynamics, docs_data.transport, docs_data.angles"`
      exits 0.
    status: pending
  - id: ac-011
    summary: verb table retires the compute debt row
    type: docs
    pass_when: |
      In .claude/notes/architecture.md § 设计铁律 4 the 分析 row's 状态 cell
      reads 已成立, the debt row for the compute shells / molpy Compute ABC is
      gone from 「本表尚未兑现的部分」, and a residual debt row names
      dielectric.py:155 from_dipole_series as the second analysis entry to be
      folded by the routed /mol:refactor.
    status: pending
out_of_scope:
  - "core/selector.py MaskPredicate — the chain's only __call__ carve-out"
  - "molrs-owned verbs (Onsager.correlation, Persist.pair_survival_tcf, LinearFit.fit, CumulativeTrapezoid, DebyeFit, RadicalVoronoi.build, VoronoiIntegration.integrate)"
  - "other families (builder 04, pack 03, io/typifier/wrapper unchanged)"
  - "routed debts: the two all-in-one dielectric recipe classes (+ from_dipole_series); the 26 _impl forwarding shells and the eight facade-split names; the 15 untested shell classes; tests/test_compute/test_dielectric.py coverage mismatch; CLAUDE.md:71 managed exemplar"
  - "compute API semantics: parameters, return types, kernels, units"
---

# Acceptance — api-verb-unification-02-compute

- **ac-001 / ac-002 / ac-003** 是静态审查项：`base.py` 消失、31 个 molpy 拥有的 shell 只剩 `compute` 一个动词、`_config`/`dump()`/`**config_kwargs` 三个面整体移除（最后一个是铁律 5 的静默水槽），且没有 `__call__ = compute` 之类的糖留下。
- **ac-004** 钉住铁律 6 的形状：契约由 molrs 拥有，用户从 `molpy.compute` 拿到的是同一个对象（identity，不是同名副本）。
- **ac-005** 钉住 voronoi 的裁定：不是改名，是把两个别名壳整个撤掉，避免后续重构对同一批公开名字做第二次改名。
- **ac-006** 是「数值不动」的唯一凭据：期望值逐字沿用，只有调用动词变了。任何期望值被改写都视为失败，不得为了让测试通过而松断言。
- **ac-007** 把逐模块镜像变成可核对的事实，同时钉住「纯 molrs re-export 不留 molpy 测试」「fixture 走 conftest、测试之间不 import」。
- **ac-008** 覆盖文档面：mkdocstrings 目标不悬空、别名表不留幽灵符号、扩展指南与迁移说明与新契约一致、voronoi 页讲 molrs 自己的动词。
- **ac-009** 是交付门。
- **ac-010** 覆盖 pytest 与 ty 都看不到的 `scripts/docs_data/`——文档曲线的生成脚本改名后不许静默炸掉。
- **ac-011** 让 sub-spec 01 的动词表在本 spec 落地时回到真实状态，并把 `from_dipole_series` 记为该行下的残余债。
