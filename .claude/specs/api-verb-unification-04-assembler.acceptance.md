---
spec: api-verb-unification-04-assembler
created: 2026-09-20
criteria:
  - id: ac-001
    summary: GraphAssembler exposes apply, with no assemble and no alias
    type: code
    pass_when: |
      src/molpy/builder/assembly/_assembler.py defines
      `def apply(self, world: Atomistic, selector: Selector) -> Atomistic`
      carrying the body and docstring contract of the former `assemble`
      (work = world.copy(); input never mutated), and no `def assemble`,
      `assemble = ` alias or deprecation shim exists anywhere under src/molpy/.
    status: pending
  - id: ac-002
    summary: No assemble identifier remains anywhere in the repo
    type: code
    pass_when: |
      A case-sensitive `\bassemble\b` sweep across src/, tests/, docs/,
      examples/, .claude/notes/ and CLAUDE.md returns zero matches referring to
      the method (calls, `:meth:` references, docstring wording such as
      "expand + assemble", table cells such as "Two assemble steps"); the noun
      forms assembler / assembles / assembled listed in Out of scope are the
      only permitted survivors.
    status: pending
  - id: ac-003
    summary: Renamed assembly test module passes
    type: runtime
    pass_when: |
      `uv run --extra dev python -m pytest tests/test_builder/test_assembly/`
      is green, with test_assembler.py calling `apply(world, selector)` throughout
      and the test formerly named test_assemble_returns_a_reacted_copy now named
      test_apply_returns_a_reacted_copy.
    status: pending
  - id: ac-004
    summary: Renamed call sites keep identical inputs and expectations
    type: runtime
    pass_when: |
      The 12 former `.assemble(` call sites (9 in
      tests/test_builder/test_assembly/test_assembler.py, 2 in
      tests/test_typifier/test_affected_region.py, 1 in
      tests/test_typifier/test_region_radii.py) differ from their pre-rename text
      only in the method name, and
      `uv run --extra dev python -m pytest tests/test_builder/test_assembly
      tests/test_typifier/test_affected_region.py tests/test_typifier/test_region_radii.py`
      is green.
    status: pending
  - id: ac-005
    summary: API-shape gates removed with behaviour cover retained
    type: code
    pass_when: |
      tests/test_builder/test_assembly/test_assembler.py and test_placer.py contain
      no `inspect.signature(` call and no `import inspect`, the three gate tests
      (test_selector_is_an_assemble_argument,
      test_constructor_has_no_charge_redistribution_switch,
      test_is_injected_into_the_assembler_constructor) are gone, and their three
      named behaviour covers
      (test_one_instance_accepts_different_pairing_rules,
      test_charged_leaving_group_cannot_silently_change_net_charge,
      test_spreads_overlapping_residue_templates) are present and passing.
    status: pending
  - id: ac-006
    summary: PolymerBuilder.build still overrides cleanly
    type: code
    pass_when: |
      src/molpy/builder/assembly/_polymer.py keeps `build(self, topology)` unchanged
      and calls `self.apply(world, TopologySelector(topology))`; `tox -e lint` reports
      no ty override / incompatible-signature error for src/molpy/builder/assembly/.
    status: pending
  - id: ac-007
    summary: Migration page carries the chain-wide verb table
    type: docs
    pass_when: |
      docs/getting-started/migration-0-14.md contains a "Verb table" section with
      rows for GraphAssembler.assemble -> apply, compute shell __call__ -> .compute,
      and "molpy.pack removed; packing is molpack" (no Packer.pack survivor is
      mentioned); the § Compute and § Packing prose owned by sub-specs 02 / 03 is
      unmodified by this spec.
    status: pending
  - id: ac-008
    summary: Docs and examples spell the new verb
    type: docs
    pass_when: |
      The 9 listed doc pages and 6 listed example scripts read `apply`, and the table
      row in docs/api/builder.md reads `apply(world, selector)`.
    status: pending
  - id: ac-009
    summary: Project knowledge and CLAUDE.md carry the landed name
    type: docs
    pass_when: |
      In .claude/notes/architecture.md § 设计铁律 4 the 图变换 row's 状态 cell reads
      已成立 and debt row 3 (`assemble`) is gone from 「本表尚未兑现的部分」 (or, if
      sub-spec 01 has not landed, the former line-155 example cites
      `GraphAssembler.apply`), and CLAUDE.md (§ `builder` module, free-form region outside the
      mol:bootstrap:managed markers) no longer names `GraphAssembler.assemble` — any
      one-line diff there having been surfaced to the operator, not landed silently.
    status: pending
  - id: ac-010
    summary: Full gate green
    type: runtime
    pass_when: |
      `uv run --extra dev python -m pytest tests/ -n auto` passes and
      `uv run --no-project --with 'tox>=4.23' --with ruff==0.16.1 --with ty==0.0.65 tox -e lint`
      (ruff format + ruff check + ty) is clean.
    status: pending
out_of_scope:
  - "PolymerBuilder.build and build_* shortcuts (recipe→structure family)"
  - "domain verbs (place, select, expand, grid, times) and class names (VirtualSiteBuilder, DrudeBuilder)"
  - "compute/ (02) and molpy.pack (03) code changes; this spec only adds the consolidated migration table"
  - "any deprecated alias for assemble; any hasattr / inspect.signature shape test"
  - "found-not-fixed: _polymer.py:93-94 double pass of topology; AssemblyFinalizer single-use subclass; remaining inspect gates in test_proximity.py / test_selector.py"
---

# Acceptance — api-verb-unification-04-assembler

- **ac-001 / ac-002** 一起定义"改名彻底"：新名字存在、旧名字在实现、调用点、文档、项目知识与 CLAUDE.md 中全部消失。"不留别名"由这两条 code 判据承担，**不**由测试承担。
- **ac-003 / ac-004** 是行为不变的证据：既有 12 个调用点除方法名外逐字不变仍然绿，说明这是纯改名。
- **ac-005** 钉住形状门禁的清除与行为替身的留存——删门禁不得连带削掉覆盖。
- **ac-006** 覆盖本裁定的技术理由（`build` 会撞 MRO），`ty` 必须干净。
- **ac-007 / ac-008 / ac-009** 覆盖用户可见面、项目知识与仓库配置的同步；ac-009 使本 sub-spec 无法在 CLAUDE.md 仍写旧名时收尾。
- **ac-010** 全量门禁。

> 本 spec 无 regression-example 判据：`.claude/notes/testing.md` § Test shape 与 CLAUDE.md § Tests 明令本仓库**没有** `regressions/` 目录。等价保障由 ac-003 / ac-004（既有单元用例逐字复用）承担。
