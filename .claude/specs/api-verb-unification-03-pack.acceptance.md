---
spec: api-verb-unification-03-pack
created: 2026-09-20
criteria:
  - id: ac-001
    summary: molpy.pack package tree fully removed from src and tests
    type: code
    pass_when: |
      Neither src/molpy/pack/ nor tests/test_pack/ exists; the six
      modules (pack/__init__.py, constraint.py, target.py,
      packer/__init__.py, packer/base.py, packer/packmol.py) and the
      three test files are gone, with no replacement module, alias, or
      deprecated shim anywhere under src/molpy/.
    status: pending
  - id: ac-002
    summary: facade no longer declares pack at any of its three sites
    type: code
    pass_when: |
      src/molpy/__init__.py contains no "pack" entry in the
      TYPE_CHECKING import block (was :29), in _LAZY_SUBMODULES
      (was :48), or in __all__ (was :317).
    status: pending
  - id: ac-003
    summary: facade __getattr__ rejects any unregistered attribute
    type: runtime
    pass_when: |
      uv run --extra dev python -m pytest tests/test_init.py passes: the
      single test asserts that a neutral name absent from
      _LAZY_SUBMODULES raises AttributeError via the __getattr__ at
      src/molpy/__init__.py:55-58. The test must not name "pack".
    status: pending
  - id: ac-004
    summary: no dangling molpy.pack reference remains in source
    type: code
    pass_when: |
      No file under src/molpy/ mentions molpy.pack or Packmol as a
      molpy symbol; specifically adapter/__init__.py:6-7 cites a
      molpy.wrapper exemplar instead of molpy.pack.packer.packmol, and
      builder/assembly/_replicas.py:3 points packing at molpack.
    status: pending
  - id: ac-005
    summary: layer-graph edges no longer name pack
    type: docs
    pass_when: |
      .claude/notes/architecture.md lines 70-91 (the custom annotation
      block) list no "pack" on any edge: it is gone from the io/,
      builder/ and wrapper/,adapter/ right-hand sides and the standalone
      "pack/ →" line is deleted; the CLAUDE.md import-direction sentence
      ("`io` may be imported by …", pre-chain line 224) no longer names pack. The /mol:map managed block
      (architecture.md:1-68) is left byte-identical.
    status: pending
  - id: ac-006
    summary: descriptive package tables drop the pack row
    type: docs
    pass_when: |
      CLAUDE.md has no pack row in Core Packages (was :199) and no
      test_pack/ in the tests tree (was :374); README.md has no pack
      table row (was :64) and its capability prose (was :28) no longer
      claims packing as a molpy package; docs/developer/
      architecture-overview.md has no pack row (was :16) and its
      wrapper row (was :20) no longer lists packmol. README.md:159 and
      the molpack mentions elsewhere are unchanged.
    status: pending
  - id: ac-007
    summary: user-facing docs route packing to molpack only
    type: docs
    pass_when: |
      docs/getting-started/external-tools.md:19 and
      docs/user-guide/15_mcp.md:616-625 name molpack only, with no
      molpy.pack path and no molpy *Constraint class; the Packing
      section of docs/getting-started/migration-0-14.md is either
      rewritten to record the removal (v0.14.0 still untagged) or
      deleted with the removal recorded on the next release page
      (v0.14.0 already tagged; recorded under the next-version heading of
      .claude/notes/release.md), and the delivery note states which
      branch was taken. docs/index.md, docs/user-guide/09_packing.md
      and docs/api/pack.md are unchanged.
    status: pending
  - id: ac-008
    summary: verb table retires the pack debt row
    type: docs
    pass_when: |
      In .claude/notes/architecture.md § 设计铁律 4 the 装填 row's 状态 cell
      reads 已成立 and the debt row for `Packer.__call__` + `molpy.pack` is
      gone from 「本表尚未兑现的部分」.
    status: pending
  - id: ac-009
    summary: managed blueprint rows are regenerated or recorded
    type: docs
    pass_when: |
      Either /mol:map has been run after landing (architecture.md managed
      block no longer lists pack at its former :24/:36/:60) or the delivery
      summary names those three stale rows and the /mol:map follow-up.
    status: pending
  - id: ac-010
    summary: lint and the full test suite pass after removal
    type: runtime
    pass_when: |
      uv run --no-project --with 'tox>=4.23' --with ruff==0.16.1 --with
      ty==0.0.65 tox -e lint passes, and uv run --extra dev python -m
      pytest tests/ -n auto reports zero failures and zero collection
      errors.
    status: pending
out_of_scope:
  - "molpack internals; no import dependency on molcrafts-molpack is added"
  - "compute/ (02) and builder/ (04) verb changes; the verb table itself (01)"
  - "core/region.py — geometry vocabulary stays where it is"
  - "the /mol:map managed block of architecture.md (:1-68, incl. :24/:36/:60)"
  - "docs/index.md, docs/user-guide/09_packing.md, docs/api/pack.md, README.md:159 — already molpack, unchanged"
  - "moving Packmol into wrapper/ (deletion removes the contradiction); deprecated shims; regressions/ directory"
---

# Acceptance — api-verb-unification-03-pack

- **ac-001 / ac-002** 是删除本身：包树消失、facade 三处声明消失。两者分开，因为漏掉 `__all__` 会让 facade 自相矛盾而 ac-001 仍然通过。
- **ac-003** 是本 spec 唯一的行为面，且**刻意不以 `pack` 命名**：测的是 `__getattr__` 的通用契约，不是某个已删名字的墓碑。
- **ac-004** 挡住"删了包但 docstring 还 `:mod:` 指过去"的悬空引用。
- **ac-005** 是最高优先级的一条：分层导入图是 `/mol:review --axis=arch` 的判据来源，而这段是 `/mol:map` **不会**重写的手工标注，漏改即长期说谎。同时明确要求**不要**去动托管块。
- **ac-006 / ac-007** 分别覆盖描述性表格与用户文档，并都写明哪些文件必须**保持不变**。
- **ac-008** 让动词表在本 spec 落地时回到真实状态。
- **ac-009** 让托管块的回填有人负责，而不是假设它会自愈。
- **ac-010** 覆盖删除的连带影响（`ty` 对 `TYPE_CHECKING` 块的解析、pytest 收集不到已删目录）。
