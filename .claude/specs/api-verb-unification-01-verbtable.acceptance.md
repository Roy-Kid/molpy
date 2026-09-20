---
spec: api-verb-unification-01-verbtable
created: 2026-09-20
criteria:
  - id: ac-001
    summary: Family-to-verb table present in 设计铁律 4 with all 11 rows
    type: docs
    pass_when: |
      `.claude/notes/architecture.md` § 设计铁律 › 4 contains one markdown
      table whose rows are exactly 构造 / 图变换 / 分析 / 装填 / 分型 /
      3D 生成 / 外部进程 / 发射 / 读写 / 选择 / 放置 — 11 rows,
      no more, no fewer.
    status: pending
  - id: ac-002
    summary: Table carries all six columns with cited symbols and paths
    type: docs
    pass_when: |
      Every row has non-empty 变换族 / 输入 → 输出 / 动词 / 成员 / 状态 / 依据
      cells; the 动词 column reads build / apply / compute / pack / typify /
      generate / run / emit / read `/` write / select / place in that row
      order; and the 成员 cells cite `builder/crystal.py:210`,
      `builder/_finalize.py:43`, `builder/virtualsite.py:69`,
      `typifier/base.py:99`, `conformer/__init__.py:42`,
      `core/selector.py:24` and `docs/api/pack.md`.
    status: pending
  - id: ac-003
    summary: Unlanded rows name their carrier sub-spec by exact slug
    type: docs
    pass_when: |
      The 分析 row's 状态 cell names `api-verb-unification-02-compute` and
      the words 已声明的债; the 装填 row's names
      `api-verb-unification-03-pack` and 已声明的债; the 图变换 row's names
      `api-verb-unification-04-assembler` and flags that `GraphAssembler`
      is still called `assemble` at `builder/assembly/_assembler.py:116`.
      No row claims 已成立 for a verb that does not exist in the tree today,
      and each of those three 状态 cells carries the retirement clause
      obliging the named sub-spec to flip the cell to 已成立 and delete its
      debt row when it lands (sub-specs 02 / 03 / 04 carry that as a task).
    status: pending
  - id: ac-004
    summary: Compute ownership ruling recorded and tied to the open decision
    type: docs
    pass_when: |
      A paragraph after the table states that the one `Compute` is
      `molrs.compute.protocol.Compute`, reached as `molpy.compute.Compute`
      by identity re-export, that `src/molpy/compute/base.py` (the ABC at
      `:18` whose abstract method at `:51` is `__call__`) is deleted by
      sub-spec 02 with no `__call__` sugar left, names the loss of `dump()`
      and `**config` (`compute/base.py:41-48`, `:65-71`) as sub-spec 02's
      owned breaking change rather than a lossless re-export, is time-scoped
      ("after 02 lands"), and cites `.claude/specs/INDEX.md` as the place
      recording this as the open decision that this chain closes.
    status: pending
  - id: ac-005
    summary: No frozen subclass count in the constitution
    type: docs
    pass_when: |
      Neither the table nor any footnote states "33" or "34" as *the* count
      of compute shells; the 分析 row's 成员 cell is worded structurally
      ("every analysis class implementing `compute()`; subclasses of the ABC
      only until 02 lands"), and footnote 5 gives the verb-based recipe
      `rg -n 'def compute\(' src/molpy/compute` (0 until 02 lands) with the
      pre-02 base-class recipe and its doctest caveat (`compute/base.py:27`;
      34 lines on 2026-09-20 meant 33 real subclasses) as the interim.
    status: pending
  - id: ac-006
    summary: pack deletion set stated; carve-out claim is time-scoped
    type: docs
    pass_when: |
      Footnote 4 (or the 装填 row) lists sub-spec 03's deletion set as the
      whole `src/molpy/pack/` package (`Packmol`, `Packer`, `Target`, all
      `*Constraint`), `tests/test_pack/`, and the facade entries at
      `src/molpy/__init__.py:29`, `:48`, `:317` and the two docstring
      cross-references `builder/assembly/_replicas.py:3` and
      `adapter/__init__.py:6-7`; and the "sole carve-out"
      sentence is explicitly conditioned on the chain having landed, with
      `Packer.__call__` (`pack/packer/base.py:69`) named as sub-spec 03's
      declared debt until then.
    status: pending
  - id: ac-007
    summary: Six footnotes present, including the build_* census and family boundary
    type: docs
    pass_when: |
      Six numbered footnotes state: (1) entry verb vs descriptive names,
      listing `virtualsite.py:82/126/205 build_sites`,
      `engine/openmm.py:303 generate_inputs`,
      `builder/polymer/sequences.py generate_sequence` and
      `adapter/rdkit.py:105 generate_3d` as outside the families;
      (2) the `build_*` census splitting shortcuts kept
      (`_polymer.py:96/100/110/114`, `builder/ambertools.py:241`,
      `builder/polymer/system.py:146`) from debt
      (`parser/moltemplate/builder.py:551`, `:967`);
      (3) `DistributionIR.build()` (`builder/polymer/distributions.py:28`)
      as a factory-style constructor deferred to a separate `/mol:refactor`;
      (4) the `__call__` policy; (5) the no-frozen-count grep;
      (6) class names out of scope and `MonomerLibrary.expand`
      (`builder/assembly/_library.py:57`), `Replicas.grid` / `.times`
      (`_replicas.py:40/89`) and `Placer.place` as domain verbs outside the
      two families. The 3D 生成 row's 成员 cell lists `Conformer.generate` only.
    status: pending
  - id: ac-008
    summary: Family boundary is input-to-output; class names out of scope
    type: docs
    pass_when: |
      The text immediately above the table states that family membership is
      decided by 输入 → 输出 and not by class name or suffix, and that the
      table renames no class; footnote 6 names both boundary cases —
      `VirtualSiteBuilder` (`builder/virtualsite.py:61`), `DrudeBuilder`
      (`:92`) and `Tip4pBuilder` (`:173`) as 图变换 despite the `Builder`
      suffix, versus `GrapheneBuilder` (`nanostructure/graphene.py:15`) and
      `CarbonTubeBuilder` (`carbon_tube.py:15`) as 构造 — and rules
      `MonomerLibrary.expand` (`builder/assembly/_library.py:57`),
      `Replicas.grid` (`_replicas.py:40`), `Replicas.times` (`:89`) and
      `Placer.place` (`_placer.py:46`) / `ResiduePlacer.place` (`:75`) domain
      verbs outside both families, keeping their current names; the 选择 row
      marks `Atomistic.select` / `CoarseGrain.select` as descriptive (naming
      owned by § Graph sink decisions, not authorised for rename here); and states the
      criterion that the core data-model API of `Atomistic`/`CoarseGrain`/`Frame`
      (`def_*`, `del_*`, `copy`, `merge`, `move`, `rotate`, `scale`, `align`,
      `replicate`, `extract_subgraph`, `to_frame`) is outside both families,
      governed by § Graph sink decisions and CLAUDE.md § What must never change
      casually. The 放置 row's 依据 cell reads 领域动词 (not 非公开面) and its
      成员 cell names `ResiduePlacer.place` and the `builder/assembly/__init__.py`
      / `builder/__init__.py` exports.
    status: pending
  - id: ac-009
    summary: Declared-debt table lists all nine debts with owners
    type: docs
    pass_when: |
      A 「本表尚未兑现的部分」 table has nine rows naming, with path:line and
      owner: compute shells → sub-spec 02; `Packer.__call__`
      (`pack/packer/base.py:69`) + `molpy.pack` → sub-spec 03; `assemble`
      (`builder/assembly/_assembler.py:116`) → sub-spec 04; the `Selector`
      collision (`core/selector.py:43` vs `builder/assembly/_selector.py:26`)
      → separate `/mol:refactor`; the free `emit` dispatcher
      (`io/emit/__init__.py:51`) + `emit_python`
      (`parser/moltemplate/py_emitter.py:42`) → separate `/mol:refactor`;
      the free `build_*` factories (`parser/moltemplate/builder.py:551`,
      `:967`) → separate `/mol:refactor`; the § 4 boundary signal still live
      at `builder/assembly/_polymer.py:93-94` → separate decision;
      `CLAUDE.md:71` (managed, two false examples) → operator
      `/mol:bootstrap`; the managed Style summary (`architecture.md:43-49`)
      → `/mol:map`. The table is followed by the sentence that it is
      bookkeeping and routing, not a list of permitted exceptions.
    status: pending
  - id: ac-010
    summary: select/mask, one-verb-per-family and map expectation recorded
    type: docs
    pass_when: |
      The 选择 row states `select` is the entry verb and
      `MaskPredicate.mask(block) -> ndarray` (`core/selector.py:24`) is the
      producer hook, not a second entry verb; a clause states a class may
      carry one verb per family it participates in, with `PolymerBuilder`
      (`build` + inherited `apply`) as the example; and one sentence asks
      the next `/mol:map` to stop restating verb rulings in the managed
      Style summary (`architecture.md:43-49`) and point its verb column at
      the table.
    status: pending
  - id: ac-011
    summary: assemble→apply ruling recorded with its LSP collision reason
    type: docs
    pass_when: |
      The same section states `GraphAssembler.assemble` becomes `apply`
      (not `build`) and names both reasons: `PolymerBuilder` already has
      `build(topology)` at `builder/assembly/_polymer.py:87` so a parent
      `build(world, selector)` would collide (LSP / `ty`), and `assemble`
      is a graph → graph transform.
    status: pending
  - id: ac-012
    summary: Old two-example line replaced; section header count fixed
    type: docs
    pass_when: |
      The line ``同族变换用**同一个动词**:`GraphAssembler.assemble`、
      `VirtualSiteBuilder.apply`。`` no longer appears in
      `.claude/notes/architecture.md`, the rule sentence
      「同族变换用**同一个动词**」 still opens the new table block, and the
      § 设计铁律 header (pre-edit line 113) reads 「六条硬约束」, matching
      its six numbered subsections.
    status: pending
  - id: ac-013
    summary: CLAUDE.md points at the table and restates no verb
    type: docs
    pass_when: |
      `CLAUDE.md` contains a pointer line in § Architecture Overview whose
      line number exceeds that of `<!-- mol:bootstrap:managed end -->` and
      which names no verb; the former bullet at lines 473–474 no longer
      contains `Lattice.build`, `PolymerBuilder.build` or
      `GraphAssembler.assemble` nor the claim "there are no free `build_*`
      / `create_*` factories", and instead points at
      `.claude/notes/architecture.md` § 设计铁律 4; and
      `grep -cE '(Lattice|PolymerBuilder|GraphAssembler)\.(build|assemble)' CLAUDE.md`
      returns 0 outside the managed block (lines 25–133 are untouched).
    status: pending
  - id: ac-014
    summary: Cross-repo note points at the table and scopes molrs rows
    type: docs
    pass_when: |
      `.claude/notes/cross-repo-spec-map.md` gains one line pointing at
      `.claude/notes/architecture.md` § 设计铁律 4 and stating that rows /
      members marked `molrs:` describe molrs and do not bind it (molpy
      adopts molrs verbs; the reverse does not hold — sink direction).
    status: pending
  - id: ac-015
    summary: Diff confined to three files; managed blocks byte-identical
    type: docs
    pass_when: |
      `git diff --stat` lists only `.claude/notes/architecture.md`,
      `CLAUDE.md` and `.claude/notes/cross-repo-spec-map.md`; every
      `architecture.md` hunk falls on pre-edit line 113 or inside
      § 设计铁律 › 4 (pre-edit lines 149–175; line 176 opens § 5); every `CLAUDE.md` hunk falls
      after pre-edit line 133; and neither
      `<!-- mol:map:managed begin…end -->` nor
      `<!-- mol:bootstrap:managed begin…end -->` appears in any hunk.
    status: pending
  - id: ac-016
    summary: Gate stays green and no code leaked into the diff
    type: runtime
    pass_when: |
      `uv run --no-project --with 'tox>=4.23' --with ruff==0.16.1 --with
      ty==0.0.65 tox -e lint` exits 0 and `uv run --extra dev python -m
      pytest tests/ -n auto` exits 0 with the same pass count as before the
      edit, and `git diff --name-only` matches no path under `src/` or
      `tests/`.
    status: pending
out_of_scope:
  - "src/ and tests/ — renames and deletions land in sub-specs 02–04"
  - "the nine declared debts themselves (Selector collision, free emit dispatcher, free build_* factories → /mol:refactor; the _polymer.py:93-94 boundary signal → separate decision; CLAUDE.md:71 → /mol:bootstrap; managed Style summary → /mol:map)"
  - "the managed block of architecture.md (refresh via /mol:map after the chain)"
  - "CLAUDE.md managed block (lines 25–133); the Core Packages `pack` row goes with sub-spec 03"
  - "docs/ — Packmol doc cleanup belongs to sub-spec 03"
  - "regressions/ directory, docs-block gates, unit tests (none apply to a docs-only change)"
---

# Acceptance — api-verb-unification-01-verbtable

ac-001 – ac-002 是表的骨架（行齐、列齐、路径真）。ac-003 是本次修订的**核心判据**：宪法不得写一句当下为假的话——三个未落地的族必须逐行点名承载 sub-spec 的**精确 slug**，任何一行对不存在的动词声称「已成立」即判失败。ac-004 关掉 `.claude/specs/INDEX.md` 里悬着的 `Compute` 二义性；ac-005 保证宪法里不再冻结一个会腐烂的计数（并记下 grep 会多数一行 doctest 的坑）；ac-006 让「唯一豁免」这句话有时间限定，否则它在 sub-spec 03 落地前就是假的。ac-008 是**族边界判据**，专为下一次改名而设：`VirtualSiteBuilder` 顶着 `Builder` 后缀却是图变换族、`GrapheneBuilder` 是构造族——同一个后缀落在两族里，正是「后缀不是判据」的证明；同时把 `MonomerLibrary.expand` / `Replicas.grid` / `Replicas.times` / `Placer.place` 裁在两族之外。ac-007、ac-009 – ac-011 覆盖脚注、债务清单、select/mask 分工、一族一动词条款与 `assemble → apply` 的理由。ac-012 兼管旧反例的删除与标题「四条 → 六条」。ac-013 要求 CLAUDE.md 两处都变成指针且一个动词名都不留——`grep -cE` 返回 0 是最硬的判法。ac-014 把表的单向性写进跨仓文档。

ac-015 是边界守卫：本 spec 最大的漂移风险是顺手把陈旧的 managed 块也刷了，判据用 pre-edit 行号把 hunk 钉死在第 113 行与 § 4（149–175）。ac-016 是唯一的 runtime 判据：markdown 改动本不入 ruff/ty 视野，门禁全绿加上 `git diff --name-only` 无 `src/` / `tests/` 命中，等于证明没有任何代码混进这次纯文档提交。

按 `.claude/notes/testing.md` § Test shape，本 spec 不产生单测、不建 `regressions/`、不加 docs-block gate；上面没有 `type: code` / `type: scientific` 判据、也没有 regression-example 判据，是这条项目规则的直接后果，不是遗漏。
