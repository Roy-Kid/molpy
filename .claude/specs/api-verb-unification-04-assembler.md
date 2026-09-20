---
slug: api-verb-unification-04-assembler
title: API 动词统一 (4/4)：GraphAssembler.assemble → apply，收尾迁移页动词表
status: approved
created: 2026-09-20
chain: api-verb-unification
depends_on:
  - api-verb-unification-01-verbtable
---

# API 动词统一 (4/4)：GraphAssembler.assemble → apply，收尾迁移页动词表

## Summary

把 `GraphAssembler.assemble(world, selector)` 更名为 `apply`，让"拿一张已有的图、返回一张改写过的新图"这件事在 molpy 里只有一个动词。此后 `StructureFinalizer.apply`、`VirtualSiteBuilder.apply`、molrs 的 `Reaction.apply` 与装配内核共用同一个名字，用户读一个方法签名就知道输入不被改动、返回值是新图。改名之外不动任何行为：签名、参数语义、警告、断言、电荷守恒检查、返回对象全部原样。实验阶段不留 `assemble` 兼容别名——保证方式是**不写这个别名**（铁律 3），不是加一个 API 形状测试。src / tests / docs / examples 的全部调用点同批改掉，并在 `docs/getting-started/migration-0-14.md` 补上覆盖整条链的动词对照表。

## Domain basis

N/A — 纯改名，无物理。

## Design

**裁定：`assemble` → `apply`，不是 `build`。** 两条理由，都已在仓库里核实：

1. **`build` 会在 MRO 上撞车。** `PolymerBuilder(GraphAssembler)`（`src/molpy/builder/assembly/_polymer.py:38`）已经有 `build(topology: CGSmilesGraphIR) -> Atomistic`（`:87`），它在 `:94` 调 `self.assemble(world, TopologySelector(topology))`。父类若叫 `build(world, selector)`，子类的 `build(topology)` 就是签名不兼容的覆写——违反 LSP，`ty` 会报出来。这不是风格问题，是类型错误。
2. **语义上它属于 `apply` 族。** `assemble` 接收一个**已经存在**的 world，`work = world.copy()`（`_assembler.py:121`）后改写并返回新图，入参全程不被修改（`_assembler.py:116-121`）。这是 graph→graph 的重写，与 `StructureFinalizer.apply(graph) -> Atomistic`（`src/molpy/builder/_finalize.py:43`）、`VirtualSiteBuilder.apply(struct) -> Atomistic`（`src/molpy/builder/virtualsite.py:69`，docstring 写的就是 "Return a new structure …; input untouched"）以及内核自己调用的 molrs `Reaction.apply_many_detailed`（`_assembler.py:153`）完全同形。`build` 表达的是 recipe→structure（`Lattice.build`、`PolymerBuilder.build`），那是另一族。类名（`VirtualSiteBuilder`、`DrudeBuilder`）与领域动词（`expand` / `place` / `grid` / `times`）不在这两族内，不改。

**已核实无新冲突。** `PolymerBuilder` 是 `GraphAssembler` 在 src 里唯一的子类；全 src 中 `def apply` 只出现在 `builder/virtualsite.py:69` 与 `builder/_finalize.py:43`（两个互不相干的类），`PolymerBuilder` 自身不定义 `apply`。改名不会引入任何覆写/MRO 冲突；`self._finalizer.apply(work)`（`_assembler.py:160`）与 `self.apply(...)` 是两个对象上的两个方法，不构成歧义。

**改动面：只有名字。** 签名 `(self, world: Atomistic, selector: Selector) -> Atomistic`、"``world`` is never mutated" 的契约、空选择的 `warnings.warn`、`_assert_disjoint` / `_assert_touched_covers_forming_bond` / `_assert_charge_conserved` 的行为与错误文案全部不动。新名字的 docstring 首行沿用 `apply` 族措辞（"Return a new …; ``world`` is never mutated."），使三处读起来是同一个契约。

**不留兼容别名，且不为此写测试。** 按 `.claude/notes/architecture.md` § 3：不加 `assemble = apply`、不加 deprecation shim。保证来自"不写"本身；`hasattr(GraphAssembler, "apply") and not hasattr(..., "assemble")` 这类断言是 API 形状门禁，被 `.claude/notes/testing.md:13-15` 与 CLAUDE.md § Tests 排除，**不写**。旧名字消失由 ac-001 / ac-002（`/mol:impl` 交付时核验的 code 判据）承担，12 个改名后的行为测试直接调用 `apply`，本身就是新名字可用的证据。

**顺带清除同类门禁（铁律：不留已见的腐）。** 本 spec 要重写的测试模块里现存三处 `inspect.signature(...).parameters` 形状门禁。它们的**正向**主张各有一个行为测试覆盖；前两条另含一个**反向**主张（「构造器上没有 `selector=` / `charges=` 这个旋钮」），那正是 CLAUDE.md § Tests 排除的 API 形状断言——**不再覆盖，直接放弃**，与 sub-spec 03 对 `"pack" not in __all__` 的处理同一口径。第三条被完整覆盖。因此删除门禁、保留行为测试：

| 门禁 | 已覆盖它的行为测试 |
|---|---|
| `tests/…/test_assembler.py:22-24` `test_selector_is_an_assemble_argument`（`selector` 在方法签名里、不在构造器里） | `:34-43` `test_one_instance_accepts_different_pairing_rules`——同一个 assembler 实例配两个不同 selector 各自成键，正是"selector 是每次调用的参数、不是构造期状态"的行为证明 |
| `tests/…/test_assembler.py:123-124` `test_constructor_has_no_charge_redistribution_switch`（构造器无 `charges=`） | `:97-105` `test_charged_leaving_group_cannot_silently_change_net_charge`——带电离去基团直接 `raise`，而不是被某个重分配开关吞掉 |
| `tests/…/test_placer.py:18-19` `test_is_injected_into_the_assembler_constructor`（构造器有 `placer=`） | `:23-28` `test_spreads_overlapping_residue_templates`——`placer=ResiduePlacer()` 真的把坐标推开了（最小间距 <1e-9 → >0.5） |

删除后两个文件的 `import inspect`（各自 `:3`）成为死导入，一并删掉（否则 ruff F401）。

**改名清单（精确盘点）。** `tests/` 内共 **13 处 `assemble` 标识符 + 2 个测试函数名**：
- 12 处 `.assemble(` 调用——`tests/test_builder/test_assembly/test_assembler.py:28,36,39,50,92,115,138,155,158`（9）、`tests/test_typifier/test_affected_region.py:366,382`（2）、`tests/test_typifier/test_region_radii.py:302`（1）；
- 1 处非调用引用——`test_assembler.py:23` 的 `inspect.signature(GraphAssembler.assemble)`（随上表整条删除）；
- 2 个测试函数名——`test_assembler.py:22` `test_selector_is_an_assemble_argument`（整条删除）、`:26` `test_assemble_returns_a_reacted_copy`（改名为 `test_apply_returns_a_reacted_copy`）。

`src/` 内：1 处真实调用（`_polymer.py:94`）、1 处 docstring 示例调用（`_replicas.py:28`）、5 处 docstring 交叉引用（`_assembler.py:5,205`，`_polymer.py:5,41,91`）。`docs/` 内 9 个页面（代码块 + 方法名正文），`examples/` 内 6 个脚本共 8 处调用。`CLAUDE.md:474` 另计（见下）。

**行号是链前快照。** 03 与 04 只依赖 01、彼此无序，且都改 `CLAUDE.md` 与 `_replicas.py`；若 03 先落地，`CLAUDE.md:474` 会因 `:199` 行被删而上移，`_replicas.py:28` 也可能漂移——一律按引文定位。

**CLAUDE.md:474 的路由（复核后更正）。** 该行位于 `<!-- mol:bootstrap:managed begin/end -->` 标记（`:25`–`:133`）**之外**的自由区（`:135` 起"Free-form additions below this line are preserved across re-runs"）。自由区被 `/mol:bootstrap` 原样保留、永不再生成——所以它不会自愈，必须手改。CLAUDE.md 是操作者拥有的配置文件：`/mol:impl` 应把这一行的 diff **呈交操作者确认**后再落，不得静默改写。（sub-spec 01 已把该行改写为指向动词表的指针；若 01 先落地，本 spec 落地时该行已无 `assemble` 字样，此项即为空操作。）

**知识同步（退役条款）。** sub-spec 01 写入 `.claude/notes/architecture.md` § 设计铁律 4 的动词表把「图变换」行标为「部分：`GraphAssembler` 现名 `assemble`，改名由 sub-spec 04 兑现」并附退役条款。本 spec 落地时把该行的状态格改为「已成立」、删除债务清单第 3 行（任务 8、ac-009）；若 01 尚未落地，则把第 155 行旧例句里的 `GraphAssembler.assemble` 改为 `GraphAssembler.apply`。

**迁移文档。** `docs/getting-started/migration-0-14.md` 目前没有动词小节。本 spec 新增一节动词对照表（该页通篇英文，标题用 "Verb table"，正文保持英文），三行覆盖整条链：`GraphAssembler.assemble` → `apply`；compute shell 的 `__call__` → `.compute`；`molpy.pack` 整包移除、打包改由 molpack 承担（sub-spec 03 删除整个包，**没有 `Packer.pack` 幸存者**）。行文以 sub-spec 01 的动词表为准；**02 负责重写 § Compute 正文，03 负责重写 § Packing 正文**，本 spec 只加这张合并表。

### Reuse decision

- `reuse` `VirtualSiteBuilder.apply`（`builder/virtualsite.py:69`）与 `StructureFinalizer.apply`（`builder/_finalize.py:43`）—— 它们就是 Closest pattern，改名后的方法沿用其命名与 docstring 契约措辞。
- `reuse` molrs `Reaction.apply` / `apply_many_detailed`（已在 `_assembler.py:153` 调用）—— 动词与下游内核对齐。
- `reuse` 既有行为测试 `test_one_instance_accepts_different_pairing_rules` / `test_charged_leaving_group_cannot_silently_change_net_charge` / `test_spreads_overlapping_residue_templates` 作为三条形状门禁的替代覆盖 —— 不新写测试。
- `new` —— 无。本 spec **不引入任何新公开符号**；公开面净变化是"一个方法换名"。
- `generalize` —— 无。

## Files to create or modify

源码（3）：

- `src/molpy/builder/assembly/_assembler.py` — `:116` 方法定义改名；模块 docstring `:5` 交叉引用；`:205` 分节注释
- `src/molpy/builder/assembly/_polymer.py` — `:94` 调用点；`:91` 的 `:meth:`assemble``；模块 docstring `:5` 与类 docstring `:41` 的 "expand + assemble" 措辞
- `src/molpy/builder/assembly/_replicas.py` — `:28` 类 docstring 示例

测试（4）：

- `tests/test_builder/test_assembly/test_assembler.py` — 删除 `:22-24`、`:123-124` 两条形状门禁与 `:3` 的 `import inspect`；9 处 `.assemble(` 改名；`:26` 测试函数名改为 `test_apply_returns_a_reacted_copy`
- `tests/test_builder/test_assembly/test_placer.py` — 删除 `:18-19` 形状门禁与 `:3` 的 `import inspect`
- `tests/test_typifier/test_affected_region.py` — `:366`、`:382`
- `tests/test_typifier/test_region_radii.py` — `:302`

文档（10）：

- `docs/api/builder.md` — `:12`、`:13` 表格行，`:127`、`:182` 代码块
- `docs/user-guide/02_assembly.md` — `:108`、`:305` 正文方法名，`:242`、`:276` 代码块
- `docs/user-guide/16_crosslinked_gel.md` — `:131`
- `docs/user-guide/topology/index.md` — `:25`（"Two assemble steps"，与 `10_dual_network.md:5` 同义）、`:48`、`:50`、`:130`、`:132`
- `docs/user-guide/topology/07_gel_exhaustive.md` — `:9`、`:21`
- `docs/user-guide/topology/08_gel_random.md` — `:16`
- `docs/user-guide/topology/09_end_linked.md` — `:18`
- `docs/user-guide/topology/10_dual_network.md` — `:5`、`:8`
- `docs/user-guide/topology/11_prepolymer_agent.md` — `:16`
- `docs/getting-started/migration-0-14.md` — 新增 "Verb table" 小节（3 行）

示例（6）：

- `examples/06_crosslinked_gel_gaff.py` — `:71`
- `examples/topology/07_gel_exhaustive.py` — `:19`
- `examples/topology/08_gel_random.py` — `:19`、`:33`
- `examples/topology/09_end_linked.py` — `:21`
- `examples/topology/10_dual_network.py` — `:1` docstring、`:27`、`:53`
- `examples/topology/11_prepolymer_agent.py` — `:32`

项目知识与配置（2）：

- `.claude/notes/architecture.md` — 动词表**图变换行**状态格改为「已成立」，删除债务清单第 3 行；若 01 尚未落地，则第 155 行的旧例句 `GraphAssembler.assemble` 由本 spec 改为 `GraphAssembler.apply`
- `CLAUDE.md` — `:474` § `builder` module 里的 `GraphAssembler.assemble` → `apply`（自由区、手改、呈交操作者确认；01 已落地则为空操作）

（无新建文件。`docs/user-guide/0{1,5,6}_*.ipynb` 是仅有的 notebook 源，其中无 `.assemble(` 调用，`05_polydisperse_systems.md:35,306` 命中的是普通名词，因此不需要 `scripts/render_notebooks.py` 重渲染。）

## Tasks

- [ ] Remove the three `inspect.signature(...).parameters` API-shape gates and their now-dead `import inspect`: `tests/test_builder/test_assembly/test_assembler.py` (`:22-24`, `:123-124`, `:3`) and `tests/test_builder/test_assembly/test_placer.py` (`:18-19`, `:3`); the behaviour tests named in the Design table already cover each claim
- [ ] Write failing tests: rename all 12 `.assemble(` call sites to `.apply(` in `tests/test_builder/test_assembly/test_assembler.py` (9), `tests/test_typifier/test_affected_region.py` (2) and `tests/test_typifier/test_region_radii.py` (1), and rename `test_assemble_returns_a_reacted_copy` to `test_apply_returns_a_reacted_copy`, keeping every input and assertion unchanged
- [ ] Rename `GraphAssembler.assemble` to `apply` in `src/molpy/builder/assembly/_assembler.py` (`:116` definition, `:5` module-docstring reference, `:205` section comment), keeping the signature, docstring contract and all assertions unchanged and adding no alias
- [ ] Update the in-tree callers and docstring references: `src/molpy/builder/assembly/_polymer.py` (`:94` → `self.apply(...)`, `:91` `:meth:` reference, `:5` and `:41` wording) and the example in `src/molpy/builder/assembly/_replicas.py:28`
- [ ] Update the doc pages to `apply`: `docs/api/builder.md`, `docs/user-guide/02_assembly.md`, `docs/user-guide/16_crosslinked_gel.md`, `docs/user-guide/topology/index.md`, `docs/user-guide/topology/07_gel_exhaustive.md`, `docs/user-guide/topology/08_gel_random.md`, `docs/user-guide/topology/09_end_linked.md`, `docs/user-guide/topology/10_dual_network.md`, `docs/user-guide/topology/11_prepolymer_agent.md`
- [ ] Update the example scripts to `apply`: `examples/06_crosslinked_gel_gaff.py`, `examples/topology/07_gel_exhaustive.py`, `examples/topology/08_gel_random.py`, `examples/topology/09_end_linked.py`, `examples/topology/10_dual_network.py`, `examples/topology/11_prepolymer_agent.py`
- [ ] Add the "Verb table" section to `docs/getting-started/migration-0-14.md` with three rows (`GraphAssembler.assemble` → `apply`; compute shell `__call__` → `.compute`; `molpy.pack` removed, packing is molpack), taking wording from sub-spec 01's verb table; do not touch the § Compute or § Packing prose (owned by sub-specs 02 / 03)
- [ ] Retire the debt in the verb table: flip the 图变换 row's 状态 cell in `.claude/notes/architecture.md` § 设计铁律 4 to 已成立 and delete debt row 3 (`assemble`); if sub-spec 01 has not landed, rewrite the former line-155 example to `GraphAssembler.apply` instead; and — if it still names `assemble` — change `CLAUDE.md:474` to point at the verb table (CLAUDE.md is operator-owned: surface that one-line diff for approval instead of landing it silently)
- [ ] Run full check + test suite

## Testing strategy

遵循 `.claude/notes/testing.md` § Test shape：仅单元测试，路径镜像 `src/`，一个测试一个行为，无 e2e、无 golden、无 `regressions/`、**无 API 形状 / source-text 门禁**。本 spec 不改数值，不新增数值验证用例，且**净减少**一个测试文件里的测试数量（删掉三条形状门禁）。

- **既有用例改名即回归**：`tests/test_builder/test_assembly/test_assembler.py`（`TestGraphAssembler`，镜像 `src/molpy/builder/assembly/_assembler.py`）的 9 个调用点换成 `.apply(`，**入参、fixture、断言逐字不变**。它们已覆盖 happy path（反应后返回新图）、`.copy()` 隔离（`:32` 断言入参 world 无键）、空选择告警、重叠 binding 拒绝、未知 map number 拒绝、电荷守恒、typifier/reach 校验、finalize 各档、frontier 结点分型、残基标签无关性。改名若破坏行为，这批测试即红。
- **跨模块调用点**：`tests/test_typifier/test_affected_region.py`、`tests/test_typifier/test_region_radii.py` 测的是 region 行为，本次只换方法名；保持绿即证明装配内核对下游 typifier 的契约未变。
- **不新增测试**：不写 `hasattr` / `inspect.signature` 之类公开面断言（形状门禁，项目法排除）。"无兼容别名"由**不写别名**保证，交付时由 ac-001 / ac-002 核验；"新名字可用"由 12 个改名后的行为测试证明。
- **删除的三条门禁各有行为替身**（见 Design 表），删除不降低覆盖：`test_one_instance_accepts_different_pairing_rules`、`test_charged_leaving_group_cannot_silently_change_net_charge`、`test_spreads_overlapping_residue_templates` 全部保留且必须绿。
- **不写的其它测试**：不测 docs/examples 文本；不为 `PolymerBuilder.build` 新增用例（行为未变，既有 polymer 测试已覆盖）；不建 `regressions/` 目录（CLAUDE.md § Tests 与 `testing.md` 明令无此目录）。
- **命令**：单档 `uv run --extra dev python -m pytest tests/test_builder/test_assembly/`；全量门禁 `uv run --extra dev python -m pytest tests/ -n auto` + `tox -e lint`（含 `ty`，覆写兼容性由它把关）。

## Out of scope

- `PolymerBuilder.build(topology)` 及全部 `build_*` 快捷方法：保留原名原签名（recipe→structure，另一族）。
- `Placer.place`、`Selector.select`、`MonomerLibrary.expand`、`Replicas.grid` / `times` 等领域动词，以及 `VirtualSiteBuilder` / `DrudeBuilder` 等**类名**：均不在两大动词族内，不改（sub-spec 01 已作此记录）。
- `compute/`（shell `__call__`→`compute`，sub-spec 02）与 `molpy.pack`（整包移除，sub-spec 03）的代码改动：本 spec 只在迁移页把结果汇成一张表。
- 不提供 `assemble` 的 deprecated 别名或过渡期（`architecture.md` § 3）。
- `docs/**` 里作为普通名词出现的 "assembler / assembles / assembled"（`05_polydisperse_systems.md:35,306`、`16_crosslinked_gel.md:106,138`、`docs/compute/*.md` 等）不改；只改方法名引用。

**已发现、本 spec 不修、已路由（铁律：命名而非沉默）：**

1. **同一份数据穿两个对象** —— `src/molpy/builder/assembly/_polymer.py:93-94`：`self._library.expand(topology)` 之后又 `TopologySelector(topology)`，`topology` 被喂给两个对象。这正是 `.claude/notes/architecture.md:174` 点名的"边界切错"信号。路由：**单独裁决一次**（`MonomerLibrary.expand` 连同配对规则一起产出，或 `TopologySelector` 从展开后的 world 推导），不在一次改名里顺手重切边界。
2. **单次使用的兼容性子类** —— `src/molpy/builder/assembly/_finalize.py:10-14` `AssemblyFinalizer(StructureFinalizer)` 类体只有一个钉死的默认值 `perceive_aromaticity=True`，却导出于 `builder/__init__.py:41,103`、`builder/assembly/__init__.py:14,39`，写进 `docs/api/builder.md:16`，且只被测了那个默认值——CLAUDE.md § Forbid 的"二次封装门面"。路由：**`/mol:refactor`**（在 `_assembler.py:99` 直接构造 `StructureFinalizer(Finalization(finalize), bonded, perceive_aromaticity=True)`，删模块与四处导出）。
3. **同类形状门禁的剩余两处** —— `tests/test_builder/test_assembly/test_proximity.py:114`、`:124` 同样是 `inspect.signature(...).parameters` 门禁（`ExplicitPairSelector` 无 `spacing=`、`SpacingSelector` 无 `pairs=`），模块不在本次改名面内；另有 `tests/test_builder/test_assembly/test_selector.py:19` 的 `inspect.isabstract(Selector)`（较轻，抽象性可用 `pytest.raises(TypeError)` 行为化，参见 `test_placer.py:14-16`）。路由：**`/mol:test`**，与本 spec 同一裁定口径处理。
