---
slug: api-verb-unification-03-pack
title: API 动词统一 (3/4)：移除 molpy.pack — 打包能力归属 molpack
status: approved
created: 2026-09-20
chain: api-verb-unification
depends_on:
  - api-verb-unification-01-verbtable
---

# API 动词统一 (3/4)：移除 molpy.pack — 打包能力归属 molpack

## Summary

`molpy.pack` 是一个零消费者的死包：`Packmol`（`src/molpy/pack/packer/packmol.py`，约 625 行 packmol 二进制的 subprocess 包装）在 `src/molpy` 内部除自身 `__init__` 导出外无人引用，`docs/`、`examples/` 与全部兄弟仓均无调用方，而仓内文档早已把打包判给外部的 `molcrafts-molpack`。本 spec 不为它的 `pack` 动词改名，而是整包删除 `src/molpy/pack/`、facade 上的惰性属性与 `tests/test_pack/`，并同步清理删除后会说谎的全部架构描述：分层导入图的四条边、四张包表、两处源码 docstring 与三处用户文档。交付后 `import molpy` 不再暴露 `pack`，`api-verb-unification` 动词表中的 `pack` 一格指向 `molpack.pack`，仓库里不留任何 deprecated shim。

## Domain basis

N/A — 删除一个 subprocess 包装与一个 numpy 罚函数模块，不引入任何方程、单位或验证目标。

## Design

**判定依据（OQ-2）。** 四条证据支持删除而非改名：

1. **零消费者。** `Packmol` 仅被自身包导出（`pack/__init__.py:11`、`packer/__init__.py:2`）。复核 `src/`、`docs/`、`tests/`、`pyproject.toml`、`zensical.toml` 与 `/home/jicli594/work/molcrafts/` 下全部兄弟仓后无新增消费者；仅有的两处"命中"——`/home/jicli594/work/molcrafts/.worktrees/molpy-0.13.1/` 与 `/home/jicli594/work/molcrafts/workspaces/lab-new/snapshots/20260920T132404Z/molpy/`——都是 molpy 自己的历史快照，不是消费方。
2. **能力在生态内已有真正的 owner，且仓内文档已经这么写了。** `docs/api/pack.md:3-4` 开篇即"packing … via **molpack** (`molcrafts-molpack`)"，`docs/developer/architecture-overview.md:16` 把 `pack` 标成 "Legacy packing helpers; new packing is **molpack**"，`README.md:159` 把 molpack 列为生态成员。也就是说 molpy 侧留着的是一个已被自家文档宣告作废的重复实现。CLAUDE.md § Forbid 的 "All-in-one façade APIs" 正对着 `Packmol.__call__`（一次跑完 写输入→起进程→读回 Frame），而 CLAUDE.md § Pattern: ForceField I/O 已经立下先例——"a reader that only forwarded to another was deleted, not kept"。（注意：CLAUDE.md:216-221 的 **Sink direction** 判据只管 molpy↔molrs，molpack 不在其射程内，故本判定不引它。）
3. **该基类本身是坏的。** `Packer`（`packer/base.py:49-81`）声明抽象 `pack`，`__call__` 委托给 `pack`；唯一子类 `Packmol` 反过来实现 `__call__` 并让 `pack` 回调 `__call__`（`packmol.py:596-607`）——倒置委托，若子类不覆写 `__call__` 就无限递归。
4. **它还在静默吞参数。** `Packmol.__call__`（`packmol.py:56-66`）接受 `workdir` / `tolerance` / `cleanup` / `pbc`，而 `Packmol.pack`（`:596-601`）只接 `targets` / `max_steps` / `seed` 并以 `return self(targets, max_steps=..., seed=...)` 丢掉其余四个；`packmol.py:81` 的 docstring 更是自认 `**kwargs: Additional arguments (ignored for now)`，`Packer.__call__`（`base.py:69-81`）同样收 `**kwargs` 却从不转发。stage `experimental` 下把这样一个无人使用、且违反 CLAUDE.md § Iron law 的形状改名保留，等于给动词表新增一处已知 rot。

**被孤立的两个模块一并删除。** `pack/constraint.py`（纯 numpy penalty 几何）与 `pack/target.py` 的消费者只有 `packmol.py:20-21` 与它们各自的测试。它们的**几何词汇并未丢失**：`src/molpy/core/region.py` 已提供 molrs 支撑的同一套可布尔组合几何（`BoxRegion` / `SphereRegion` / `AndRegion` / `OrRegion` / `NotRegion`，基于 `molrs.Cuboid` / `molrs.Sphere`）。更直接的证据是 `constraint.py:5` **本来就 `from molpy.core.region import BoxRegion, SphereRegion`**，并在 `:75` / `:105` / `:144` / `:174` 把它赋给 `self.region` —— 而 `self.region` 全模块 **0 次读取**，是纯死字段。该模块相对 `core.region` 多出来的只有 penalty / dpenalty 打分，而打分是打包优化器的内部事务，属于 molpack。`docs/api/pack.md` 记录的 `Target` / restraint 也是 molpack 的类型，不是这两个模块。

**Facade 收口。** `src/molpy/__init__.py` 有三处 `pack`，全部删除：`:29`（`TYPE_CHECKING` 导入块）、`:48`（`_LAZY_SUBMODULES`）、`:317`（`__all__`）。删除后 `mp.pack` 由 `__getattr__`（`__init__.py:55-58`）抛 `AttributeError`。

**分层导入图必须同改（删除后不会自愈）。** 这是本 spec 最容易漏、后果最重的一处——图里还写着 `pack` 的边，会让后续 `/mol:review --axis=arch` 拿一个不存在的包去判边：

- `.claude/notes/architecture.md:80`（`io/ → builder, pack, engine, typifier`）、`:83`（`builder/ → pack`）、`:84`（`pack/ → (application code only)` 整行删除）、`:86`（`wrapper/, adapter/ → builder, pack, engine`）。**关键**：这段（`:70-91`）标注了 "custom annotation — preserved across /mol:map runs"，`/mol:map` **不会**重写它，必须手改。
- `CLAUDE.md:224` 的 import-direction 段落散文："`io` may be imported by `builder`, `pack`, `engine` and `typifier`" —— 去掉 `pack`（按引文定位，不按行号）。

**同文件中 `/mol:map` 托管、本 spec 不手改的三行。** `.claude/notes/architecture.md:24`（包清单 `- **pack** — Packmol packer …`）、`:36`（`- **pack**: Packmol, Target, constraint types`）、`:60`（Layer roles 表 `| pack | Spatial placement |`）都在 `<!-- mol:map:managed -->` 块内（该块结束于 `:68`），由 `/mol:map` 重新生成。**不要手工编辑这三行**——下次 `/mol:map` 运行时自愈；手改只会与生成器打架。

**其余描述性表格（非托管，必须手改）。**

- `CLAUDE.md:199` —— Core Packages 表行 `| `pack` | Packing: Packmol, constraints, density targets |`，删行。
- `CLAUDE.md:374` —— 测试树那行里列着 `test_pack/`，删掉该条目。
- `README.md:64` —— `| **`pack`** | Packmol-based packing with density targets |`，删行；`README.md:28` 的散文能力清单里去掉 `packing`（或改述为经 molpack）。`README.md:159` 已列 molpack，保留不动。
- `docs/developer/architecture-overview.md:16` —— `| `pack` | Legacy packing helpers; …|` 删行；`:20` —— `| `wrapper` | Subprocess boundaries to external CLI tools (antechamber, packmol, …) |` 去掉 `packmol`：这句**本来就是假的**（packmol 包装从不在 `wrapper/` 里，见下），删包之后更假。

**删除后会悬空的源码与用户文档引用。**

- `src/molpy/adapter/__init__.py:6-7` —— docstring 把 `:mod:`molpy.pack.packer.packmol`` 当作"wrapper 侧的那一个范例"。这句同时暴露了一个旧错：一个约 600 行的 subprocess 包装住在 `pack/` 而非 `wrapper/`，与 CLAUDE.md:201（`wrapper` = External CLIs）矛盾。改写为指向真正的 wrapper 范例（`molpy.wrapper` 的 AmberTools 系列），adapter/wrapper 两种桥接模式的叙述保留。
- `src/molpy/builder/assembly/_replicas.py:3` —— "Packing production boxes belongs in :mod:`molpy.pack`"，改指 molpack。
- `docs/getting-started/external-tools.md:19` —— 默认路径表 `| Pack a box | `molpy.pack` → molpack |`，改为只写 molpack。
- `docs/user-guide/15_mcp.md:616-625` —— `molmcp_outline(path="molpy/pack")` 的示范输出，列表本身还是混的（`Molpack` / `InsideBoxRestraint` 属 molpack，`OutsideBoxConstraint` / `InsideSphereConstraint` / `MinDistanceConstraint` / `Target` 属被删的 `molpy.pack`）。路径与输出改写为 `molpack`。
- `docs/getting-started/migration-0-14.md:102-108` —— "## Packing" 一节描述的正是这批 `*Constraint` 的 0.14 行为变更。**前置条件**：`.claude/notes/release.md:21` 记 v0.14.0 为 `untagged`，故本 spec 落地时若 0.14.0 **仍未打 tag**，就地把该节改写为"整包移除、打包走 molpack"，与同文件 `:117` 的 "Removed without replacement" 口径一致；若落地时 0.14.0 **已打 tag**，则该页已是冻结的历史记录，改为**删除**这一节，并把移除记到 `.claude/notes/release.md` 的下一个版本标题下（`docs/` 没有下一版页面，且 `architecture.md:147` 禁止手写 changelog）。tag 状态按 CLAUDE.md § Release with molrs 的手工检查（`git tag` + 索引）判定，实现者在交付说明里写明走了哪条。

**托管块的回填要有人跑。** `architecture.md:24`/`:36`/`:60` 由 `/mol:map` 再生成，但再生成不是自动的（该块头 `:7` 仍写着 `Generated 2026-08-04 for release 0.12`）：本 spec 的最后一项任务是在落地后运行 `/mol:map`，或至少把这三行记进交付说明。

**行号是链前快照。** 03 与 04 只依赖 01、彼此无序，且都改 `CLAUDE.md` 与 `_replicas.py`；若 04 先落地，本文的行号可能漂移——一律按引文定位。

**不需要改的（已复核）。** `docs/index.md:166-168` 与 `docs/user-guide/09_packing.md:8` 说的是外部 molpack 包 / Packmol 程序本身，措辞已正确；`docs/api/pack.md` 整页是 molpack，无 mkdocstrings 指向 `molpy.pack`；`zensical.toml:82,152` 的两条导航指向这两个 molpack 页面，保留；`pyproject.toml` 用 `[tool.setuptools.packages.find]` 自动发现，删目录即可。

**Reuse decision.** 纯删除，不新增任何公开符号。逐项结论：`reuse molpy.core.region`（几何词汇的既有归属，本 spec 不动它，也不把 `constraint.py` 的任何形状搬过去）；`reuse molpack`（打包能力的生态归属，由文档指向，molpy 不产生 import 依赖）；无 `new` 符号——唯一新增文件是一个测试。

**形状合规。** 本 spec 不新增公开符号，CLAUDE.md § Shape check 无适用项；净效果是移除一个 all-in-one 外部管线 facade，方向与 § Forbid 一致。

## Files to create or modify

- `src/molpy/pack/__init__.py` (delete)
- `src/molpy/pack/constraint.py` (delete)
- `src/molpy/pack/target.py` (delete)
- `src/molpy/pack/packer/__init__.py` (delete)
- `src/molpy/pack/packer/base.py` (delete)
- `src/molpy/pack/packer/packmol.py` (delete)
- `tests/test_pack/test_constraint.py` (delete)
- `tests/test_pack/test_packmol.py` (delete)
- `tests/test_pack/test_target.py` (delete)
- `src/molpy/__init__.py` — 删除 `:29`、`:48`、`:317` 三处 `pack`
- `src/molpy/adapter/__init__.py` — 改写 docstring `:6-7`
- `src/molpy/builder/assembly/_replicas.py` — 改写 docstring `:3`
- `tests/test_init.py` (new) — `__getattr__` 契约的单元测试
- `.claude/notes/architecture.md` — 改 `:80`、`:83`、`:86`，删 `:84`（`:24`/`:36`/`:60` 属 `/mol:map` 托管块，不手改）；动词表「装填」行改为已成立、删除债务清单第 2 行
- `CLAUDE.md` — 删 `:199` 表行、改 `:224` 导入方向散文（"`io` may be imported by `builder`, `pack`, …"）、改 `:374` 测试树
- `README.md` — 删 `:64` 表行、改 `:28` 散文
- `docs/developer/architecture-overview.md` — 删 `:16` 表行、改 `:20`（去掉 packmol）
- `docs/getting-started/external-tools.md` — 改写 `:19`
- `docs/user-guide/15_mcp.md` — 改写 `:616-625`
- `docs/getting-started/migration-0-14.md` — 按 tag 状态改写或删除 `:102-108`
- `.claude/notes/release.md` — 仅在 v0.14.0 已打 tag 的分支下：在下一版本标题下记录移除

## Tasks

- [ ] Write a failing unit test in `tests/test_init.py` (new)：断言 `molpy` 对**不在 `_LAZY_SUBMODULES` 中的中性未知名**抛 `AttributeError`（`src/molpy/__init__.py:55-58` 的 `__getattr__` 契约）
- [ ] Delete `src/molpy/pack/`（6 个模块）与 `tests/test_pack/`（3 个测试文件）
- [ ] Remove the `pack` lazy entry from `src/molpy/__init__.py`（`:29` TYPE_CHECKING、`:48` `_LAZY_SUBMODULES`、`:317` `__all__`）
- [ ] Rewrite the stale in-source docstring references：`src/molpy/adapter/__init__.py:6-7` 改指 `molpy.wrapper` 范例，`src/molpy/builder/assembly/_replicas.py:3` 改指 molpack
- [ ] Drop `pack` from the layer-graph edges：`.claude/notes/architecture.md:80,83,86` 去掉 `pack`、`:84` 整行删除，`CLAUDE.md:224` 导入方向散文去掉 `pack`（不动 `architecture.md` 的 `/mol:map` 托管块）
- [ ] Drop the `pack` rows from the descriptive tables：`CLAUDE.md:199,374`、`README.md:28,64`、`docs/developer/architecture-overview.md:16` 及 `:20`（同时移除 `wrapper` 行里的 `packmol`）
- [ ] Rewrite the user-facing doc references to molpack：`docs/getting-started/external-tools.md:19`、`docs/user-guide/15_mcp.md:616-625`、`docs/getting-started/migration-0-14.md:102-108`（按 v0.14.0 tag 状态就地改写或删节）
- [ ] Flip the 装填 row of the verb table in `.claude/notes/architecture.md` § 设计铁律 4 to 已成立 and delete debt row 2 (`Packer.__call__` + `molpy.pack`), as the table's retirement clause requires
- [ ] Run `/mol:map` after landing (or record the three stale managed rows `architecture.md:24/:36/:60` in the delivery summary) so the blueprint stops describing a deleted package
- [ ] Run full check + test suite

## Testing strategy

依 `.claude/notes/testing.md`：单元测试、一测一行为、路径镜像源码。

- **删除面。** `tests/test_pack/` 三个文件随实现一起删除（被测对象不存在，测试即无主）。不保留任何改写后的替身测试。
- **唯一新增测试，且不以 `pack` 命名。** 新建 `tests/test_init.py`（镜像 `src/molpy/__init__.py`），**只有一个**测试，断言 `__getattr__`（`src/molpy/__init__.py:55-58`）对一个中性的、不在 `_LAZY_SUBMODULES` 里的名字抛 `AttributeError`：

  ```python
  def test_unknown_attribute_raises():
      with pytest.raises(AttributeError):
          molpy.not_a_submodule
  ```

  这测的是 facade 的**契约**（未注册的名字一律 `AttributeError`），对任何被删名字都成立，删除 `pack` 只是它的一个实例。
- **明确不写的测试。** 不写 `"pack" not in mp.__all__`、不写 `mp.pack` 专门抛错：那是对"某个名字不存在"的导出清单闸（`.claude/notes/testing.md:12-15` 把 import-all / source-text 闸排除在套件外），也等于为一个已删名字长期背一条向后兼容守卫（CLAUDE.md § Iron law：不要把已知 rot 以守卫形式留下）。**移除这件事由 `migration-0-14.md` 与版本 tag 记录，不由测试记录**；静态面由 ac-002 把关。
- **既有模式。** 负向导出断言在本仓虽有先例（`tests/test_engine/test_lammps_relax.py:56`），本 spec 刻意不沿用，理由同上。
- **单测绿灯口径。** `uv run --extra dev python -m pytest tests/test_init.py`。
- **无域校验。** 本 spec 不引入物理量、方程或参考值。
- **无 regression 示例。** `.claude/notes/testing.md:13` 明确把 `regressions/` 目录列为"不在本套件内"，仓库当前确无该目录；本 spec 不创建，交付验证由单测 + 全量 `check` / `pytest` 承担（ac-003 / ac-010）。

## Out of scope

- molpack 内部的一切（API、命名、restraint 形状）——本 spec 只把文档指向它，不产生 import 依赖，也不把 `molcrafts-molpack` 加进任何 extras。
- `compute/`（sub-spec 02）与 `builder/`（sub-spec 04）的动词整改。
- `api-verb-unification` 动词表本身（sub-spec 01 `api-verb-unification-01-verbtable`）；本 spec 只保证 `pack` 一格落到 `molpack.pack` 时表内引用不再指向 molpy 符号。
- `src/molpy/core/region.py`：几何词汇的既有归属，不动、不吸收 `constraint.py` 的任何形状。
- `.claude/notes/architecture.md` 的 `/mol:map` 托管块（`:1-68`，含 `:24`、`:36`、`:60` 三处 `pack`）：由 `/mol:map` 重新生成，本 spec 不手改。
- `docs/index.md:166-168`、`docs/user-guide/09_packing.md`、`docs/api/pack.md`、`README.md:159`：已复核，指的是外部 molpack / Packmol 程序，正确，不改。
- 把 `Packmol` 搬进 `wrapper/` 以修正 CLAUDE.md:201 的分类矛盾：删除使该矛盾消失，不需要搬迁。
- 任何 deprecated shim、别名或 `DeprecationWarning` 过渡期：stage `experimental` + 零消费者，直接删。
- 新建 `regressions/` 目录（被项目测试法明文排除）。
