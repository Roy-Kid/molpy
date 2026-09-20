---
slug: api-verb-unification-01-verbtable
title: API 动词统一 (1/4)：把「变换族 → 动词」表写入设计铁律 4
status: approved
created: 2026-09-20
chain: api-verb-unification
---

# API 动词统一 (1/4)：把「变换族 → 动词」表写入设计铁律 4

## Summary

`.claude/notes/architecture.md` § 设计铁律 › 4 已经写着「同族变换用**同一个动词**」，但只给了两个例子、没有表，于是这条铁律无法被 review 或 agent 判定：谁和谁算一族、一族该用哪个动词，全靠读者猜。更糟的是那两个例子本身就是反例——`GraphAssembler.assemble`（`builder/assembly/_assembler.py:116`）与 `VirtualSiteBuilder.apply`（`builder/virtualsite.py:69`）同属「图 → 图」一族却用了两个动词，规则的示例在打规则自己的脸。本 sub-spec 只改知识文档：把链上已裁定的十一族写成完整的表，**逐行标注它现在是否成立**——分析族在 sub-spec 02 落地前是已声明的债，装填族在 sub-spec 03 落地前是已声明的债，图变换族的 `GraphAssembler.apply` 由 sub-spec 04 兑现；补六条脚注、`Compute` 归属裁定、`__call__` 政策与一份「已声明的债」清单；顺手把 § 设计铁律 标题的「四条」改成「六条」，并把 CLAUDE.md 里两处复述动词的地方改成指针。落地后每一次重命名都能指着表里的一行说话，`/mol:review --axis=arch` 也有了可判定的依据。本 spec 不碰 `src/`、不碰 `tests/`。

## Domain basis

不适用。本 spec 不声明任何物理：它写的是仓内命名宪法，无方程、无单位、无文献。

## Design

**宪法定位。** 仓内没有 `.claude/notes/law.md`；实际的宪法是 CLAUDE.md § Design preferences 加 `.claude/notes/architecture.md` § 设计铁律 1–6。本 spec 不新增任何符号、不写一行 Python，铁律 1（不硬编码字段名）、2（体系判断留 molpy）、5（不 fallback）在此不被触发；它要做的正是把**铁律 4** 从一句口号补成一张可判定的表。铁律 3（experimental 阶段不背向后兼容）是后续 sub-spec 能直接改名删包、不留 shim 的授权来源，表里因此不含任何 deprecated 别名列。铁律 6（molrs 是实现细节）与表的「依据」列同向，且给出本表的**单向性**：凡 molrs 已经定下的动词，molpy 一律沿用；本表对 molrs 只作**描述**，不构成约束——sink direction 是单向的。

**诚实性原则（本次修订的核心）。** 宪法不得写一句当下为假的话。表里有三族在写入时**尚未成立**：分析族（全部子类仍是 `__call__`，`compute/base.py:51` 的抽象方法就是 `__call__`，全树零个 `def compute`）、装填族（`molpy.pack` 整包仍在，`Packer.__call__` 在 `pack/packer/base.py:69` 与 `Packer.pack`（`:50`）并存）、图变换族（`assemble` 尚未改名）。处理方式不是回避，而是**逐行标注承载者**：每一行写明「已成立」或「落地前为已声明的债（sub-spec NN）」，并在表后给一份集中的债务清单。这样表从第一天起就是真的——它描述的是「裁定」，并明确区分哪些裁定已兑现。

**族的边界定义（写在表的正上方）。** 一个方法属于哪一族，**只看它吃什么、吐什么**，不看它所在类的名字或后缀。本表**不改任何类名**，类名不在本链的范围内：

- `VirtualSiteBuilder`（`builder/virtualsite.py:61`）、`DrudeBuilder`（`:92`）、`Tip4pBuilder`（`:173`）顶着 `Builder` 后缀，但它们吃一个已有的 `Atomistic`、吐一个被改写的 `Atomistic`——**图变换族**，动词 `apply`。
- `GrapheneBuilder`（`nanostructure/graphene.py:15`）、`CarbonTubeBuilder`（`carbon_tube.py:15`）吃参数、吐新结构——**构造族**，动词 `build`。

同一个后缀落在两族里，正说明后缀不是判据。下一次改名照「输入 → 输出」判，不照名字判。

**写入位置。** 表进 `.claude/notes/architecture.md` § 设计铁律 › 4「OOP,且是真 OOP」（第 149 行起），替换现有的这一条：

```
- 同族变换用**同一个动词**:`GraphAssembler.assemble`、`VirtualSiteBuilder.apply`。
  不要一个类叫 `build`、一个叫 `apply`、一个叫 `run`。
```

保留该条的首句（规则本身），把两个例子换成边界定义 + 表 + 脚注 + 裁定 + 债务清单。§ 4 其余部分（假 OOP 的禁令、门面 vs 真类的对照表、边界切错的信号）逐字不动；§ 1/2/3/5/6 与 managed 块（`<!-- mol:map:managed begin/end -->`，第 3–68 行）一个字节都不动。唯一的例外是第 113 行的标题行（见下）。

**要写入的表（逐字；`molrs:` 前缀 = 对 molrs 现状的描述，不约束 molrs）：**

| 变换族 | 输入 → 输出 | 动词 | 成员 | 状态 | 依据 |
|---|---|---|---|---|---|
| 构造 | 配方/IR/参数 → 新结构 | `build` | `PolymerBuilder.build(topology)`(`builder/assembly/_polymer.py:87`)、`Lattice.build(region)`(`builder/crystal.py:210`)、`GrapheneBuilder.build`(`nanostructure/graphene.py:56`)、`CarbonTubeBuilder.build`(`nanostructure/carbon_tube.py:65`)、`AmberPolymerBuilder.build`(`polymer/ambertools/amber_builder.py:166`) | 已成立 | `molrs:` 侧已是 `build`(`molrs.builder.*`、`NeighborList.build`);molpy 永不改写 molrs 的动词(sink direction) |
| 图变换 | 已有图 → 被改写的图 | `apply` | `StructureFinalizer.apply`(`builder/_finalize.py:43`)、`VirtualSiteBuilder.apply`(`builder/virtualsite.py:69`)、`GraphAssembler.apply`、`molrs:Reaction.apply` | **部分**:`GraphAssembler` 现名 `assemble`(`builder/assembly/_assembler.py:116`),改名由 **sub-spec 04**(`api-verb-unification-04-assembler`)兑现;**04 落地时必须把本格改为「已成立」并删除债务清单第 3 行** | 签名形状相同;`molrs:Reaction.apply` 是既有成员(描述性) |
| 分析 | frames/arrays → Result | `compute` | `molpy.compute` 下每一个实现 `compute()` 的分析类(02 落地前它们是 `Compute` ABC 的子类,落地后不继承任何基类;计数见脚注 5,勿写死)+ `molrs:` 每个 kernel | **落地前为已声明的债(sub-spec 02**,`api-verb-unification-02-compute`**)**:当前全部子类实现 `__call__`(抽象方法在 `compute/base.py:51`),全树零个 `def compute`;**02 落地时必须把本格改为「已成立」并删除债务清单第 1 行** | `molrs:molrs.compute.protocol.Compute` 只认 `compute`;`protocol.py` 模块 docstring 明说 `__call__` / `dump()` 不在契约内 |
| 装填 | targets → Frame | `pack` | 外部 `molpack`(`docs/api/pack.md`) | **落地前为已声明的债(sub-spec 03**,`api-verb-unification-03-pack`**)**:`molpy.pack` 整包仍在,删除集见脚注 4;**03 落地时必须把本格改为「已成立」并删除债务清单第 2 行** | molpy 及全部兄弟仓零消费者,文档已指向 molpack |
| 分型 | 图 → 带类型的图 | `typify` | `molrs:Typifier`;`typifier/base.py:99` 是参考范式(继承 molrs 基类、保留动词、加一个 hook,`typify` 在 `:106`) | 已成立 | molrs 所有(描述性) |
| 3D 生成 | 图 → 坐标 | `generate` | `Conformer.generate`(`conformer/__init__.py:42`)——**仅此一个**,其余同词干的名字见脚注 1 | 已成立 | `molrs:molrs.conformer.Conformer` 所有(描述性) |
| 外部进程 | 文件/体系 → 文件/轨迹 | `run` | `Wrapper.run`(`wrapper/base.py:74`)、`Engine.run`(`engine/base.py:212`)、`molrs:LBFGS.run` | 已成立 | 「跑一个东西」,不是数据变换;不动 |
| 发射 | Frame → 引擎输入 | `emit` | `Emitter.emit`(`io/emit/__init__.py:32`)及各格式实现(`io/emit/{lammps,gromacs,openmm,xml}.py`) | 已成立(自由函数形态的债见债务清单) | 不动 |
| 读写 | 路径 ↔ 对象 | `read` / `write` | `io` | 已成立 | 不动 |
| 选择 | context/struct → 子集 | `select` | 入口动词 `select`:`Selector.select`(`builder/assembly/_selector.py:35`)、`VirtualSiteBuilder.select`(`builder/virtualsite.py:78`)、`Atomistic.select`(`core/atomistic.py:270`) / `CoarseGrain.select`(`core/cg.py:189`)(核心数据模型上的两处只作描述,其命名归 § Graph sink decisions,本表不授权改名);`MaskPredicate.mask(block) -> ndarray`(`core/selector.py:24`)是**生产者 hook**,不是第二个入口动词 | 已成立(名字撞车的债见债务清单) | 不动 |
| 放置 | struct → 就地坐标 | `place` | `Placer.place`(`builder/assembly/_placer.py:46`)、`ResiduePlacer.place`(`:75`);经 `builder/assembly/__init__.py:16` 与 `builder/__init__.py:42,44` 导出 | 已成立 | 领域动词,不在两族之内,见脚注 6 |

**`Compute` 归属裁定（逐字，紧随表后）：** 仓里同时活着两个 `Compute` 概念——`molpy.compute.base.Compute`（`compute/base.py:18`，一个以 `__call__` 为抽象方法的 ABC）与 `molrs.compute.protocol.Compute`（只认 `compute`）。**唯一的 `Compute` 是 molrs 那个**，sub-spec 02 落地后用户经 identity re-export 以 `molpy.compute.Compute` 取得（铁律 6 的第一档手段；今天 `compute/__init__.py:19` 还是 `from .base import Compute`）。molpy 的 ABC 并非空壳——它带 `__init__(**config)`（`compute/base.py:41-48`）与 `dump()`（`:65-71`），按铁律 6 的三档表那本该是「继承」档；本链**裁定不继承而是删除**：`dump()` 在 `src/`、`tests/` 中零调用（唯一引用是 `docs/developer/extending-compute.md:69`），`**config` 只喂 `dump()`，二者都是没有消费者的公开面。这是 **sub-spec 02 名下的一次 breaking change**（`dump()` 与 `**config` 从公开面消失，`extending-compute.md` 随之改写），不是无损再导出；02 的 ac-003 记录它。不留 `__call__` 糖、不留转发壳。`.claude/specs/INDEX.md:12-14` 把这件事记为 release-0-14 遗留的 open decision（"the callable `compute.base.Compute` shells versus the molrs `Compute` Protocol (verb unification)"）——**本链就是它的关闭动作**，sub-spec 02 落地后该条从 open 转 closed。

**一族一动词，不是一类一动词（新增一句）：** 一个类可以在它参与的**每一个**族里各有一个动词——`PolymerBuilder` 自己有 `build`（构造族，`_polymer.py:87`），并从 `GraphAssembler` 继承 `apply`（图变换族）。被禁止的是**同一族里两个动词**，不是一个类上有两个动词。

**`assemble → apply` 裁定（逐字）：** `GraphAssembler.assemble` → `apply` 而非 `build`：`PolymerBuilder(GraphAssembler)`（`_polymer.py:38`）已有 `build(topology)`（`:87`，在 `:94` 调 `self.assemble`），父类再挂 `build(world, selector)` 会撞签名（LSP，`ty` 报错）；而 `assemble` 吃一个已有的 world、吐一个被改写的 world——图 → 图，本就属 `apply` 族。

**脚注（逐字）：**

1. 本表约束的是变换的**入口动词**;抽象模板 hook 与描述性方法名不在族内,保留原名:`VirtualSiteBuilder.build_sites`(`builder/virtualsite.py:82`/`:126`/`:205`)、`engine/openmm.py:303 generate_inputs`、`builder/polymer/sequences.py:34`/`:74`/`:138`/`:182 generate_sequence`、`adapter/rdkit.py:105 generate_3d`。(`pack` 侧的 `generate_input_*` 随 sub-spec 03 整包消失,不必单独裁决。)
2. **`build_*` 全树普查。** 同词干的一行捷径,保留:`PolymerBuilder.build_sequence`/`build_linear`/`build_ring`/`build_star`(`_polymer.py:96`/`:100`/`:110`/`:114`)、`builder/ambertools.py:241 build_polymer`、`builder/polymer/system.py:146 build_chain`。模板 hook,保留:`build_sites`(见脚注 1)。**债**:`parser/moltemplate/builder.py:551 build_forcefield` 与 `:967 build_system` 是自由函数形态的工厂,正撞 CLAUDE.md § Forbid「Factory functions as the primary constructor story」与铁律 4「自由函数不是公开面」——见债务清单,单独 `/mol:refactor`,不属本链。
3. `DistributionIR.build()`(`builder/polymer/distributions.py:28`)不吃数据——它是 factory 式构造器,不是变换族成员;改名或删除是另一次 `/mol:refactor`(维护者裁定 OQ-3 记为「只留脚注」)。
4. **`__call__` 政策。** 一族一个动词;`__call__` 只出现在「这个对象**就是**一个函数」的地方。本链全部落地后,唯一豁免是 `core/selector.py` 的 `MaskPredicate.__call__(block) -> Block`(`:26`)与 `mask(block) -> ndarray`(`:24`)并存——返回类型不同,且谓词要用 `& | ~` 组合(`:30`–`:40`)。其余把 `__call__` 当作同一操作第二个名字的,一律删:(a) 全部 compute 壳,由 sub-spec 02 改名为 `compute`;(b) `Packer.__call__`(`pack/packer/base.py:69`)与 `Packmol.__call__`(`pack/packer/packmol.py:56`),随 sub-spec 03 **整包删除**——删除集为 `src/molpy/pack/` 全包(`Packmol`、`Packer`、`Target`、全部 `*Constraint`)、`tests/test_pack/`、以及 `src/molpy/__init__.py:29`/`:48`/`:317` 三处 facade 登记,并**改写**(不是删除)两处会悬空的 docstring 交叉引用 `src/molpy/builder/assembly/_replicas.py:3` 与 `src/molpy/adapter/__init__.py:6-7`。**在 03 落地之前,「唯一豁免」这句话尚不为真**——那之前 `Packer.__call__` 是 sub-spec 03 名下的已声明的债,见债务清单。
5. **不写死计数。** 本表不记「N 个 compute 壳」这类会腐烂的数字;需要时现场数,按动词而不按基类:`rg -n 'def compute\(' src/molpy/compute`(02 落地前该命令给 0——那时壳还叫 `__call__`,用 `rg -n 'class \w+\(Compute\)' src/molpy/compute` 数基类形式,2026-09-20 给 34 行,其中 33 行是真实子类,另一行是 `compute/base.py:27` docstring 里的示例——**读 grep 输出要减掉 doctest**)。
6. **族的边界只由「输入 → 输出」决定,与类名/后缀无关,本表也不改任何类名;核心数据模型 API 整体在两族之外。** 本表只裁方法名,不裁类名:`VirtualSiteBuilder`/`DrudeBuilder`/`Tip4pBuilder` 带 `Builder` 后缀却是图变换族成员(动词 `apply`),`CarbonTubeBuilder`/`GrapheneBuilder` 是构造族;类名改不改是另一次决定。`MonomerLibrary.expand(topology)`(`builder/assembly/_library.py:57`)、`Replicas.grid`/`.times`(`builder/assembly/_replicas.py:40`/`:89`)与 `Placer.place`/`ResiduePlacer.place` 是领域动词,不在构造/图变换两族之内。**判据**(不是例子清单):`Atomistic`/`CoarseGrain`/`Frame` 的核心数据模型 API——`def_*`/`del_*`/`copy`/`merge`(`core/atomistic.py:462`)/`move`/`rotate`/`scale`/`align`(`:534`)/`replicate`(`:559`)/`extract_subgraph`(`:377`)/`to_frame`——虽然也是「已有图 → 改写后的图」,但它们是数据模型自身的动作,由 § Graph sink decisions(锁定)与 CLAUDE.md § What must never change casually 管辖,**不属于图变换族**;图变换族只收「以整张图为输入、产出新图的领域变换入口」(finalize、virtual sites、assembly、reaction)。

**已声明的债（逐字，表后独立小节「本表尚未兑现的部分」，9 行）：**

| 债 | 位置 | 归属 |
|---|---|---|
| 分析族全部实现 `__call__`,无 `def compute`;molpy 自带 `Compute` ABC 与 molrs Protocol 并存 | `compute/base.py:18`/`:51` + 全部子类 | **sub-spec 02** `api-verb-unification-02-compute` |
| `Packer.__call__` 与 `Packer.pack` 并存;`molpy.pack` 整包仍在 | `pack/packer/base.py:50`/`:69`、`pack/packer/packmol.py:56` | **sub-spec 03** `api-verb-unification-03-pack` |
| 图变换族入口仍叫 `assemble` | `builder/assembly/_assembler.py:116` | **sub-spec 04** `api-verb-unification-04-assembler` |
| `Selector` 一名两义:`core/selector.py:43 Selector = MaskPredicate`(掩码谓词)与 `builder/assembly/_selector.py:26 class Selector(ABC)`(装配选择器)是两个无关的类型 | 同左 | 单独 `/mol:refactor`,**不在本链** |
| 模块级自由函数 `emit(name, …)` 注册表 dispatcher 与 `emit_python`,违铁律 4「自由函数不是公开面」 | `io/emit/__init__.py:51`、`parser/moltemplate/py_emitter.py:42` | 单独 `/mol:refactor`,**不在本链** |
| 自由 `build_*` 工厂(见脚注 2) | `parser/moltemplate/builder.py:551`/`:967` | 单独 `/mol:refactor`,**不在本链** |
| § 4「边界切错的信号」(`architecture.md:173-174`)在 `_polymer.py:93-94` 上仍在鸣响:`expand(topology)` 之后又 `TopologySelector(topology)`,同一份数据穿两个对象 | `builder/assembly/_polymer.py:93-94` | 单独裁决(让 `MonomerLibrary.expand` 连同配对规则一起产出,或 `TopologySelector` 从展开后的 world 推导),**不在本链**;sub-spec 04 已记为 found-not-fixed |
| CLAUDE.md managed 块 `:71` 复述动词且两例皆假(`NeighborList.build` 实为 `__call__`→02 后为 `compute`;`ForceField.energy` 不存在) | `CLAUDE.md:71`(managed) | 链后由操作者 `/mol:bootstrap` 重生成 |
| managed 块 § Style summary(`architecture.md:43-49`)复述命名/构造约定,与本表两个副本 | `architecture.md:43-49`(managed) | 链后 `/mol:map`:Style summary 不再复述动词裁定,动词一栏只指向 § 设计铁律 4(其余 naming/construction/errors 栏照常) |

本清单是**记账与路由**,不是放行:每一行要么有承载 sub-spec、要么有明确的下一次 `/mol:refactor` / `/mol:bootstrap` / `/mol:map`;不得把任何一行读作「本表允许的例外」。

**给下一次 `/mol:map` 的一句话（写在 § 4 末尾，managed 块之外）：** 本文件 managed 块里的 § Style summary（`architecture.md:43-49`）仍在复述命名/构造约定；下一次 `/mol:map` 生成的 Style summary 不得复述动词裁定，动词一栏只指向 § 设计铁律 4（librarian 需要的 naming/construction/errors 其余栏照常），避免同一套约定有两个副本（债务清单第 9 行）。

**CLAUDE.md 的两处改动。** 二者是同一件事——让 CLAUDE.md 停止复述动词、改为指向本表：

1. `## Architecture Overview` 段内、import-direction 段落之后加一行纯指针（**不列任何动词名**，否则它本身就会随 02/03/04 腐烂）：
   `Transformation verbs: one verb per transformation family; the binding table (families, members, declared debt, owning sub-spec) is `.claude/notes/architecture.md` § 设计铁律 4. Do not restate verbs here.`
2. `### builder module` 里 `:473-474` 这一条（两行）复述了动词、且点名了即将改名的 `assemble`：
   `- Construction is a method on the owning type (`Lattice.build`, `PolymerBuilder.build`, `GraphAssembler.assemble`); there are no free `build_*` / `create_*` factories`
   改写为指针：
   `- Construction and transformation verbs: look them up in the family→verb table (`.claude/notes/architecture.md` § 设计铁律 4), which also records the declared debt and its owning sub-spec.`
   **顺带修掉一句当下为假的话**：原句断言「there are no free `build_*` / `create_*` factories」，而 `parser/moltemplate/builder.py:551`/`:967` 正是两个自由 `build_*` 工厂（脚注 2 与债务清单已记）。新句不再作这个断言，改写后在 sub-spec 04 前后都成立。

两处均在 `<!-- mol:bootstrap:managed end -->`（CLAUDE.md:133）**之后**，属自由区，`/mol:bootstrap` 重跑会原样保留。CLAUDE.md 是操作者拥有的配置：与 sub-spec 04 同一口径，`/mol:impl` 落地时把这两处 diff 呈交操作者确认后再落。

**跨仓指针。** `.claude/notes/cross-repo-spec-map.md` 加一行，把表的单向性写在跨仓文档里（molrs 侧的人会读这个文件）：表中标 `molrs:` 的行是对 molrs 现状的**描述**，不约束 molrs；molpy 沿用 molrs 的动词，反向不成立（sink direction）。

**Reuse decision.** 本 sub-spec 不设计任何新符号——它不写 `src/`，「复用 vs 新建」在符号层面无从谈起。等价判断落在文档结构上：**复用**既有的 `.claude/notes/architecture.md` § 设计铁律 › 4 小节（把表补进这条已经存在的规则里），**不**新开 `.claude/notes/api-verbs.md`；CLAUDE.md 与 cross-repo-spec-map.md 各只留一行指针，**不**复制表——单一事实源在 notes。这与 CLAUDE.md 既有的「Import direction (full table in `.claude/notes/architecture.md`)」同一范式，新指针紧随其后，读起来与现有那行同构。librarian 的放置结论（表落 architecture.md:155；CLAUDE.md managed 块不得承载；notes.md 只放窄决策）按此采纳。

## Files to create or modify

- `.claude/notes/architecture.md` — (a) § 设计铁律 › 4 内第 155–156 行的两例改写为表 + `Compute` 归属裁定 + 一族一动词条款 + `assemble → apply` 裁定 + 六条脚注 + 债务清单 + 给 `/mol:map` 的一句话；(b) 第 113 行标题行「四条硬约束」→「六条硬约束」。managed 块（3–68）与 § 1/2/3/5/6 不动
- `CLAUDE.md` — `## Architecture Overview` 段（第 ~228 行后）加一行指针；`### builder module` 第 473–474 行整条改写为指针（不再出现 `Lattice.build` / `PolymerBuilder.build` / `GraphAssembler.assemble` 任何一个）并删去当下为假的「no free `build_*`」断言。managed 块（25–133）不动
- `.claude/notes/cross-repo-spec-map.md` — 加一行指向动词表、并声明 `molrs:` 行为描述性、不约束 molrs

## Tasks

- [ ] Write the family-boundary clause (family decided by input→output, class names out of scope and unchanged) plus the 变换族 → 动词 table (11 rows, 6 columns incl. 状态) into `.claude/notes/architecture.md` § 设计铁律 › 4, replacing the two-example line 155–156, marking each row 已成立 or 已声明的债 with its owning sub-spec (02 / 03 / 04), each debt cell carrying the retirement clause ("NN 落地时改为已成立并删除债务行 k"), and prefixing molrs-descriptive members with `molrs:`
- [ ] Write the `Compute` ownership ruling below the table (one `Compute` = `molrs.compute.protocol.Compute` via identity re-export; `compute/base.py` deleted by sub-spec 02; cite `.claude/specs/INDEX.md:12-14` as the open decision this chain closes)
- [ ] Write footnotes 1–6 below the ruling (template hooks + `generate_*` exclusions; the full `build_*` census; `DistributionIR.build`; the `__call__` policy with sub-spec 03's exact deletion set and the carve-out's truth condition; the no-frozen-count grep with its doctest caveat; class names and the domain verbs `expand` / `grid` / `times` / `place` outside the two families)
- [ ] Write the 「本表尚未兑现的部分」 debt table (9 rows: compute shells → 02, `Packer.__call__` + `molpy.pack` → 03, `assemble` → 04, `Selector` name collision, free `emit` dispatcher + `emit_python`, free `build_*` factories → separate `/mol:refactor`, the `_polymer.py:93-94` boundary signal → separate decision, `CLAUDE.md:71-72` → `/mol:bootstrap`, the managed Style summary → `/mol:map`) with the retirement clause on rows 1–3 and the 「记账与路由,不是放行」 sentence
- [ ] Write the `assemble → apply` ruling, the one-verb-per-family-not-per-class clause, and the sentence asking the next `/mol:map` to reduce the managed Style summary (`architecture.md:43-49`) to a pointer at the table
- [ ] Fix the section header at `.claude/notes/architecture.md:113` from 「四条硬约束」 to 「六条硬约束」
- [ ] Update `CLAUDE.md` to stop restating verbs: add the pointer line in § Architecture Overview after the import-direction paragraph, and rewrite the whole bullet at lines 473–474 to point at the table, naming no verb (dropping the false "no free `build_*` / `create_*` factories" claim)
- [ ] Add the one-line pointer to `.claude/notes/cross-repo-spec-map.md` marking `molrs:` rows descriptive and non-binding on molrs
- [ ] Verify the diff touches only these three files, leaves both managed blocks byte-identical, confines `architecture.md` hunks to line 113 and § 设计铁律 › 4, and contains no file under `src/` or `tests/`
- [ ] Run full check + test suite

## Testing strategy

本 spec 只改三份知识文档，不新增、不修改任何 `src/molpy/` 符号，因此**不写单测**——按 `.claude/notes/testing.md` § Test shape，`tests/` 里只放针对单个函数/方法的单测，且明令排除 source-text / docs-block gates；给一段 markdown 写断言正是被排除的那一类。同理**不建 `regressions/` 目录**：该 note 把 `regressions/` 直接列在「不在 suite 里」的清单上（`testing.md:16`），项目规则压过 spec 模板的默认要求，本 spec 因此没有 regression-example 任务、也没有对应的 regression 判据（ac-016 是门禁判据，不是 regression 判据）。

验收改由 `type: docs` 判据承担，全部可由第三方对着三个文件逐条肉眼判定：表的十一行是否齐、六列是否齐、三族的债务标注是否点名了正确的 sub-spec slug、`Compute` 归属裁定是否引了 INDEX.md、计数是否已从宪法里拿掉、脚注与债务清单是否齐、CLAUDE.md 两处是否都不再复述动词、跨仓指针是否声明了单向性。

回归保护由 diff 边界提供，不由测试提供：

- 边界检查 — `architecture.md` 的 `<!-- mol:map:managed begin/end -->` 块与 CLAUDE.md 的 `<!-- mol:bootstrap:managed begin/end -->` 块前后字节一致；§ 设计铁律 1/2/3/5/6、§ Import-direction rules、§ Graph sink decisions 不变；`architecture.md` 的改动只落在第 113 行与 § 4 内。
- 溢出检查 — diff 里不出现任何 `src/` 或 `tests/` 下的文件（本链的代码改动全在 02/03/04）。
- 门禁不回退 — `uv run --no-project --with 'tox>=4.23' --with ruff==0.16.1 --with ty==0.0.65 tox -e lint` 与 `uv run --extra dev python -m pytest tests/ -n auto` 仍全绿（与改动前同一结果；markdown 改动本不入 ruff/ty 视野，绿灯即证明没有代码泄漏进本次改动）。

## Out of scope

- **不改任何类名。** `VirtualSiteBuilder` / `DrudeBuilder` / `Tip4pBuilder` 的 `Builder` 后缀与它们所属的图变换族不一致，这是本表明确接受的现状——族由「输入 → 输出」判，后缀不是判据；类名重命名若真有必要，是另一次 `/mol:refactor`。
- **不改领域动词。** `MonomerLibrary.expand`、`Replicas.grid` / `.times`、`Placer.place` 已裁定在两族之外，保留原名，不在本链也不在后续 sub-spec 的改名范围内。
- **不动 `src/` 与 `tests/`。** 实际改名与删除由 sub-spec 02（`api-verb-unification-02-compute`）、03（`api-verb-unification-03-pack`）、04（`api-verb-unification-04-assembler`）执行；本 spec 只立表、只记债。
- **不修债务清单里的任何一笔债。** `Selector` 一名两义（`core/selector.py:43` vs `builder/assembly/_selector.py:26`）、自由 `emit` dispatcher（`io/emit/__init__.py:51`）与 `emit_python`（`py_emitter.py:42`）、自由 `build_forcefield`/`build_system`（`parser/moltemplate/builder.py:551`/`:967`）——本 spec **记录并路由**到各自的 `/mol:refactor`，不在本链修；这是「不背默认忽略」的记账动作，不是放行。
- **不刷新 architecture.md 的 managed 块。** 它已陈旧（仍列已删除的 `Workflow`、缺 `md/` 与 `integrations/`，§ Style summary 43–49 与新表重复），由本链全部落地后的 `/mol:map` 统一重生成；本 spec 只在 § 4 里留下对它的期待，一个字节都不进那个块。
- **不改 CLAUDE.md managed 块内任何内容**（第 25–133 行）；`## Core Packages` 表里的 `pack` 行由 sub-spec 03 连同整包删除时处理。
- **不改 `docs/`。** `docs/api/pack.md` 已指向 molpack；其文档清理属 sub-spec 03。
- **不建 `regressions/`、不加 docs 门禁测试、不写单测**（理由见 Testing strategy）。
