---
slug: api-verb-unification-02-compute
title: API 动词统一 (2/4)：compute 家族动词统一到 compute()
status: approved
created: 2026-09-20
chain: api-verb-unification
depends_on:
  - api-verb-unification-01-verbtable
---

# API 动词统一 (2/4)：compute 家族动词统一到 compute()

## Summary

`molpy.compute` 目前有两套调用习惯：molpy 自己拥有的 shell 把数据输入挂在 `__call__`（`RDF(...)(frames, neighbors)`），而同一个命名空间里 re-export 的 molrs kernel 用 `.compute(...)`——`compute/__init__.py:14` 的模块 docstring 自己就在教 `EinsteinConductivity().compute(M, …)`。本规格统一 **molpy 自己拥有的那一族**：17 个模块里的 31 个 shell 一律改名为 `def compute(...)`，`voronoi.py` 的两个 shell 因为包装的 molrs 类根本不用 `compute` 这个动词而直接退化为 identity re-export；molpy 自己声明的 `compute/base.py`（`Compute(ABC)` + `_config` + `dump()`）整个删除，`Compute` 改为 molrs 那个 `@runtime_checkable` Protocol 的身份 re-export，使 `molpy.compute.Compute is molrs.compute.Compute`。这**不**意味着 `molpy.compute` 从此只有一个动词：molrs 拥有的符号保留 molrs 的动词（`Onsager.correlation`、`Persist.pair_survival_tcf`、`LinearFit.fit`、`CumulativeTrapezoid`、`DebyeFit`，以及本规格新暴露的 `RadicalVoronoi.build` / `VoronoiIntegration.integrate`），它们是 sub-spec 01 动词表里 `molrs:` 的描述性行，不受本规格约束。数值不动：shell 转发的 kernel、参数与顺序都不变，既有单元测试的输入和期望值原样保留，只换动词。

## Domain basis

N/A — 本规格不引入、不修改任何物理或数值：shell 以相同参数转发同一 kernel，`voronoi.py` 转 re-export 后调用的是原来就在跑的那个 molrs 对象。数值保真由「既有测试的输入与期望值逐字保留」承载（ac-006）。

## Design

**契约的唯一所有者是 molrs。** `molrs.compute.protocol.Compute` 是 `@runtime_checkable` 的 `Protocol`，只有一个方法 `compute(self, *args, **kwargs)`，`protocol.py` 的模块 docstring 已明说 `__call__` 和 `dump()` 不属于该契约。molpy 不再并行声明一份：`src/molpy/compute/base.py` 删除，`src/molpy/compute/__init__.py:19` 的 `from .base import Compute` 换成 `from molrs.compute import Compute`，`"Compute"` 留在 `__all__`（架构铁律 6：用户需要的每个 molrs 符号都必须能从 `molpy` 拿到；`compute/__init__.py:26-29,72-104` 已是这个形状）。

**shell 的形状。** 参照 `src/molpy/typifier/base.py:99-124` —— 继承 molrs 的基类、沿用 molrs 的动词（`typify`）、只加一个 hook，不在 molpy 复述契约。这里唯一的差别是 molrs 侧是 `Protocol` 而不是 ABC，所以 shell 不写继承：31 个 shell 从 `class RDF(Compute):` 变成 `class RDF:`，靠**结构一致**（定义了 `compute`）满足 `isinstance(shell, Compute)`。显式继承 Protocol 会把它退化成 ABC，等于又在 molpy 侧复述一遍契约。

**每个 shell 的改动是机械的三步**：`def __call__` → `def compute`（签名、参数名、返回标注一字不改）；删掉 `super().__init__(**config)`；删掉 `from .base import Compute` 与基类声明。构造参数照旧存成实例属性；`_config` 与 `dump()` 随 `base.py` 消失——全仓 `dump()` 的唯一调用点是 `docs/developer/extending-compute.md:69`，`src/` 与 `tests/` 无调用者。这是本规格名下的一次 **breaking change**（`dump()` 与 `**config` 从公开面消失，sub-spec 01 的 `Compute` 归属裁定已记录），不是无损再导出。shell 包装的 kernel 本来就暴露 `.compute(...)`（`src/molpy/compute/msd.py:33`），多数改名是纯直通。

**`**config_kwargs` 一并删除。** `src/molpy/compute/dielectric.py:130`、`:462` 收 `**config_kwargs`，而它唯一的去处就是 `super().__init__`（`:140`、`:471`）。`super().__init__` 一走，这个参数就变成**静默水槽**——用户拼错关键字不再报错（铁律 5：契约被违反 ⟹ 立刻 raise，永远不要悄悄换一条路走）。全仓 grep 只有这四行，无调用者传过它，因此在同一个任务里从两个签名中删掉。

**`voronoi.py` 两个 shell 改为 identity re-export（选项 a）。** `RadicalVoronoi.__call__` 转发 `self._inner.build(...)`（`voronoi.py:53`），`VoronoiIntegration.__call__` 转发 `self._inner.integrate(...)`（`:74`）——molrs 的属主用的动词是 `build` / `integrate`（`molrs/_lib.pyi:2667-2671`、`:2769-2783`，构造器都无参），不是 `compute`。如果只改名，这两个类就只剩「给 molrs 的动词起个别名」这一个存在理由，而后续的 re-export 重构会对同一批公开名字造成**第二次** breaking rename。所以现在就做掉：`voronoi.py` 退化为纯 re-export 模块（与 `pmsd.py` 同形），`RadicalVoronoi = molrs.compute.voronoi.RadicalVoronoi`、`VoronoiIntegration = molrs.compute.voronoi.VoronoiIntegration`，与该模块已有的 `VoronoiCells` / `voronoi_domains` / `voronoi_voids` 三行并列。行为逐字相同（原 shell 就是无参构造 + 单行转发），没有测试引用它们，文档改动落在同一个 docs 任务里，因此不占额外任务格。`docs/compute/voronoi.md:69-70` 现在写着「despite what its docstring says, you call the object itself; there is no `.compute` method」——这句道歉随之消失。

**放弃的东西，写在明处。** `base.py:51-52` 的 `@abstractmethod __call__` 提供了「子类没实现就不能实例化」的构造期护栏；换成 Protocol 后只有在显式 `isinstance(…, Compute)` 时做一次存在性检查，没有构造期护栏。今天全仓 `isinstance(…, Compute)` 出现 0 次，这层保护实际上没人用；experimental 阶段（铁律 3）接受这个交换，换来契约不在 molpy 侧被复述第二遍。

**测试树同构。** `tests/test_compute/test_orientations.py` 横跨 `pmft`/`neighborlist`/`order`/`spatial` 四个模块，`test_conductivity.py` 横跨 `pmsd`+`result`，都违反 `.claude/notes/testing.md:10-11` 的逐模块镜像。既然要改就拆：orientations → `test_order.py` / `test_spatial.py` / `test_pmft.py`（`NeighborList` 只是搭 nlist 的手段，其自身行为已由 `test_neighborlist.py` 覆盖）；conductivity → 只留 `test_result.py`。共享的 `orientations` 拓扑块构造器走 `tests/test_compute/conftest.py` 的 fixture，不做测试间 import（`testing.md:22`）。同时按 `testing.md:19-21`（纯 molrs re-export 模块没有 molpy 测试）删除三处 molrs 数值重复：`test_conductivity.py` 的 `test_einstein_conductivity_raw_msd`（`pmsd.py:8` 是纯 re-export）、`tests/test_compute/test_onsager.py`（`onsager.py:11`）、`tests/test_compute/test_persist.py`（`persist.py:10`）。**不**新建 `test_pmsd.py`。

**文档数据脚本同批改。** `.claude/notes/docs-style.md:134` 规定 `docs/compute/` 下每条曲线都来自 `scripts/docs_data/`；那里有 27 处以 `obj(...)` 调壳的站点（`structure.py:55,56,103,138,139,164`；`aggregate.py:39,40,96,97,100,101,135,136,153,155`；`order.py:45,50,88,89,110,111`；`dynamics.py:28,54,171`；`transport.py:35`；`angles.py:48`），pytest 与 `ty check src/molpy/` 都看不到它们，改名后全部 `TypeError`。它们随任务 9 一起改为 `.compute(...)`（`dynamics.py:54` 改为 `.build(...)`），ac-010 单独验收。

**动词表回填。** sub-spec 01 写入 `.claude/notes/architecture.md` § 设计铁律 4 的动词表把「分析」行标为「落地前为已声明的债（sub-spec 02）」并附退役条款；本规格落地时把该行改为「已成立」、删除债务清单第 1 行（任务 7、ac-011）。

**Reuse decision**（按 codebase 实证裁定）：

- `reuse molrs.compute.protocol.Compute` —— 契约直接复用，molpy 不新建任何基类/Protocol。
- `reuse molrs.compute.voronoi.RadicalVoronoi` / `VoronoiIntegration` —— 直接 identity re-export，不保留 molpy 类。
- `reuse` 现有 kernel 的 `.compute(...)` 入口 —— shell 新方法体与旧 `__call__` 体逐字相同。
- `reuse tests/test_compute/conftest.py` 的 `random_periodic_frame` / `frame_coords_snapshot` —— 三个新测试文件用它们建帧，不复制私有 helper。
- 无 `new` 符号：本规格不引入任何新公开符号。

## Files to create or modify

- `src/molpy/compute/base.py` —— 删除
- `src/molpy/compute/__init__.py`
- `src/molpy/compute/cluster.py`
- `src/molpy/compute/decomposition.py`
- `src/molpy/compute/density.py`
- `src/molpy/compute/dielectric.py`
- `src/molpy/compute/diffraction.py`
- `src/molpy/compute/distribution.py`
- `src/molpy/compute/environment.py`
- `src/molpy/compute/hbond.py`
- `src/molpy/compute/msd.py`
- `src/molpy/compute/neighborlist.py`
- `src/molpy/compute/order.py`
- `src/molpy/compute/pmft.py`
- `src/molpy/compute/rdf.py`
- `src/molpy/compute/reorientation.py`
- `src/molpy/compute/shape.py`
- `src/molpy/compute/spatial.py`
- `src/molpy/compute/van_hove.py`
- `src/molpy/compute/voronoi.py`
- `tests/test_compute/conftest.py`
- `tests/test_compute/test_rdf.py`
- `tests/test_compute/test_neighborlist.py`
- `tests/test_compute/test_diffraction.py`
- `tests/test_compute/test_distribution.py`
- `tests/test_compute/test_reorientation.py`
- `tests/test_compute/test_orientations.py` —— 删除
- `tests/test_compute/test_conductivity.py` —— 删除
- `tests/test_compute/test_onsager.py` —— 删除
- `tests/test_compute/test_persist.py` —— 删除
- `tests/test_compute/test_order.py` (new)
- `tests/test_compute/test_spatial.py` (new)
- `tests/test_compute/test_pmft.py` (new)
- `tests/test_compute/test_result.py` (new)
- `tests/test_compute/test_init.py` (new)
- `docs/api/compute.md`
- `docs/developer/extending-compute.md`
- `docs/developer/molrs-backend.md`
- `docs/compute/index.md`
- `docs/compute/cluster.md`
- `docs/compute/decomposition.md`
- `docs/compute/density.md`
- `docs/compute/diffraction.md`
- `docs/compute/distribution.md`
- `docs/compute/environment.md`
- `docs/compute/hbond.md`
- `docs/compute/msd.md`
- `docs/compute/neighborlist.md`
- `docs/compute/order.md`
- `docs/compute/pmft.md`
- `docs/compute/rdf.md`
- `docs/compute/reorientation.md`
- `docs/compute/shape.md`
- `docs/compute/spatial.md`
- `docs/compute/van_hove.md`
- `docs/compute/voronoi.md`
- `docs/index.md`
- `docs/getting-started/migration-0-14.md`
- `scripts/docs_data/structure.py`、`aggregate.py`、`order.py`、`dynamics.py`、`transport.py`、`angles.py` — 27 处 `obj(...)` 调用改为 `.compute(...)`（`dynamics.py:54` 的 `RadicalVoronoi()(...)` 改为 `.build(...)`）
- `.claude/notes/architecture.md` — 动词表「分析」行改为已成立、删除债务清单第 1 行

## Tasks

- [ ] Write failing unit tests: switch the verb to `shell.compute(...)` in `tests/test_compute/test_rdf.py`, `test_neighborlist.py`, `test_diffraction.py`, `test_distribution.py`, `test_reorientation.py`, and drop their `from molpy.compute.base import Compute` imports together with the `issubclass(..., Compute)` asserts the dying ABC guarded (inputs and expected values unchanged)
- [ ] Split `tests/test_compute/test_orientations.py` into `tests/test_compute/test_order.py`, `test_spatial.py`, `test_pmft.py` against the `compute` verb, move the shared `orientations`-block builder into `tests/test_compute/conftest.py` as a fixture, and delete `test_orientations.py`
- [ ] Prune the mirror per `.claude/notes/testing.md:19-21`: move `DielectricResult.fit_debye` into `tests/test_compute/test_result.py`, then delete `tests/test_compute/test_conductivity.py`, `test_onsager.py` and `test_persist.py` (all three assert molrs numerics through pure re-exports `pmsd.py:8`, `onsager.py:11`, `persist.py:10`)
- [ ] Write failing unit test `tests/test_compute/test_init.py` asserting `molpy.compute.Compute is molrs.compute.Compute`
- [ ] Rename `__call__` to `compute` and drop `super().__init__(**config)`, the `from .base import Compute` import and the `Compute` base declaration across the 31 shells in the 17 modules `cluster, decomposition, density, dielectric, diffraction, distribution, environment, hbond, msd, neighborlist, order, pmft, rdf, reorientation, shape, spatial, van_hove` (including the docstring examples at `msd.py:23`, `neighborlist.py:32`), and delete the now-sinkless `**config_kwargs` parameter at `dielectric.py:130` and `:462`
- [ ] Convert `src/molpy/compute/voronoi.py` into a pure re-export module binding `RadicalVoronoi` and `VoronoiIntegration` directly to `molrs.compute.voronoi`, deleting both molpy shell classes
- [ ] Delete `src/molpy/compute/base.py`, repoint `src/molpy/compute/__init__.py:19` to `from molrs.compute import Compute` (keeping `"Compute"` in `__all__`), replace the dangling `::: molpy.compute.base` target at `docs/api/compute.md:74` with `::: molpy.compute.Compute`, delete the stale `JACF` / `PMSDCompute` alias sentence and table at `docs/api/compute.md` (from "Historical tame names" on `:29` through the table ending `:35`, keeping the words "class." that close the `:28-29` sentence; neither symbol exists in `src/`; 铁律 3 leaves no alias), replace the false claim that no all-in-one recipe class exists with the true statement that `IonicConductivity` / `DielectricSusceptibility` currently exist, violate the raw-Compute → Fit → scale rule at `:20-26`, are routed to `/mol:refactor` for removal and must not be used in new code, rewrite the package docstring `src/molpy/compute/__init__.py:3-4` ("call it on frames" → `compute(...)`), and flip the 分析 row of the verb table in `.claude/notes/architecture.md` § 设计铁律 4 to 已成立 and delete debt row 1
- [ ] Rewrite the contract prose against `compute()` with no `Compute` base class and no `dump()`: `docs/developer/extending-compute.md` (`:12-20` table, the worked example, `:69`, `:72-85`), the "configure, then call" paragraph at `docs/compute/index.md:33-42`, and the Compute bullet at `docs/getting-started/migration-0-14.md:73-75` (replace it, do not append)
- [ ] Update the shell call sites to `obj.compute(...)` in `docs/compute/{index,cluster,decomposition,density,diffraction,distribution,environment,hbond,msd,neighborlist,order,pmft,rdf,reorientation,shape,spatial,van_hove}.md`, to `.build(...)` / `.integrate(...)` in `docs/compute/voronoi.md` (dropping the "you call the object itself" note at `:69-70`), in `docs/index.md:216-217` plus `docs/developer/molrs-backend.md:98-170` — where the "small `Compute` operator that forwards arguments" sentence at `:130-132` must not be re-emitted unqualified — and at the 27 call sites in `scripts/docs_data/{structure,aggregate,order,dynamics,transport,angles}.py` (`.compute(...)`; `dynamics.py:54` → `.build(...)`)
- [ ] Run full check + test suite

## Testing strategy

按 `.claude/notes/testing.md`，只写单元测试，路径逐模块镜像 `src/molpy/compute/`，一个测试一个行为，输入手写。

- **动词切换就是验证本身**：既有的逐 shell 测试保留原输入与原期望值，只把 `obj(...)` 换成 `obj.compute(...)`。测试能跑通即证明改名后的方法就是原来那个方法——这是「数值不动」的检验方式。
- **不加逐模块 `isinstance` 断言**：`molrs.compute.Compute` 是 presence-only 的 Protocol，`isinstance(shell, Compute)` 除了「有个 `compute` 属性」什么也不保证，而这一点已被上一条的真实调用证明。既有的 6 条 `issubclass(X, Compute)` 断言守的是「继承 molpy 那个 ABC」，ABC 按设计消失，断言随之删除——不是放宽断言，是被断言的结构不再存在。
- **唯一保留的契约断言**：`tests/test_compute/test_init.py`（镜像 `src/molpy/compute/__init__.py`）里一条 `molpy.compute.Compute is molrs.compute.Compute`，这是铁律 6 的门卫。测试树里 `Compute` 只在这一处被 import，且从 `molpy.compute` 导入，任何文件都不再穿透 `molpy.compute.base`。
- **Edge case 原样搬家**：`test_orientations.py` 三条 `pytest.raises(TypeError)` 的「拒绝外部 director/orientation 数组」契约随拆分进入 `test_order.py` / `test_pmft.py` / `test_spatial.py`，断言不变；`test_rdf.py` 的 free-box `ValueError` 断言不变；`test_conductivity.py` 的 `DielectricResult.fit_debye`（tau ≈ 5.0, rel=0.2）进入 `test_result.py`，数值不变。
- **删除而非搬迁的**：`test_einstein_conductivity_raw_msd`、`test_onsager.py`、`test_persist.py` —— 它们测的是 molrs 的数值，经由 molpy 的纯 re-export 模块（`testing.md:19-21`）。不以 identity 断言替代：一条 `is` 断言对纯 re-export 模块没有独立价值，铁律 6 的门卫由 `test_init.py` 一处承担。
- **`voronoi.py` 无测试**：转为 re-export 后按 `testing.md:19-21` 不建 molpy 测试；其 identity 形状由静态验收项 ac-005 覆盖。
- **不做的事**（`testing.md:10-15` 与本链裁定）：不建 `regressions/` 目录；不做源码文本断言（"compute/ 下不再有 `def __call__`" 是静态审查的 `type: code` 验收项，不是测试）；不做遍历全部 shell 的 import-all 参数化门；不做 molrs 数值 parity；不引入第三方科学软件。
- **单点绿灯**：`uv run --extra dev python -m pytest tests/test_compute/test_<mod>.py`；全量门 `uv run --extra dev python -m pytest tests/ -n auto`。

## Out of scope

- **`core/selector.py` 的 `MaskPredicate`**：本链唯一的 `__call__` 豁免（对象本身就是函数），`docs/tutorials/06_selector.md:40-58` 不动。
- **molrs 拥有的动词**：`Onsager.correlation`、`Persist.pair_survival_tcf`、`LinearFit.fit`、`CumulativeTrapezoid`、`DebyeFit`、`RadicalVoronoi.build`、`VoronoiIntegration.integrate` 保持原样，它们是 sub-spec 01 动词表里 `molrs:` 的描述性行，不是 molpy 可以改的契约。
- **其他家族的动词统一**（`io` / `builder` / `typifier` / `wrapper`）：本链其他 sub-spec。
- **已发现、未修、已路由的债（铁律：不静默留存）**：
  1. **all-in-one 门面类（违规的是代码，不是文档）**：`src/molpy/compute/dielectric.py:78 DielectricSusceptibility` 与 `:418 IonicConductivity` 在一次调用里做完 unwrap → 收集偶极 → `EinsteinConductivity().compute` → `LinearFit().fit` → SI 前因子（`:480-560`，前因子在 `:538-546`），正是 CLAUDE.md § Forbid「All-in-one façade APIs」禁止的形状，也与 `.claude/notes/architecture.md:49`、`docs/api/compute.md:28`、`docs/developer/extending-compute.md:5-10` 所写的「raw Compute → Fit → 由调用者做 SI 缩放」直接冲突。**冲突的解法是删类，不是删规则**——本规格把这两个类的 `__call__` 与其他 shell 一样改名为 `compute`，但不因此承认它们合法。同一个类还有第二个分析入口：`dielectric.py:155 from_dipole_series(M, …)` 是**实例方法**、返回 `DielectricSusceptibilityResult`（`:161`，经 `:179-185` 的 `_spectra_from_series`），`from_*` 的名字是误用——任务 5 之后 `DielectricSusceptibility` 仍暴露 `compute(trajectory)` 与 `from_dipole_series(M, …)` 两个入口，ac-002 只数 `def compute(`，看不见它。**路由**：单独的 `/mol:refactor`（拆成 primitive + 文档里的组合示例，把 `from_dipole_series` 折进 primitive）；**退役条款**：动词表「分析」行在本规格落地时改为「已成立」（ac-011）的含义是「molpy 拥有的分析入口动词是 `compute`」，`from_dipole_series` 作为该行下的已声明残余债写进债务清单，由那次重构删除。在那之前 `docs/api/compute.md:28` 按任务 7 改为陈述现状（两类存在、违规、已路由）。
  2. **26 个纯转发壳**（持 `_impl`/`_inner`、方法体只有一行 `return self._impl.<verb>(...)`，零增补，违反架构铁律 6「转发门面：永不」，正确形状是 identity re-export）。以今日树为准的行号：`cluster.py:57,:72`；`decomposition.py:30,:51`；`density.py:36,:55`；`diffraction.py:41`；`distribution.py:43,:66,:85,:110`；`environment.py:36`；`hbond.py:52`；`msd.py:33`；`order.py:54,:71,:89,:110`；`pmft.py:41`；`reorientation.py:55`；`shape.py:57,:68,:85,:102`；`spatial.py:76`；`van_hove.py:53`。（同族共 28 处，其中 `voronoi.py:53,:74` 由本规格当场修掉，见 Design；`rdf.py:52` 带 box 校验与 list 归一、`cluster.py:45-46` 带 `keys` 分支、`neighborlist.py:41-48` 与 `dielectric.py` 带真实行为，均不属此列。）**门面上已有分裂**：`src/molpy/__init__.py:223,227-231,235,245,249` 把 molrs 的类以**同名**挂在包根，而 `molpy.compute.<同名>` 是 molpy 的壳——运行时 `molpy.VanHove is molpy.compute.VanHove` 为 `False`，同样分裂的还有 `SpatialDistribution`、`AngleDistribution`、`CombinedDistribution`、`DihedralDistribution`、`DistanceDistribution`、`HBonds`、`LegendreReorientation`（`RadicalVoronoi`/`VoronoiIntegration` 的同一分裂由本规格的选项 (a) 修掉；`HBondCriterion` 已是同一对象）。本规格之后这八对名字动词相同、对象不同——重构的目标因此是 **identity 收敛**（`molpy.X is molpy.compute.X is molrs 类`），不只是动词卫生。**无镜像测试的模块**：`cluster.py`（3 壳）、`decomposition.py`（2）、`density.py`（2）、`environment.py`（1）、`hbond.py`（1）、`msd.py`（1）、`shape.py`（4）、`van_hove.py`（1）——15 个 molpy 拥有的类今天零测试，不在 `testing.md:19-21` 的豁免内；它们的改名只能由 ac-002/ac-003 的静态审查与 `scripts/docs_data/` 的调用（ac-010）核验。重构把它们变成 re-export 后 `testing.md:19-21` 才适用。**路由**：单独 sub-spec / `/mol:refactor`——它改的是公开符号的归属而非动词，且会把本规格推过 10 任务上限。对这 26 个名字，重构将是第二次 breaking 改名，铁律 3（experimental）接受。
  3. **`tests/test_compute/test_dielectric.py` 的覆盖错位**：该文件只测 molrs 的 `Dielectric` re-export（`compute/__init__.py:26`），属 `testing.md:19-21` 的重复；而它镜像的 `dielectric.py` 里真正属于 molpy 的两个类今天**零测试**。本规格不动它——删掉会让一个有真实行为的模块失去镜像文件。**路由**：与债 1 的 `/mol:refactor` 同行，那次重构必须用这两个 molpy 类的单元测试**替换**该文件的内容（不是追加），否则错位的 molrs 数值重复会一并存活。
  4. **CLAUDE.md § Prefer「OOP by default」以 `NeighborList.build` 作范例**（`CLAUDE.md:71`），而本规格之后 `molpy.compute.NeighborList` 带的是 `compute`。该行位于 bootstrap 托管块内（`CLAUDE.md:25-133`，范例在 `:71`），不在本规格手改。**路由**：操作者重跑 `/mol:bootstrap`（sub-spec 01 的债务清单第 8 行）。
- **`compute` 族的 API 语义**：参数名、返回类型、kernel 选择、单位一律不动。
