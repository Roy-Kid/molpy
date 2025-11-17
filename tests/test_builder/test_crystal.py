"""
测试晶格构建器模块
"""

from __future__ import annotations

import numpy as np
import pytest

from molpy.builder.crystal import (
    Site,
    Lattice,
    BlockRegion,
    CrystalBuilder,
)
from molpy import Atomistic


class TestSite:
    """测试 Site 数据类"""

    def test_site_creation(self):
        """测试创建基位点"""
        site = Site(label="A", species="Ni", frac=(0.0, 0.0, 0.0))
        assert site.label == "A"
        assert site.species == "Ni"
        assert site.frac == (0.0, 0.0, 0.0)
        assert site.charge == 0.0
        assert site.meta is None

    def test_site_with_charge(self):
        """测试带电荷的位点"""
        site = Site(label="Na", species="Na", frac=(0.0, 0.0, 0.0), charge=1.0)
        assert site.charge == 1.0

    def test_site_with_metadata(self):
        """测试带元数据的位点"""
        meta = {"tag": "corner"}
        site = Site(label="A", species="C", frac=(0.0, 0.0, 0.0), meta=meta)
        assert site.meta == meta


class TestLattice:
    """测试 Lattice 类"""

    def test_lattice_creation(self):
        """测试创建晶格"""
        a1 = np.array([1.0, 0.0, 0.0])
        a2 = np.array([0.0, 1.0, 0.0])
        a3 = np.array([0.0, 0.0, 1.0])
        lattice = Lattice(a1=a1, a2=a2, a3=a3, basis=[])

        assert np.allclose(lattice.a1, a1)
        assert np.allclose(lattice.a2, a2)
        assert np.allclose(lattice.a3, a3)
        assert len(lattice.basis) == 0

    def test_lattice_cell_property(self):
        """测试 cell 属性返回正确的矩阵"""
        a1 = np.array([3.0, 0.0, 0.0])
        a2 = np.array([0.0, 3.0, 0.0])
        a3 = np.array([0.0, 0.0, 3.0])
        lattice = Lattice(a1=a1, a2=a2, a3=a3, basis=[])

        cell = lattice.cell
        assert cell.shape == (3, 3)
        assert np.allclose(cell[0], a1)
        assert np.allclose(cell[1], a2)
        assert np.allclose(cell[2], a3)

    def test_add_site(self):
        """测试添加基位点"""
        lattice = Lattice(
            a1=np.array([1.0, 0.0, 0.0]),
            a2=np.array([0.0, 1.0, 0.0]),
            a3=np.array([0.0, 0.0, 1.0]),
            basis=[],
        )
        site = Site(label="A", species="C", frac=(0.0, 0.0, 0.0))
        lattice.add_site(site)
        assert len(lattice.basis) == 1
        assert lattice.basis[0] == site

    def test_frac_to_cart_single(self):
        """测试分数坐标到笛卡尔坐标转换（单个点）"""
        a1 = np.array([3.0, 0.0, 0.0])
        a2 = np.array([0.0, 4.0, 0.0])
        a3 = np.array([0.0, 0.0, 5.0])
        lattice = Lattice(a1=a1, a2=a2, a3=a3, basis=[])

        frac = np.array([0.5, 0.5, 0.5])
        cart = lattice.frac_to_cart(frac)
        expected = np.array([1.5, 2.0, 2.5])
        assert np.allclose(cart, expected)

    def test_frac_to_cart_multiple(self):
        """测试分数坐标到笛卡尔坐标转换（多个点）"""
        a1 = np.array([2.0, 0.0, 0.0])
        a2 = np.array([0.0, 2.0, 0.0])
        a3 = np.array([0.0, 0.0, 2.0])
        lattice = Lattice(a1=a1, a2=a2, a3=a3, basis=[])

        frac = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.5, 0.5, 0.5]]
        )
        cart = lattice.frac_to_cart(frac)
        expected = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [1.0, 1.0, 1.0]]
        )
        assert np.allclose(cart, expected)

    def test_cubic_sc(self):
        """测试简单立方晶格"""
        lat = Lattice.cubic_sc(a=2.0, species="Cu")

        # 检查晶格向量
        assert np.allclose(lat.a1, [2.0, 0.0, 0.0])
        assert np.allclose(lat.a2, [0.0, 2.0, 0.0])
        assert np.allclose(lat.a3, [0.0, 0.0, 2.0])

        # 检查基位点
        assert len(lat.basis) == 1
        assert lat.basis[0].species == "Cu"
        assert lat.basis[0].frac == (0.0, 0.0, 0.0)

    def test_cubic_bcc(self):
        """测试体心立方晶格"""
        lat = Lattice.cubic_bcc(a=3.0, species="Fe")

        # 检查晶格向量
        assert np.allclose(lat.a1, [3.0, 0.0, 0.0])

        # 检查基位点数量
        assert len(lat.basis) == 2

        # 检查位点位置
        fracs = [site.frac for site in lat.basis]
        assert (0.0, 0.0, 0.0) in fracs
        assert (0.5, 0.5, 0.5) in fracs

        # 检查物种
        for site in lat.basis:
            assert site.species == "Fe"

    def test_cubic_fcc(self):
        """测试面心立方晶格"""
        lat = Lattice.cubic_fcc(a=3.52, species="Ni")

        # 检查基位点数量
        assert len(lat.basis) == 4

        # 检查位点位置
        fracs = [site.frac for site in lat.basis]
        assert (0.0, 0.0, 0.0) in fracs
        assert (0.5, 0.5, 0.0) in fracs
        assert (0.5, 0.0, 0.5) in fracs
        assert (0.0, 0.5, 0.5) in fracs

        # 检查物种
        for site in lat.basis:
            assert site.species == "Ni"

    def test_rocksalt(self):
        """测试岩盐结构"""
        lat = Lattice.rocksalt(a=5.64, species_a="Na", species_b="Cl")

        # 检查基位点数量：4个Na + 4个Cl
        assert len(lat.basis) == 8

        # 计数每种物种
        na_count = sum(1 for site in lat.basis if site.species == "Na")
        cl_count = sum(1 for site in lat.basis if site.species == "Cl")
        assert na_count == 4
        assert cl_count == 4

        # 检查Na的位置（fcc位点）
        na_fracs = [site.frac for site in lat.basis if site.species == "Na"]
        assert (0.0, 0.0, 0.0) in na_fracs
        assert (0.5, 0.5, 0.0) in na_fracs

        # 检查Cl的位置（偏移的fcc位点）
        cl_fracs = [site.frac for site in lat.basis if site.species == "Cl"]
        assert (0.5, 0.0, 0.0) in cl_fracs
        assert (0.0, 0.5, 0.0) in cl_fracs


class TestBlockRegion:
    """测试 BlockRegion 类"""

    def test_lattice_coord_system(self):
        """测试晶格坐标系统"""
        region = BlockRegion(0, 10, 0, 10, 0, 10, coord_system="lattice")

        points = np.array(
            [
                [5.0, 5.0, 5.0],  # 内部
                [0.0, 0.0, 0.0],  # 边界（包含）
                [10.0, 10.0, 10.0],  # 边界（包含）
                [11.0, 5.0, 5.0],  # 外部
                [-1.0, 5.0, 5.0],  # 外部
            ]
        )

        mask = region.contains_mask(points)
        assert mask[0] == True
        assert mask[1] == True
        assert mask[2] == True
        assert mask[3] == False
        assert mask[4] == False

    def test_cartesian_coord_system(self):
        """测试笛卡尔坐标系统"""
        region = BlockRegion(0, 10, 0, 10, 0, 10, coord_system="cartesian")

        points = np.array(
            [
                [5.0, 5.0, 5.0],
                [0.0, 0.0, 0.0],
                [10.0, 10.0, 10.0],
                [11.0, 5.0, 5.0],
            ]
        )

        mask = region.contains_mask(points)
        assert mask[0] == True
        assert mask[1] == True
        assert mask[2] == True
        assert mask[3] == False

    def test_partial_overlap(self):
        """测试部分重叠的点"""
        region = BlockRegion(-5, 5, -5, 5, -5, 5, coord_system="lattice")

        points = np.array(
            [
                [0.0, 0.0, 0.0],  # 内部
                [5.0, 0.0, 0.0],  # x边界
                [0.0, 6.0, 0.0],  # y外部
                [-5.0, -5.0, -5.0],  # 角落（包含）
            ]
        )

        mask = region.contains_mask(points)
        assert mask[0] == True
        assert mask[1] == True
        assert mask[2] == False
        assert mask[3] == True


class TestCrystalBuilder:
    """测试 CrystalBuilder 类"""

    def test_simple_cubic_small(self):
        """测试生成小的简单立方结构"""
        lat = Lattice.cubic_sc(a=2.0, species="Cu")
        region = BlockRegion(0, 2, 0, 2, 0, 2, coord_system="lattice")
        builder = CrystalBuilder(lat)

        structure = builder.build_block(region)

        # 检查返回类型
        assert isinstance(structure, Atomistic)

        # 检查原子数量：2x2x2 = 8 个原子
        assert len(list(structure.atoms)) == 8

        # 检查物种
        symbols = structure.symbols
        assert all(s == "Cu" for s in symbols)

    def test_simple_cubic_with_explicit_ranges(self):
        """测试使用显式范围生成简单立方"""
        lat = Lattice.cubic_sc(a=3.0, species="Fe")
        region = BlockRegion(0, 5, 0, 5, 0, 5, coord_system="lattice")
        builder = CrystalBuilder(lat)

        # 使用显式范围
        structure = builder.build_block(
            region, i_range=range(0, 2), j_range=range(0, 2), k_range=range(0, 2)
        )

        # 2x2x2 = 8 个原子
        assert len(list(structure.atoms)) == 8

        # 检查坐标
        positions = structure.xyz
        assert positions.shape == (8, 3)

        # 检查第一个原子在原点
        assert np.allclose(positions[0], [0.0, 0.0, 0.0])

    def test_bcc_structure(self):
        """测试体心立方结构"""
        lat = Lattice.cubic_bcc(a=2.0, species="Fe")
        region = BlockRegion(0, 2, 0, 2, 0, 2, coord_system="lattice")
        builder = CrystalBuilder(lat)

        structure = builder.build_block(region)

        # 2x2x2 晶胞 × 2个基位点 = 16 个原子
        assert len(list(structure.atoms)) == 16

        symbols = structure.symbols
        assert all(s == "Fe" for s in symbols)

    def test_fcc_structure(self):
        """测试面心立方结构"""
        lat = Lattice.cubic_fcc(a=3.52, species="Ni")
        region = BlockRegion(0, 2, 0, 2, 0, 2, coord_system="lattice")
        builder = CrystalBuilder(lat)

        structure = builder.build_block(region)

        # 2x2x2 晶胞 × 4个基位点 = 32 个原子
        assert len(list(structure.atoms)) == 32

        symbols = structure.symbols
        assert all(s == "Ni" for s in symbols)

    def test_rocksalt_structure(self):
        """测试岩盐结构"""
        lat = Lattice.rocksalt(a=5.64, species_a="Na", species_b="Cl")
        region = BlockRegion(0, 2, 0, 2, 0, 2, coord_system="lattice")
        builder = CrystalBuilder(lat)

        structure = builder.build_block(region)

        # 2x2x2 晶胞 × 8个基位点 = 64 个原子
        assert len(list(structure.atoms)) == 64

        # 检查物种分布
        symbols = structure.symbols
        na_count = sum(1 for s in symbols if s == "Na")
        cl_count = sum(1 for s in symbols if s == "Cl")
        assert na_count == 32
        assert cl_count == 32

    def test_empty_basis(self):
        """测试空基的情况"""
        lat = Lattice(
            a1=np.array([1.0, 0.0, 0.0]),
            a2=np.array([0.0, 1.0, 0.0]),
            a3=np.array([0.0, 0.0, 1.0]),
            basis=[],
        )
        region = BlockRegion(0, 2, 0, 2, 0, 2, coord_system="lattice")
        builder = CrystalBuilder(lat)

        structure = builder.build_block(region)

        # 空基应该返回空结构
        assert len(list(structure.atoms)) == 0
        assert isinstance(structure, Atomistic)

    def test_cartesian_region_without_ranges_raises_error(self):
        """测试笛卡尔坐标系统不提供范围时抛出错误"""
        lat = Lattice.cubic_sc(a=2.0, species="Cu")
        region = BlockRegion(0, 10, 0, 10, 0, 10, coord_system="cartesian")
        builder = CrystalBuilder(lat)

        with pytest.raises(ValueError, match="i_range.*j_range.*k_range"):
            builder.build_block(region)

    def test_cartesian_region_with_ranges(self):
        """测试笛卡尔坐标系统提供显式范围"""
        lat = Lattice.cubic_sc(a=2.0, species="Cu")
        region = BlockRegion(0, 3, 0, 3, 0, 3, coord_system="cartesian")
        builder = CrystalBuilder(lat)

        # 显式提供范围
        structure = builder.build_block(
            region, i_range=range(0, 2), j_range=range(0, 2), k_range=range(0, 2)
        )

        # 应该创建 2x2x2 = 8 个原子
        # 但只有在 (0,3) x (0,3) x (0,3) 范围内的才会被保留
        # 由于晶格常数为2，所以 (0,0,0), (2,0,0), (0,2,0), (0,0,2), (2,2,0), (2,0,2), (0,2,2), (2,2,2)
        # 所有这些都在笛卡尔区域内
        assert len(list(structure.atoms)) == 8

    def test_positions_are_correct(self):
        """测试生成的原子位置正确"""
        lat = Lattice.cubic_sc(a=2.0, species="Cu")
        region = BlockRegion(0, 2, 0, 2, 0, 2, coord_system="lattice")
        builder = CrystalBuilder(lat)

        structure = builder.build_block(region)
        positions = structure.xyz

        # 期望的位置
        expected_positions = np.array(
            [
                [0, 0, 0],
                [2, 0, 0],
                [0, 2, 0],
                [0, 0, 2],
                [2, 2, 0],
                [2, 0, 2],
                [0, 2, 2],
                [2, 2, 2],
            ],
            dtype=float,
        )

        # 排序以便比较
        positions_sorted = positions[
            np.lexsort((positions[:, 2], positions[:, 1], positions[:, 0]))
        ]
        expected_sorted = expected_positions[
            np.lexsort(
                (
                    expected_positions[:, 2],
                    expected_positions[:, 1],
                    expected_positions[:, 0],
                )
            )
        ]

        assert np.allclose(positions_sorted, expected_sorted)

    def test_box_is_set_correctly(self):
        """测试生成的结构的盒子信息正确"""
        lat = Lattice.cubic_sc(a=3.0, species="Cu")
        region = BlockRegion(0, 2, 0, 2, 0, 2, coord_system="lattice")
        builder = CrystalBuilder(lat)

        structure = builder.build_block(region)

        # 检查是否设置了 box
        assert "box" in structure
        from molpy import Box

        box = structure["box"]
        assert isinstance(box, Box)

        # 检查盒子矩阵
        expected_cell = lat.cell
        assert np.allclose(box.matrix, expected_cell)
