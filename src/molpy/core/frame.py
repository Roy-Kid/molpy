from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from typing import Any, Sequence, Literal

from xarray import DataTree
from .arraydict import ArrayDict  # legacy helper
from .utils import TagApplyer
from .box import Box


def _dataset_from_dicts(dicts: list[dict[str, Any]]) -> xr.Dataset:
    if not dicts:
        return xr.Dataset()
    keys = set().union(*[d.keys() for d in dicts])
    data = {k: ("index", np.asarray([d.get(k) for d in dicts])) for k in keys}
    return xr.Dataset(data)


class Frame(MutableMapping[str, xr.Dataset]):
    """Data container backed by ``xarray`` and ``DataTree``."""

    box: Box | None = None

    def __new__(cls, data: Mapping[str, Any] | None = None, *, style: str = "atomic") -> "Frame":
        if cls is Frame and style == "atomic":
            return super().__new__(AllAtomFrame)
        return super().__new__(cls)

    def __init__(self, data: Mapping[str, Any] | None = None, *args, **kwargs) -> None:
        data = data or {}
        datasets: dict[str, xr.Dataset] = {}
        for key, value in data.items():
            datasets[key] = self._to_dataset(value)
        self.tree = DataTree(datasets)

    # ------------------------------------------------------------------
    # Mapping interface
    # ------------------------------------------------------------------
    def _to_dataset(self, value: Any) -> xr.Dataset:
        if isinstance(value, xr.Dataset):
            return value
        if isinstance(value, ArrayDict):
            value = value.to_dict()
        if isinstance(value, Mapping):
            arrays = {k: ("index", np.asarray(v)) for k, v in value.items()}
            return xr.Dataset(arrays)
        raise TypeError("Frame values must be mapping or xarray.Dataset")

    def __getitem__(self, key: str) -> xr.Dataset:
        return self.tree[key].ds

    def __setitem__(self, key: str, value: Any) -> None:
        self.tree[key] = self._to_dataset(value)

    def __delitem__(self, key: str) -> None:
        del self.tree[key]

    def __iter__(self):
        return iter(self.tree.keys())

    def __len__(self) -> int:
        if "atoms" in self.tree:
            return self.tree["atoms"].ds.sizes.get("index", 0)
        return len(self.tree.keys())

    # Additional helpers ------------------------------------------------
    def keys(self):
        return self.tree.keys()

    def values(self):
        return (self.tree[k].ds for k in self.tree.keys())

    def items(self):
        for k in self.tree.keys():
            yield k, self.tree[k].ds

    # ------------------------------------------------------------------
    def copy(self) -> "Frame":
        copied = {k: v.copy() for k, v in self.items()}
        new = self.__class__(copied)
        new.box = None if self.box is None else Box.from_box(self.box)
        return new

    # ------------------------------------------------------------------
    @classmethod
    def from_frames(cls, others: Sequence["Frame"]) -> "Frame":
        base = others[0].copy()
        for other in others[1:]:
            for key in other.keys():
                if key in base:
                    base[key] = xr.concat([base[key], other[key]], dim="index")
                else:
                    base[key] = other[key].copy()
        return base

    @classmethod
    def from_structs(cls, structs):
        frame = cls()

        atom_dicts: list[dict[str, Any]] = []
        bond_dicts: list[dict[str, Any]] = []
        bond_index = []
        tager = TagApplyer()
        for struct in structs:
            if "bonds" in struct:
                topo = struct.get_topology()
                bond_dicts.extend([bond.to_dict() for bond in struct.bonds])
                bond_index.append(topo.bonds + len(atom_dicts))

            atom_dicts.extend([atom.to_dict() for atom in struct.atoms])
            tager.update_dollar_counter()
            tager.apply_tags(atom_dicts)
            tager.apply_tags(bond_dicts)

        frame["atoms"] = _dataset_from_dicts(atom_dicts)
        frame["bonds"] = _dataset_from_dicts(bond_dicts)
        if bond_dicts:
            bond_index = np.concatenate(bond_index)
            bds = frame["bonds"].copy()
            bds["i"] = ("index", bond_index[:, 0])
            bds["j"] = ("index", bond_index[:, 1])
            frame["bonds"] = bds
        return frame


class AllAtomMixin(MutableMapping[Literal["atoms", "bonds", "angles", "dihedrals", "impropers"], xr.Dataset]):
    """Mixin providing atomistic helpers for :class:`Frame`."""

    def split(self, masks: Sequence[bool] | Sequence[int] | np.ndarray) -> list["Frame"]:
        frames: list[Frame] = []
        masks = np.asarray(masks)
        if masks.dtype == bool:
            unique_masks = [masks]
        else:
            unique_masks = [masks == i for i in np.unique(masks)]

        for mask in unique_masks:
            frame = self.__class__()
            frame["atoms"] = self["atoms"].isel(index=mask)
            atom_ids = frame["atoms"].get("id")
            if atom_ids is not None:
                atom_ids = atom_ids.values
            if "bonds" in self:
                bds = self["bonds"]
                bi = bds.get("i")
                bj = bds.get("j")
                if bi is not None and bj is not None and atom_ids is not None:
                    bond_mask = np.logical_and(np.isin(bi, atom_ids), np.isin(bj, atom_ids))
                    frame["bonds"] = bds.isel(index=bond_mask)
            if "angles" in self:
                ads = self["angles"]
                ai = ads.get("i")
                aj = ads.get("j")
                ak = ads.get("k")
                if ai is not None and aj is not None and ak is not None and atom_ids is not None:
                    m = np.isin(ai, atom_ids) & np.isin(aj, atom_ids) & np.isin(ak, atom_ids)
                    frame["angles"] = ads.isel(index=m)
            if "dihedrals" in self:
                dds = self["dihedrals"]
                ai = dds.get("i")
                aj = dds.get("j")
                ak = dds.get("k")
                al = dds.get("l")
                if ai is not None and aj is not None and ak is not None and al is not None and atom_ids is not None:
                    m = (
                        np.isin(ai, atom_ids)
                        & np.isin(aj, atom_ids)
                        & np.isin(ak, atom_ids)
                        & np.isin(al, atom_ids)
                    )
                    frame["dihedrals"] = dds.isel(index=m)
            if "impropers" in self:
                ids = self["impropers"]
                ai = ids.get("i")
                aj = ids.get("j")
                ak = ids.get("k")
                al = ids.get("l")
                if ai is not None and aj is not None and ak is not None and al is not None and atom_ids is not None:
                    m = (
                        np.isin(ai, atom_ids)
                        & np.isin(aj, atom_ids)
                        & np.isin(ak, atom_ids)
                        & np.isin(al, atom_ids)
                    )
                    frame["impropers"] = ids.isel(index=m)
            frames.append(frame)
        return frames

    def to_struct(self):
        from .struct import Entities, Struct
        import molpy as mp

        struct = Struct()
        atoms = self["atoms"]
        n_atoms = atoms.dims.get("index", 0)
        for i in range(n_atoms):
            atom_dict = {k: atoms[k].values[i] for k in atoms.data_vars}
            struct.def_atom(**atom_dict)

        if "bonds" in self:
            struct["bonds"] = Entities()
            bonds = self["bonds"]
            for i in range(bonds.dims.get("index", 0)):
                bond = {k: bonds[k].values[i] for k in bonds.data_vars}
                bi = bond.pop("i")
                bj = bond.pop("j")
                itom = struct["atoms"].get_by(lambda a: a["id"] == bi)
                jtom = struct["atoms"].get_by(lambda a: a["id"] == bj)
                struct["bonds"].add(mp.Bond(itom, jtom, **bond))

        if "angles" in self:
            struct["angles"] = Entities()
            ads = self["angles"]
            for i in range(ads.dims.get("index", 0)):
                angle = {k: ads[k].values[i] for k in ads.data_vars}
                ai, aj, ak = angle.pop("i"), angle.pop("j"), angle.pop("k")
                itom = struct["atoms"].get_by(lambda a: a["id"] == ai)
                jtom = struct["atoms"].get_by(lambda a: a["id"] == aj)
                ktom = struct["atoms"].get_by(lambda a: a["id"] == ak)
                struct["angles"].add(mp.Angle(itom, jtom, ktom, **angle))

        if "dihedrals" in self:
            struct["dihedrals"] = Entities()
            dds = self["dihedrals"]
            for i in range(dds.dims.get("index", 0)):
                dihedral = {k: dds[k].values[i] for k in dds.data_vars}
                ai, aj, ak, al = (
                    dihedral.pop("i"),
                    dihedral.pop("j"),
                    dihedral.pop("k"),
                    dihedral.pop("l"),
                )
                itom = struct["atoms"].get_by(lambda a: a["id"] == ai)
                jtom = struct["atoms"].get_by(lambda a: a["id"] == aj)
                ktom = struct["atoms"].get_by(lambda a: a["id"] == ak)
                ltom = struct["atoms"].get_by(lambda a: a["id"] == al)
                struct["dihedrals"].add(
                    mp.Dihedral(itom, jtom, ktom, ltom, **dihedral)
                )

        return struct


class AllAtomFrame(Frame, AllAtomMixin):
    """Frame specialised for atomistic systems."""

    def __init__(self, data: Mapping[str, Any] | None = None, *args, **kwargs) -> None:
        data = data or {}
        for key in ["atoms", "bonds", "angles", "dihedrals", "impropers"]:
            data.setdefault(key, xr.Dataset())
        super().__init__(data, *args, **kwargs)

