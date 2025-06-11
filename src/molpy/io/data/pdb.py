from pathlib import Path

from .base import DataWriter, DataReader
import numpy as np

import molpy as mp
from collections import defaultdict
import xarray as xr


class PDBReader(DataReader):

    def __init__(self, file: Path):
        super().__init__(path=file)

    @staticmethod
    def sanitizer(line: str) -> str:
        return line.strip()

    def read(self, frame):

        with open(self._path, "r") as f:

            lines = filter(
                lambda line: line.startswith("ATOM")
                or line.startswith("CONECT")
                or line.startswith("HETATM"),
                map(PDBReader.sanitizer, f),
            )

            atoms = {
                "id": [],
                "name": [],
                "resName": [],
                "chainID": [],
                "resSeq": [],
                "xyz": [],
                "element": [],
            }

            bonds = []

            for line in lines:

                if line.startswith("ATOM") or line.startswith("HETATM"):

                    atoms["id"].append(int(line[6:11]))
                    atoms["name"].append(line[12:16].strip())
                    atoms["resName"].append(line[17:20].strip())
                    atoms["chainID"].append(line[21])
                    atoms["resSeq"].append(int(line[22:26]))
                    atoms["xyz"].append(
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    )
                    atoms["element"].append(line[76:78].strip())

                elif line.startswith("CONECT"):

                    bond_indices = list(map(lambda l: int(l), line.split()[1:]))
                    i = bond_indices[0]
                    for j in bond_indices[1:]:
                        bonds.append([i, j] if i < j else [j, i])

            bonds = np.unique(np.array(bonds), axis=0)

            if len(set(atoms["name"])) != len(atoms["name"]):
                atom_name_counter = defaultdict(int)
                for i, name in enumerate(atoms["name"]):
                    atom_name_counter[name] += 1
                    if atom_name_counter[name] > 1:
                        atoms["name"][i] = f"{name}{atom_name_counter[name]}"

            atoms_ds = {k: ("index", np.asarray(v)) for k, v in atoms.items()}
            frame["atoms"] = xr.Dataset(atoms_ds)
            frame["box"] = mp.Box()

            if len(bonds):
                frame["bonds"] = xr.Dataset(
                    {
                        "i": ("index", bonds[:, 0]),
                        "j": ("index", bonds[:, 1]),
                    }
                )


            return frame


class PDBWriter(DataWriter):

    def __init__(
        self,
        path: Path,
    ):
        super().__init__(path=path)

    def write(self, frame):

        with open(self._path, "w") as f:
            atoms = frame["atoms"]
            name = frame.get("name", "MOL")
            f.write(f"REMARK  {name}\n")
            n_atoms = atoms.dims.get("index", 0)
            for i in range(n_atoms):
                serial = i + 1
                altLoc = atoms.get("altLoc", xr.DataArray("")).values[i]
                unique_name = atoms.get("unique_name", atoms["name"]).values[i]

                resName = atoms.get("resName", xr.DataArray("UNK")).values[i]
                chainID = atoms.get("chainID", xr.DataArray("A")).values[i]
                resSeq = atoms.get("resSeq", xr.DataArray(1)).values[i]
                iCode = atoms.get("iCode", xr.DataArray("")).values[i]
                elem = atoms.get("element", xr.DataArray("X")).values[i]
                x, y, z = atoms["xyz"].values[i]

                f.write(
                    f"{'HETATM':6s}{serial:>5d} {unique_name.upper():>4s}{altLoc:1s}{resName:>3s} {chainID:1s}{resSeq:>4d}{iCode:1s}   "
                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{' '*22}{elem:>2s}  \n"
                )

            bonds = defaultdict(list)
            if "bonds" in frame:
                bds = frame["bonds"]
                n_bonds = bds.dims.get("index", 0)
                for i in range(n_bonds):
                    bi = int(bds["i"].values[i] + 1)
                    bj = int(bds["j"].values[i] + 1)
                    bonds[bi].append(bj)
                    bonds[bj].append(bi)

                for i, js in bonds.items():
                    if len(js) > 4:
                        raise ValueError(
                            f"PDB only supports up to 4 bonds, but atom {i} has {len(js)} bonds: {js}"
                        )
                    f.write(
                        f"{'CONECT':6s}{i:>5d}{''.join([f'{j:>5d}' for j in js])}\n"
                    )

            f.write("END\n")
