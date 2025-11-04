from typing import Iterable, Optional

from molpy.core.aa.base import Atomistic
from molpy.core.aa import Atom
from molpy.core.aa import Bond
from molpy.core.assembly import Assembly
from molpy.core.entity import Entity
from molpy.core.wrappers.monomer import Monomer


class RDKit:
	
    def __init__(self):
        self._require_rdkit()

    def _require_rdkit(self):
        try:
            import rdkit
        except ImportError as e:
            raise ImportError(
                "RDKit is required for this functionality. Please install RDKit."
            ) from e
		
    @staticmethod
    def atomistic_to_rdkit(atomistic: Atomistic):
        """Convert an Assembly to an RDKit Chem.Mol.

        Notes
        -----
        - Uses AA Atom/Bond buckets when registered.
        - If 3D positions are present under key "pos", they are embedded as a conformer.
        - Only basic bond orders (1,2,3) are mapped; aromaticity/stereo are not inferred here.
        """
        from rdkit import Chem
        from rdkit.Chem import rdchem
        from rdkit.Geometry import rdGeometry as Geometry

        rw = Chem.RWMol()
        atom_index: dict[Entity, int] = {}
        atoms = list(atomistic.atoms)

        # add atoms
        for ent in atoms:
            sym = _element_symbol(ent)
            rd_atom = Chem.Atom(sym if sym != "X" else 0)
            charge = ent.get("charge")
            if isinstance(charge, int):
                rd_atom.SetFormalCharge(charge)
            isotope = ent.get("isotope")
            if isinstance(isotope, int):
                rd_atom.SetIsotope(isotope)
            hcount = ent.get("hcount")
            if isinstance(hcount, int):
                rd_atom.SetNumExplicitHs(hcount)
            idx = rw.AddAtom(rd_atom)
            atom_index[ent] = idx

        # add bonds
        for bond in atomistic.bonds:
            eps = getattr(bond, "endpoints", None)
            if not isinstance(eps, list) or len(eps) != 2:
                continue
            a, b = eps[0], eps[1]
            if a not in atom_index or b not in atom_index:
                # skip links to non-atom entities
                continue
            btype = _bond_type_from(rdchem, bond.get("order"), bond.get("kind"))
            rw.AddBond(atom_index[a], atom_index[b], btype)

        mol = rw.GetMol()

        # coordinates if present
        def _has_pos(e: Entity) -> bool:
            pos = e.get("pos")
            return isinstance(pos, list) and len(pos) == 3

        has_any_pos = any(_has_pos(ent) for ent in atoms)
        if has_any_pos and len(atoms) == mol.GetNumAtoms():
            conf = Chem.Conformer(len(atoms))
            for i, ent in enumerate(atoms):
                pos = ent.get("pos")
                if isinstance(pos, list) and len(pos) == 3:
                    x, y, z = float(pos[0]), float(pos[1]), float(pos[2])
                    conf.SetAtomPosition(i, Geometry.Point3D(x, y, z))
            mol.AddConformer(conf, assignId=True)

        # sanitize
        Chem.SanitizeMol(mol)
        return mol

    @staticmethod
    def monomer_to_rdkit(monomer: Monomer):
        """Convert a Monomer (wrapper of Assembly) to an RDKit Chem.Mol."""
        return RDKit.atomistic_to_rdkit(monomer.unwrap())




# def _iter_atoms(asm: Assembly) -> Iterable[Entity]:
# 	# Prefer strongly-typed AA atoms if registered
# 	try:
# 		bucket = asm.entities.bucket(Atom)
# 		return bucket
# 	except KeyError:
# 		pass
# 	# Fallback: any entity carrying element/symbol
# 	atoms: list[Entity] = []
# 	for ecls in asm.entities.classes():
# 		for ent in asm.entities.bucket(ecls):
# 			if isinstance(ent, Entity) and ("element" in ent or "symbol" in ent):
# 				atoms.append(ent)
# 	return atoms


# def _iter_bonds(asm: Assembly) -> Iterable[Bond]:
# 	# Prefer AA Bond bucket
# 	try:
# 		return asm.links.bucket(Bond)  # type: ignore[return-value]
# 	except KeyError:
# 		pass
# 	# Fallback: any two-endpoint link
# 	bonds: list[Bond] = []
# 	for lcls in asm.links.classes():
# 		for link in asm.links.bucket(lcls):
# 			eps = getattr(link, "endpoints", [])
# 			if isinstance(eps, list) and len(eps) == 2:
# 				bonds.append(link)  # type: ignore[arg-type]
# 	return bonds


def _element_symbol(atom: Entity) -> str:
	sym = atom.get("element") or atom.get("symbol") or "X"
	s = str(sym)
	if len(s) == 1:
		return s.upper()
	return s[0].upper() + s[1:].lower()


def _bond_type_from(rdchem, order: int | None, kind: str | None):
	if kind == ":":
		return rdchem.BondType.SINGLE
	if order == 2 or kind == "=":
		return rdchem.BondType.DOUBLE
	if order == 3 or kind == "#":
		return rdchem.BondType.TRIPLE
	return rdchem.BondType.SINGLE


