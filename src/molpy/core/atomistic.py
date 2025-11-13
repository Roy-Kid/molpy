from .entity import ConnectivityMixin, Entity, Assembly, Link, MembershipMixin, SpatialMixin, Entities
from typing import Any, cast, Self


class Atom(Entity):
    """Atom entity (expects optional keys like {"type": "C", "pos": [...]})"""
    def __repr__(self) -> str:
        identifier: str
        if "symbol" in self.data:
            identifier = str(self.data["symbol"])
        elif "type" in self.data:
            identifier = str(self.data["type"])
        else:
            identifier = str(id(self))
        return f"<Atom: {identifier}>"


class Bond(Link):
    def __init__(self, a: Atom, b: Atom, /, **attrs: Any):
        super().__init__([a, b], **attrs)

    def __repr__(self) -> str:
        return f"<Bond: {self.itom} - {self.jtom}>"

    @property
    def itom(self) -> Atom:
        return self.endpoints[0]

    @property
    def jtom(self) -> Atom:
        return self.endpoints[1]

class Angle(Link):
    def __init__(self, a: Atom, b: Atom, c: Atom, /, **attrs: Any):
        super().__init__([a, b, c], **attrs)

    def __repr__(self) -> str:
        return f"<Angle: {self.itom} - {self.jtom} - {self.ktom}>"

    @property
    def itom(self) -> Atom:
        return self.endpoints[0]

    @property
    def jtom(self) -> Atom:
        return self.endpoints[1]

    @property
    def ktom(self) -> Atom:
        return self.endpoints[2]


class Dihedral(Link):
    """Dihedral (torsion) angle between four atoms"""
    def __init__(self, a: Atom, b: Atom, c: Atom, d: Atom, /, **attrs: Any):
        super().__init__([a, b, c, d], **attrs)

    def __repr__(self) -> str:
        return f"<Dihedral: {self.itom} - {self.jtom} - {self.ktom} - {self.ltom}>"

    @property
    def itom(self) -> Atom:
        return self.endpoints[0]

    @property
    def jtom(self) -> Atom:
        return self.endpoints[1]

    @property
    def ktom(self) -> Atom:
        return self.endpoints[2]

    @property
    def ltom(self) -> Atom:
        return self.endpoints[3]


class Atomistic(Assembly, MembershipMixin, SpatialMixin, ConnectivityMixin):

    def __init__(self, **props) -> None:
        super().__init__(**props)
        # Call __post_init__ if it exists (for template pattern)
        if hasattr(self, '__post_init__'):
            # Get the method from the actual class, not from parent
            for klass in type(self).__mro__:
                if klass is Atomistic:
                    break
                if '__post_init__' in klass.__dict__:
                    klass.__dict__['__post_init__'](self, **props)
                    break
        self.entities.register_type(Atom)
        self.links.register_type(Bond)
        self.links.register_type(Angle)
        self.links.register_type(Dihedral)

    def merge(self, other: Assembly) -> Self:
        """
        Transfer all entities and links from another assembly into this one.
        
        Args:
            other: Assembly to merge into self
        
        Returns:
            Self for method chaining
        
        Example:
            >>> asm1.merge(asm2)  # Transfers asm2 into asm1
            >>> # asm2 should not be used after this!
        """
        return super().merge(other)

    @property
    def atoms(self) -> Entities[Atom]:
        return self.entities.bucket(Atom)

    @property
    def bonds(self) -> Entities[Bond]:
        return self.links.bucket(Bond)
    
    @property
    def angles(self) -> Entities[Angle]:
        return self.links.bucket(Angle)
    
    @property
    def dihedrals(self) -> Entities[Dihedral]:
        return self.links.bucket(Dihedral)

    @property
    def symbols(self) -> list[str]:
        atoms = list(self.atoms)
        return [str(a.get("symbol", "")) for a in atoms]
    
    @property
    def xyz(self):
        """
        Get atomic positions as numpy array.
        
        Returns:
            Nx3 array of atomic coordinates, or list of lists if numpy not available
        """
        atoms = list(self.atoms)
        positions = []
        for atom in atoms:
            pos = atom.get('xyz', atom.get('pos', [0.0, 0.0, 0.0]))
            positions.append(pos)
        
        try:
            import numpy as np
            return np.array(positions)
        except ImportError:
            return positions
    
    @property
    def positions(self):
        """Alias for xyz property."""
        return self.xyz
    
    def __repr__(self) -> str:
        """
        Return a concise representation of the atomistic structure.
        
        Shows:
        - Number of atoms (with element composition)
        - Number of bonds
        - Bounding box if positions available
        """
        atoms = self.atoms
        bonds = self.bonds
        
        # Count atoms by symbol
        from collections import Counter
        symbols = [a.get("symbol", "?") for a in atoms]
        symbol_counts = Counter(symbols)
        
        # Format composition
        if len(symbol_counts) <= 5:
            composition = " ".join(f"{sym}:{cnt}" for sym, cnt in sorted(symbol_counts.items()))
        else:
            composition = f"{len(symbol_counts)} types"
        
        # Check if we have positions
        has_coords = any("pos" in a or "xyz" in a for a in atoms)
        
        parts = [f"Atomistic"]
        parts.append(f"{len(atoms)} atoms ({composition})")
        parts.append(f"{len(bonds)} bonds")
        
        if has_coords:
            parts.append("with coords")
        
        return f"<{', '.join(parts)}>"
    
    def add_atom(self, **attrs: Any) -> Atom:
        atom = Atom(**attrs)
        self.entities.add(atom)
        return atom
    
    def add_bond(self, a: Atom, b: Atom, **attrs: Any) -> Bond:
        bond = Bond(a, b, **attrs)
        self.links.add(bond)
        return bond

    # ========== Spatial Operations (return self for chaining) ==========
    
    def move(self, delta: list[float], *, entity_type: type[Entity]=Atom) -> 'Atomistic':
        """Move all entities by delta. Returns self for chaining."""
        super().move(delta, entity_type=entity_type)
        return self

    def rotate(self, axis: list[float], angle: float, about: list[float] | None = None, *, entity_type: type[Entity]=Atom) -> 'Atomistic':
        """Rotate entities around axis. Returns self for chaining."""
        super().rotate(axis, angle, about=about, entity_type=entity_type)
        return self

    def scale(self, factor: float, about: list[float] | None = None, *, entity_type: type[Entity]=Atom) -> 'Atomistic':
        """Scale entities by factor. Returns self for chaining."""
        super().scale(factor, about=about, entity_type=entity_type)
        return self

    def align(
        self,
        a: Entity,
        b: Entity,
        *,
        a_dir: list[float] | None = None,
        b_dir: list[float] | None = None,
        flip: bool = False,
        entity_type: type[Entity] = Atom,
    ) -> 'Atomistic':
        """Align entities. Returns self for chaining."""
        super().align(a, b, a_dir=a_dir, b_dir=b_dir, flip=flip, entity_type=entity_type)
        return self

    # ========== System Composition ==========
    
    def __iadd__(self, other: 'Atomistic') -> 'Atomistic':
        """
        Merge another Atomistic into this one (in-place).
        
        Example:
            system = Water()
            system += Water().move([5, 0, 0])
            system += Methane().move([10, 0, 0])
        """
        self.merge(other)
        return self
    
    def __add__(self, other: 'Atomistic') -> 'Atomistic':
        """
        Create a new Atomistic by merging two systems.
        
        Example:
            combined = water1 + water2
        """
        result = self.copy()
        result.merge(other)
        return result
    
    def replicate(self, n: int, transform=None) -> 'Atomistic':
        """
        Create n copies and merge them into a new system.
        
        Args:
            n: Number of copies to create
            transform: Optional callable(copy, index) -> None to transform each copy
        
        Example:
            # Create 10 waters in a line
            waters = Water().replicate(10, lambda mol, i: mol.move([i*5, 0, 0]))
            
            # Create 3x3 grid
            grid = Methane().replicate(9, lambda mol, i: mol.move([i%3*5, i//3*5, 0]))
        """
        result = type(self)()  # Empty system of same type
        
        for i in range(n):
            replica = self.copy()
            if transform is not None:
                transform(replica, i)
            result.merge(replica)
        
        return result

    def __len__(self) -> int:
        return len(self.atoms)

    def get_topo(self, gen_angle: bool=False, gen_dihe=False):

        vertrices = {}
        for i, atom in enumerate(self.entities[Atom]):
            vertrices[atom] = i
        edges = []
        for bond in self.links[Bond]:
            edges.append((vertrices[bond.itom], vertrices[bond.jtom]))
        atoms = list(vertrices.keys())
        from .topology import Topology
        topo = Topology(len(vertrices), edges=edges)
        if gen_angle:
            for angle in topo.angles:
                angle = angle.tolist()
                self.links.add(Angle(atoms[angle[0]], atoms[angle[1]], atoms[angle[2]]))
        if gen_dihe:
            for dihe in topo.dihedrals:
                dihe = dihe.tolist()
                self.links.add(Dihedral(atoms[dihe[0]], atoms[dihe[1]], atoms[dihe[2]], atoms[dihe[3]]))

        return topo
