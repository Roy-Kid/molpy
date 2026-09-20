"""CL&P ionic-liquid force-field typifier.

CL&P stays in molpy (OPLS-AA moved to the native core). It is an OPLS-AA *overlay*: the
built-in force field is ``oplsaa.xml`` with ``clp.xml`` layered on top (layer 1),
so CL&P atom types (imidazolium ring, alkyl chain, and the BF4/PF6/NTf2/FSI/dca
anions) override the OPLS base while OPLS remains the fallback.

Atom typing itself is SMARTS-based and that matcher is owned by the native core, so the
only thing this typifier's ``match`` does is ask the native core for the atom types and
hand them to :class:`~molpy.typifier.forcefield.ForceFieldParams` — exactly the
"SMARTS owned by the native core, CL&P parameters stay a molpy overlay" split.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, override

from molrs.ff.typifier import OPLSAATypifier

from molpy.core import fields
from molpy.core.atomistic import Atomistic
from molpy.typifier.base import Match, Typifier
from molpy.typifier.forcefield import ForceFieldParams

if TYPE_CHECKING:
    from collections.abc import Mapping

    from molpy.core.forcefield import ForceField
    from molpy.typifier.base import Annotation


@lru_cache(maxsize=1)
def _clp_molrs_typifier() -> OPLSAATypifier:
    """Shared, stateless the native core SMARTS typifier over the CL&P overlay (compiling
    the SMARTS engine once instead of per :class:`ClpTypifier` construction)."""
    # strict=False so molrs's own bonded matching never errors — only the
    # atom-level type/class it assigns via the CL&P SMARTS defs is harvested.
    from molpy.data import get_forcefield_path

    return OPLSAATypifier(
        Path(get_forcefield_path("clp.xml")).read_text(encoding="utf-8"), strict=False
    )


@lru_cache(maxsize=1)
def _load_clp_forcefield() -> ForceField:
    """Parse ``oplsaa.xml`` + ``clp.xml`` once natively (OPLS unit path + merge)."""
    from molpy.data import get_forcefield_path
    from molpy.io.forcefield.xml import read_xml_forcefield

    ff = read_xml_forcefield(get_forcefield_path("oplsaa.xml"))
    return read_xml_forcefield(get_forcefield_path("clp.xml"), ff, layer=1)


@lru_cache(maxsize=2)
def _default_clp_params(strict: bool) -> ForceFieldParams:
    """``ForceFieldParams`` over the built-in overlay is expensive (TypeClassIndex
    walks every OPLS atom type). Cache one instance per ``strict`` flag."""
    return ForceFieldParams(_load_clp_forcefield(), strict=strict)


class ClpTypifier(Typifier[Atomistic]):
    """CL&P ionic-liquid typifier — the native core SMARTS atom typing + molpy parameters.

    Args:
        forcefield: The CL&P-over-OPLS overlay; the built-in one by default.
        strict: Raise on an atom no SMARTS pattern matches, or a bonded term the
            force field does not parameterise.
    """

    def __init__(
        self, forcefield: ForceField | None = None, *, strict: bool = True
    ) -> None:
        self._strict = strict
        self._smarts = _clp_molrs_typifier()
        if forcefield is None:
            self.ff = _load_clp_forcefield()
            self._params = _default_clp_params(strict)
        else:
            self.ff = forcefield
            self._params = ForceFieldParams(forcefield, strict=strict)

    @override
    def match(self, graph: Atomistic) -> Match:
        return self._params.match(graph, self._atom_types(graph))

    def _atom_types(self, graph: Atomistic) -> list[Mapping[str, Annotation]]:
        """Ask the native core which CL&P type (and class) each atom carries."""
        typed = self._smarts.typify(graph).to_frame()["atoms"]
        names = typed[fields.TYPE]
        classes = typed["class"] if "class" in typed else [None] * len(names)

        out: list[Mapping[str, Annotation]] = []
        for atom, name, class_name in zip(graph.atoms, names, classes, strict=True):
            if name in (None, ""):
                if self._strict:
                    raise ValueError(f"CL&P: no atom type matched for atom {atom}")
                out.append({})
                continue
            annotation: dict[str, Annotation] = {fields.TYPE: str(name)}
            if class_name not in (None, ""):
                annotation["class"] = str(class_name)
            out.append(annotation)
        return out

    @staticmethod
    def load_forcefield() -> ForceField:
        """Load the built-in CL&P force field as an OPLS-AA overlay."""
        return _load_clp_forcefield()


__all__ = ["ClpTypifier"]
