from __future__ import annotations

"""
Flavor-aware IR definitions for SMILES, BigSMILES, and G-BigSMILES.

Includes:
- Core types (Span, IRDiagnostic, IRBond, IRRingBond, IRBranch, IRAtom)
- SMILES IR (IRSmilesDot, IRSmilesMolecule, IRSmilesDocument)
- BigSMILES IR (IRBig* types) and compatibility aliases used by existing parser/converters
- G-BigSMILES IR (parameters/templates minimal superset)
- Serialization helpers: to_dict, to_json, from_dict, from_json
- Degradation utilities: degrade_to_smiles_big, degrade_to_smiles_gbig

Notes
- To minimize churn, we expose compatibility aliases mapping historical names
  (e.g. IRMolecule, IRDocument) to BigSMILES counterparts. This keeps existing
  parser/convert modules working while we introduce explicit flavored documents.
- Dict/JSON helpers are flavor-aware and include a top-level "flavor" key.
"""

from dataclasses import asdict, dataclass, field, is_dataclass
from typing import Any, Literal, Mapping, Callable
import json


# ----------------------- Core (shared) -----------------------


@dataclass(slots=True)
class Span:
    start: int
    end: int
    line: int
    column: int


@dataclass(slots=True)
class IRDiagnostic:
    level: Literal["error", "warning", "info"]
    message: str
    span: Span | None


@dataclass(slots=True)
class IRBond:
    kind: Literal["-", "=", "#", ";", ":", "/", "\\"]


@dataclass(slots=True)
class IRRingBond:
    bond: IRBond | None
    index: int
    span: Span | None


@dataclass(slots=True)
class IRBranch:
    bond: IRBond | None
    content: list[Any] = field(default_factory=list)
    span: Span | None = None


@dataclass(slots=True)
class IRAtom:
    kind: Literal["bracket", "aliphatic", "aromatic", "star"]
    symbol: str | None
    isotope: int | None
    chiral: str | None
    hcount: int | None
    charge: int | None
    clazz: int | None
    rings: list[IRRingBond] = field(default_factory=list)
    branches: list[IRBranch] = field(default_factory=list)
    span: Span | None = None
    bond_to_prev: IRBond | None = None


# Forward ref for type checkers
IRCorePart = IRAtom  # minimal, extended per flavor below


# ----------------------- SMILES -----------------------


@dataclass(slots=True)
class IRSmilesDot:
    """Plain SMILES dot separator (no system fields)."""

    span: Span | None


IRSmilesPart = IRAtom | IRSmilesDot


@dataclass(slots=True)
class IRSmilesMolecule:
    parts: list[IRSmilesPart] = field(default_factory=list)
    span: Span | None = None


@dataclass(slots=True)
class IRSmilesDocument:
    flavor: Literal["smiles"] = "smiles"
    molecules: list[IRSmilesMolecule] = field(default_factory=list)
    features: set[str] = field(default_factory=set)
    errors: list[IRDiagnostic] = field(default_factory=list)


# ----------------------- BigSMILES -----------------------


@dataclass(slots=True)
class IRBigGeneration:
    values: list[float] = field(default_factory=list)
    span: Span | None = None


@dataclass(slots=True)
class IRBigIndexExpr:
    left: "IRBigIndexExpr | int"
    op: Literal["~", "&"] | None
    right: "IRBigIndexExpr | int | None"
    unary: Literal["!"] | None = None
    span: Span | None = None


@dataclass(slots=True)
class IRBigIndexContext:
    expr: IRBigIndexExpr
    kv: dict[str, str] = field(default_factory=dict)
    span: Span | None = None


@dataclass(slots=True)
class IRBigNonCovalent:
    label: int | None
    context: IRBigIndexContext | None
    span: Span | None


@dataclass(slots=True)
class IRBigConnector:
    mode: Literal["simple", "ladder", "non_covalent"]
    connector: Literal["$", "<", ">", ":"] | None
    index: int | None
    generation: IRBigGeneration | None = None
    inner: list["IRBigConnector"] | None = None
    noncovalent: IRBigNonCovalent | None = None
    span: Span | None = None


@dataclass(slots=True)
class IRBigDistribution:
    kind: Literal[
        "flory_schulz",
        "schulz_zimm",
        "gauss",
        "uniform",
        "log_normal",
        "poisson",
    ]
    p1: float | None
    p2: float | None
    span: Span | None


@dataclass(slots=True)
class IRBigSystemDot:
    """System-level dot in Big/G-BigSMILES."""

    system_size: float | None
    span: Span | None


@dataclass(slots=True)
class IRBigStochasticBlock:
    left_terminal: IRBigConnector
    right_terminal: IRBigConnector
    units: list["IRBigSmilesMolecule"] = field(default_factory=list)
    end_group: list["IRBigSmilesMolecule"] | None = None
    distribution: IRBigDistribution | None = None
    span: Span | None = None


@dataclass(slots=True)
class IRBigFragmentRef:
    name: str
    span: Span | None


@dataclass(slots=True)
class IRBigFragmentDef:
    name: str
    molecules: list["IRBigSmilesMolecule"] = field(default_factory=list)
    span: Span | None = None


IRBigPart = IRAtom | IRBigConnector | IRBigStochasticBlock | IRBigFragmentRef | IRBigSystemDot


@dataclass(slots=True)
class IRBigSmilesMolecule:
    parts: list[IRBigPart] = field(default_factory=list)
    generation: IRBigSystemDot | None = None
    span: Span | None = None


@dataclass(slots=True)
class IRBigSmilesDocument:
    flavor: Literal["bigsmiles"] = "bigsmiles"
    molecules: list[IRBigSmilesMolecule] = field(default_factory=list)
    fragments: list[IRBigFragmentDef] = field(default_factory=list)
    features: set[str] = field(default_factory=set)
    errors: list[IRDiagnostic] = field(default_factory=list)


# ----------------------- G-BigSMILES (superset) -----------------------


@dataclass(slots=True)
class IRGParameter:
    name: str
    value: str | int | float | None
    span: Span | None


@dataclass(slots=True)
class IRGConstraint:
    expr: str
    span: Span | None


@dataclass(slots=True)
class IRGTemplateUnit:
    unit: IRBigSmilesMolecule
    params: list[IRGParameter] = field(default_factory=list)
    constraints: list[IRGConstraint] = field(default_factory=list)
    span: Span | None = None


IRGPart = IRBigPart | IRGTemplateUnit


@dataclass(slots=True)
class IRGBigSmilesMolecule:
    parts: list[IRGPart] = field(default_factory=list)
    generation: IRBigSystemDot | None = None
    span: Span | None = None


@dataclass(slots=True)
class IRGBigSmilesDocument:
    flavor: Literal["gbigsmiles"] = "gbigsmiles"
    molecules: list[IRGBigSmilesMolecule] = field(default_factory=list)
    fragments: list[IRBigFragmentDef] = field(default_factory=list)
    features: set[str] = field(default_factory=set)
    errors: list[IRDiagnostic] = field(default_factory=list)


# ----------------------- Compatibility aliases -----------------------

# These keep existing code working while we migrate call sites.
# Old unified IR names map to BigSMILES equivalents.
IRGeneration = IRBigGeneration
IRIndexExpr = IRBigIndexExpr
IRIndexContext = IRBigIndexContext
IRNonCovalent = IRBigNonCovalent
IRBondDescriptor = IRBigConnector
IRStochasticDistribution = IRBigDistribution
IRDot = IRBigSystemDot
IRStochasticObject = IRBigStochasticBlock
IRFragmentRef = IRBigFragmentRef
IRFragmentDef = IRBigFragmentDef
IRPart = IRBigPart

# Compatibility IRMolecule/IRDocument (legacy unified IR used by current parser)
@dataclass(slots=True)
class IRMolecule:
    parts: list[IRPart] = field(default_factory=list)
    generation: "IRDotGeneration | None" = None
    span: Span | None = None


@dataclass(slots=True)
class IRDocument:
    molecules: list[IRMolecule] = field(default_factory=list)
    fragments: list[IRFragmentDef] = field(default_factory=list)
    features: set[str] = field(default_factory=set)
    errors: list[IRDiagnostic] = field(default_factory=list)

# Legacy helper used by existing parser to wrap system dot
@dataclass(slots=True)
class IRDotGeneration:
    dot: IRDot
    span: Span | None = None


# ----------------------- Serialization helpers -----------------------


def _strip_nones(obj: Any) -> Any:
    """Recursively drop None fields from dataclass dicts and lists."""
    if is_dataclass(obj) and not isinstance(obj, type):
        return _strip_nones(asdict(obj))
    if isinstance(obj, dict):
        out: dict[str, Any] = {}
        for k, v in obj.items():
            if v is None:
                continue
            sv = _strip_nones(v)
            out[k] = sv
        return out
    if isinstance(obj, list):
        return [_strip_nones(x) for x in obj]
    if isinstance(obj, set):
        return sorted(obj)
    return obj


def to_dict(document: IRSmilesDocument | IRBigSmilesDocument | IRGBigSmilesDocument) -> dict[str, Any]:
    data = _strip_nones(document)
    # Ensure flavor is present and features serialized as list
    if "flavor" not in data:
        if isinstance(document, IRSmilesDocument):
            data["flavor"] = "smiles"
        elif isinstance(document, IRGBigSmilesDocument):
            data["flavor"] = "gbigsmiles"
        else:
            data["flavor"] = "bigsmiles"
    return data  # already stripped


def to_json(document: IRSmilesDocument | IRBigSmilesDocument | IRGBigSmilesDocument, *, indent: int | None = None) -> str:
    return json.dumps(to_dict(document), indent=indent)


def from_dict(data: Mapping[str, Any]) -> IRSmilesDocument | IRBigSmilesDocument | IRGBigSmilesDocument:
    flavor = str(data.get("flavor") or "bigsmiles")
    if flavor == "smiles":
        return _smiles_doc_from_dict(data)
    if flavor == "gbigsmiles":
        return _gbig_doc_from_dict(data)
    return _big_doc_from_dict(data)


def from_json(payload: str) -> IRSmilesDocument | IRBigSmilesDocument | IRGBigSmilesDocument:
    return from_dict(json.loads(payload))


# ---------- builders: shared small helpers ----------


def _span_from_dict(data: Mapping[str, Any] | None) -> Span | None:
    if not data:
        return None
    return Span(
        start=int(data["start"]),
        end=int(data["end"]),
        line=int(data["line"]),
        column=int(data["column"]),
    )


def _bond_from_dict(data: Mapping[str, Any] | None) -> IRBond | None:
    if not data:
        return None
    return IRBond(kind=data["kind"])  # type: ignore[arg-type]


def _ring_bond_from_dict(data: Mapping[str, Any]) -> IRRingBond:
    return IRRingBond(
        bond=_bond_from_dict(data.get("bond")),
        index=int(data["index"]),
        span=_span_from_dict(data.get("span")),
    )


def _branch_from_dict_core(
    data: Mapping[str, Any],
    part_builder: Callable[[Mapping[str, Any]], Any],
) -> IRBranch:
    return IRBranch(
        bond=_bond_from_dict(data.get("bond")),
        content=[part_builder(item) for item in data.get("content", [])],
        span=_span_from_dict(data.get("span")),
    )


def _atom_from_dict_core(data: Mapping[str, Any]) -> IRAtom:
    return IRAtom(
        kind=data["kind"],
        symbol=data.get("symbol"),
        isotope=data.get("isotope"),
        chiral=data.get("chiral"),
        hcount=data.get("hcount"),
        charge=data.get("charge"),
        clazz=data.get("clazz"),
        rings=[_ring_bond_from_dict(item) for item in data.get("rings", [])],
        branches=[],  # branches handled by flavor-specific builder to allow part types
        span=_span_from_dict(data.get("span")),
        bond_to_prev=_bond_from_dict(data.get("bond_to_prev")),
    )


# ---------- SMILES builders ----------


def _smiles_branch_from_dict(data: Mapping[str, Any]) -> IRBranch:
    return _branch_from_dict_core(data, _smiles_part_from_dict)


def _smiles_part_from_dict(data: Mapping[str, Any]) -> IRSmilesPart:
    if set(data.keys()) == {"span"}:
        return IRSmilesDot(span=_span_from_dict(data.get("span")))
    atom = _atom_from_dict_core(data)
    # fill branches if present
    atom.branches = [_smiles_branch_from_dict(item) for item in data.get("branches", [])]
    return atom


def _smiles_molecule_from_dict(data: Mapping[str, Any]) -> IRSmilesMolecule:
    return IRSmilesMolecule(
        parts=[_smiles_part_from_dict(item) for item in data.get("parts", [])],
        span=_span_from_dict(data.get("span")),
    )


def _smiles_doc_from_dict(data: Mapping[str, Any]) -> IRSmilesDocument:
    return IRSmilesDocument(
        flavor="smiles",
        molecules=[_smiles_molecule_from_dict(item) for item in data.get("molecules", [])],
        features=set(data.get("features", [])),
        errors=[_diagnostic_from_dict(item) for item in data.get("errors", [])],
    )


# ---------- BigSMILES builders ----------


def _big_generation_from_dict(data: Mapping[str, Any] | None) -> IRBigGeneration | None:
    if not data:
        return None
    return IRBigGeneration(values=[float(v) for v in data.get("values", [])], span=_span_from_dict(data.get("span")))


def _big_index_expr_value(value: Any) -> IRBigIndexExpr | int:
    if isinstance(value, Mapping):
        return _big_index_expr_from_dict(value)
    return int(value)


def _big_index_expr_from_dict(data: Mapping[str, Any]) -> IRBigIndexExpr:
    right_raw = data.get("right")
    right: IRBigIndexExpr | int | None
    if right_raw is None:
        right = None
    elif isinstance(right_raw, Mapping):
        right = _big_index_expr_from_dict(right_raw)
    else:
        right = int(right_raw)

    left_raw = data.get("left")
    left = _big_index_expr_value(left_raw) if left_raw is not None else 0

    return IRBigIndexExpr(
        left=left,
        op=data.get("op"),
        right=right,
        unary=data.get("unary"),
        span=_span_from_dict(data.get("span")),
    )


def _big_index_context_from_dict(data: Mapping[str, Any] | None) -> IRBigIndexContext | None:
    if not data:
        return None
    return IRBigIndexContext(
        expr=_big_index_expr_from_dict(data["expr"]),
        kv={str(k): str(v) for k, v in data.get("kv", {}).items()},
        span=_span_from_dict(data.get("span")),
    )


def _big_non_covalent_from_dict(data: Mapping[str, Any] | None) -> IRBigNonCovalent | None:
    if not data:
        return None
    return IRBigNonCovalent(
        label=data.get("label"),
        context=_big_index_context_from_dict(data.get("context")),
        span=_span_from_dict(data.get("span")),
    )


def _big_connector_from_dict(data: Mapping[str, Any]) -> IRBigConnector:
    inner_data = data.get("inner")
    inner = None
    if inner_data is not None:
        inner = [_big_connector_from_dict(item) for item in inner_data]
    return IRBigConnector(
        mode=data["mode"],
        connector=data.get("connector"),
        index=data.get("index"),
        generation=_big_generation_from_dict(data.get("generation")),
        inner=inner,
        noncovalent=_big_non_covalent_from_dict(data.get("noncovalent")),
        span=_span_from_dict(data.get("span")),
    )


def _big_dot_from_dict(data: Mapping[str, Any]) -> IRBigSystemDot:
    return IRBigSystemDot(
        system_size=(float(data["system_size"]) if data.get("system_size") is not None else None),
        span=_span_from_dict(data.get("span")),
    )


def _big_part_from_dict(data: Mapping[str, Any]) -> IRBigPart:
    if "left_terminal" in data:
        return _big_stochastic_block_from_dict(data)
    if "mode" in data and data.get("mode") in {"simple", "ladder", "non_covalent"}:
        return _big_connector_from_dict(data)
    if "system_size" in data:
        return _big_dot_from_dict(data)
    if set(data.keys()) == {"name", "span"}:
        return _big_fragment_ref_from_dict(data)
    return _atom_from_dict_core(data)


def _big_molecule_from_dict(data: Mapping[str, Any]) -> IRBigSmilesMolecule:
    # fill branches for atoms after parts are created (requires part builder)
    parts = [_big_part_from_dict(item) for item in data.get("parts", [])]
    # rebuild branches on atoms using big-part context
    for part in parts:
        if isinstance(part, IRAtom):
            branches_data = next((d.get("branches", []) for d in data.get("parts", []) if isinstance(d, Mapping) and d.get("span") == (part.span and asdict(part.span))), [])  # best-effort
            part.branches = [_branch_from_dict_core(bd, _big_part_from_dict) for bd in branches_data]
    return IRBigSmilesMolecule(
        parts=parts,
        generation=_big_dot_from_dict(data["generation"]) if data.get("generation") else None,
        span=_span_from_dict(data.get("span")),
    )


def _big_fragment_ref_from_dict(data: Mapping[str, Any]) -> IRBigFragmentRef:
    return IRBigFragmentRef(name=data["name"], span=_span_from_dict(data.get("span")))


def _big_fragment_def_from_dict(data: Mapping[str, Any]) -> IRBigFragmentDef:
    return IRBigFragmentDef(
        name=data["name"],
        molecules=[_big_molecule_from_dict(item) for item in data.get("molecules", [])],
        span=_span_from_dict(data.get("span")),
    )


def _big_distribution_from_dict(data: Mapping[str, Any] | None) -> IRBigDistribution | None:
    if not data:
        return None
    return IRBigDistribution(
        kind=data["kind"],
        p1=(float(data["p1"]) if data.get("p1") is not None else None),
        p2=(float(data["p2"]) if data.get("p2") is not None else None),
        span=_span_from_dict(data.get("span")),
    )


def _big_stochastic_block_from_dict(data: Mapping[str, Any]) -> IRBigStochasticBlock:
    end_group_data = data.get("end_group")
    end_group = None
    if end_group_data is not None:
        end_group = [_big_molecule_from_dict(item) for item in end_group_data]
    return IRBigStochasticBlock(
        left_terminal=_big_connector_from_dict(data["left_terminal"]),
        units=[_big_molecule_from_dict(item) for item in data.get("units", [])],
        end_group=end_group,
        right_terminal=_big_connector_from_dict(data["right_terminal"]),
        distribution=_big_distribution_from_dict(data.get("distribution")),
        span=_span_from_dict(data.get("span")),
    )


def _big_doc_from_dict(data: Mapping[str, Any]) -> IRBigSmilesDocument:
    return IRBigSmilesDocument(
        flavor="bigsmiles",
        molecules=[_big_molecule_from_dict(item) for item in data.get("molecules", [])],
        fragments=[_big_fragment_def_from_dict(item) for item in data.get("fragments", [])],
        features=set(data.get("features", [])),
        errors=[_diagnostic_from_dict(item) for item in data.get("errors", [])],
    )


# ---------- G-Big builders ----------


def _g_param_from_dict(data: Mapping[str, Any]) -> IRGParameter:
    return IRGParameter(name=data["name"], value=data.get("value"), span=_span_from_dict(data.get("span")))


def _g_constraint_from_dict(data: Mapping[str, Any]) -> IRGConstraint:
    return IRGConstraint(expr=data["expr"], span=_span_from_dict(data.get("span")))


def _g_template_unit_from_dict(data: Mapping[str, Any]) -> IRGTemplateUnit:
    return IRGTemplateUnit(
        unit=_big_molecule_from_dict(data["unit"]),
        params=[_g_param_from_dict(p) for p in data.get("params", [])],
        constraints=[_g_constraint_from_dict(c) for c in data.get("constraints", [])],
        span=_span_from_dict(data.get("span")),
    )


def _g_part_from_dict(data: Mapping[str, Any]) -> IRGPart:
    if "unit" in data and "params" in data:
        return _g_template_unit_from_dict(data)
    return _big_part_from_dict(data)


def _gbig_molecule_from_dict(data: Mapping[str, Any]) -> IRGBigSmilesMolecule:
    parts = [_g_part_from_dict(item) for item in data.get("parts", [])]
    return IRGBigSmilesMolecule(
        parts=parts,
        generation=_big_dot_from_dict(data["generation"]) if data.get("generation") else None,
        span=_span_from_dict(data.get("span")),
    )


def _gbig_doc_from_dict(data: Mapping[str, Any]) -> IRGBigSmilesDocument:
    return IRGBigSmilesDocument(
        flavor="gbigsmiles",
        molecules=[_gbig_molecule_from_dict(item) for item in data.get("molecules", [])],
        fragments=[_big_fragment_def_from_dict(item) for item in data.get("fragments", [])],
        features=set(data.get("features", [])),
        errors=[_diagnostic_from_dict(item) for item in data.get("errors", [])],
    )


def _diagnostic_from_dict(data: Mapping[str, Any]) -> IRDiagnostic:
    return IRDiagnostic(level=data["level"], message=data["message"], span=_span_from_dict(data.get("span")))


# ----------------------- Degradation utilities -----------------------


def _warn(msg: str, span: Span | None) -> IRDiagnostic:
    return IRDiagnostic(level="warning", message=msg, span=span)


def degrade_to_smiles_big(doc: IRBigSmilesDocument, strict: bool) -> tuple[IRSmilesDocument | None, list[IRDiagnostic]]:
    """Attempt to degrade BigSMILES to plain SMILES.

    - Drops stochastic blocks and descriptors in strict mode (error) or non-strict (warning).
    - Keeps atoms and covalent bonds; maps descriptor connector ':' to a bond ':' warning.
    - System dots are dropped with a warning.
    """
    diags: list[IRDiagnostic] = []
    smiles_mols: list[IRSmilesMolecule] = []

    for mol in doc.molecules:
        parts: list[IRSmilesPart] = []
        for part in mol.parts:
            if isinstance(part, IRAtom):
                parts.append(part)
            elif isinstance(part, IRBigConnector):
                if strict:
                    diags.append(_warn("Dropped connector while degrading to SMILES", part.span))
                else:
                    # Represent as a ':'-bonded pseudo dot (best effort)
                    diags.append(_warn("Mapping connector to ':' bond is not representable in SMILES; dropping", part.span))
            elif isinstance(part, IRBigStochasticBlock):
                msg = "Stochastic block not representable in SMILES"
                if strict:
                    diags.append(_warn(msg, part.span))
                else:
                    diags.append(_warn(msg + "; contents dropped", part.span))
            elif isinstance(part, IRBigSystemDot):
                diags.append(_warn("System dot dropped while degrading to SMILES", part.span))
            else:
                # fragment refs are not supported inline here
                diags.append(_warn("Fragment/reference dropped while degrading to SMILES", getattr(part, "span", None)))
        smiles_mols.append(IRSmilesMolecule(parts=parts, span=mol.span))

    result = IRSmilesDocument(flavor="smiles", molecules=smiles_mols, features={"smiles"}, errors=diags)
    return result, diags


def degrade_to_smiles_gbig(doc: IRGBigSmilesDocument, strict: bool) -> tuple[IRSmilesDocument | None, list[IRDiagnostic]]:
    """Degrade G-BigSMILES via BigSMILES path (drops G-specific constructs)."""
    diags: list[IRDiagnostic] = []
    # First coerce to BigSMILES by dropping template units and parameters
    big = IRBigSmilesDocument(
        flavor="bigsmiles",
        molecules=[],
        fragments=doc.fragments,
        features=doc.features | {"bigsmiles"},
        errors=list(doc.errors),
    )
    for mol in doc.molecules:
        flat_parts: list[IRBigPart] = []
        for p in mol.parts:
            if isinstance(p, IRGTemplateUnit):
                # Inline unit parts when degrading
                flat_parts.extend(p.unit.parts)
            else:
                flat_parts.append(p)
        big.molecules.append(
            IRBigSmilesMolecule(parts=flat_parts, generation=mol.generation, span=mol.span)
        )
    diags.append(_warn("Dropped G-BigSMILES template parameters/constraints", None))
    return degrade_to_smiles_big(big, strict)


__all__ = [
    # Core
    "Span",
    "IRDiagnostic",
    "IRBond",
    "IRRingBond",
    "IRBranch",
    "IRAtom",
    # SMILES
    "IRSmilesDot",
    "IRSmilesPart",
    "IRSmilesMolecule",
    "IRSmilesDocument",
    # BigSMILES
    "IRBigGeneration",
    "IRBigIndexExpr",
    "IRBigIndexContext",
    "IRBigNonCovalent",
    "IRBigConnector",
    "IRBigDistribution",
    "IRBigSystemDot",
    "IRBigStochasticBlock",
    "IRBigFragmentRef",
    "IRBigFragmentDef",
    "IRBigPart",
    "IRBigSmilesMolecule",
    "IRBigSmilesDocument",
    # G-BigSMILES
    "IRGParameter",
    "IRGConstraint",
    "IRGTemplateUnit",
    "IRGPart",
    "IRGBigSmilesMolecule",
    "IRGBigSmilesDocument",
    # Compatibility aliases
    "IRGeneration",
    "IRIndexExpr",
    "IRIndexContext",
    "IRNonCovalent",
    "IRBondDescriptor",
    "IRStochasticDistribution",
    "IRDot",
    "IRStochasticObject",
    "IRFragmentRef",
    "IRFragmentDef",
    "IRPart",
    # Legacy unified IR
    "IRDotGeneration",
    "IRMolecule",
    "IRDocument",
    # Serialization & degradation
    "to_dict",
    "to_json",
    "from_dict",
    "from_json",
    "degrade_to_smiles_big",
    "degrade_to_smiles_gbig",
]
