from __future__ import annotations

import math
import random
import string
from dataclasses import dataclass, field
from typing import Optional, Sequence

from molpy.core.aa import Atom, Bond
from molpy.core.aa.base import Atomistic
from molpy.core.entity import Assembly
from molpy.core.entity import Entity
from molpy.core.wrappers.monomer import Monomer

from .smiles_ir import (
    IRAtom,
    IRBond,
    IRSmilesDot,
    IRSmilesDocument,
    IRSmilesMolecule,
    IRBigSmilesDocument,
    IRBigSmilesMolecule,
    IRGBigSmilesDocument,
    IRBigSystemDot,
    IRBigConnector,
    IRBigIndexContext,
    IRBigIndexExpr,
    IRBigStochasticBlock,
    IRBigDistribution,
)

__all__ = [
    "DistributionSpec",
    "PolymerChain",
    "StochasticSequence",
    "GBigSmilesMolecule",
    "GBigSmilesResult",
    "ir_to_atomistic",
    "ir_to_bigsmiles",
    "ir_to_gbigsmiles",
]


@dataclass(slots=True)
class DistributionSpec:
    """
    Lightweight description of a stochastic distribution declared in G-BigSMILES.
    """

    kind: str
    parameters: tuple[float | None, float | None]


@dataclass(slots=True)
class PolymerChain:
    """
    Concrete polymer chain instantiated from a stochastic specification.
    """

    labels: list[str]
    degree: int


@dataclass(slots=True)
class StochasticSequence:
    """
    Realization of a stochastic object, including sampled chains and metadata.
    """

    unit_labels: list[str]
    unit_weights: list[float]
    chains: list[PolymerChain] = field(default_factory=list)
    distribution: DistributionSpec | None = None
    left_terminal: dict[str, object] | None = None
    right_terminal: dict[str, object] | None = None
    end_group_labels: list[str] = field(default_factory=list)
    target_total: int | None = None


@dataclass(slots=True)
class GBigSmilesMolecule:
    """
    Container for per-molecule G-BigSMILES data.
    """

    sequences: list[StochasticSequence] = field(default_factory=list)
    system_size: int | None = None
    bond_descriptors: list[dict[str, object] | None] = field(default_factory=list)


@dataclass(slots=True)
class GBigSmilesResult:
    """
    Aggregate output for G-BigSMILES conversion.
    """

    monomers: dict[str, Monomer]
    molecules: list[GBigSmilesMolecule] = field(default_factory=list)


class _BigSmilesRegistry:
    """
    Tracks monomers generated during conversion and allocates unique labels.
    """

    def __init__(self) -> None:
        self.monomers: dict[str, Monomer] = {}
        self._label_index = 0

    def next_label(self) -> str:
        alphabet = string.ascii_uppercase
        idx = self._label_index
        if idx < len(alphabet):
            label = alphabet[idx]
        else:
            letter = alphabet[idx % len(alphabet)]
            cycle = idx // len(alphabet)
            label = f"{letter}{cycle}"
        self._label_index += 1
        return label

    def register_molecule(self, molecule: IRBigSmilesMolecule) -> str:
        label = self.next_label()
        self.monomers[label] = _molecule_to_monomer(molecule)
        return label

    def register_named(self, name: str, monomer: Monomer) -> None:
        self.monomers[name] = monomer


@dataclass(slots=True)
class _StochasticComponent:
    unit_labels: list[str]
    unit_weights: list[float]
    end_group_labels: list[str]
    distribution: IRBigDistribution | None
    left_descriptor: dict[str, object] | None
    right_descriptor: dict[str, object] | None


class _DistributionSampler:
    """
    Samples integer degrees of polymerization for supported distributions.
    """

    def __init__(self, distribution: IRBigDistribution, rng: random.Random) -> None:
        self._dist = distribution
        self._rng = rng

    def sample(self) -> int:
        kind = self._dist.kind
        if kind == "schulz_zimm":
            return self._sample_schulz_zimm()
        if kind == "flory_schulz":
            return self._sample_flory_schulz()
        if kind == "gauss":
            return self._sample_gauss()
        if kind == "uniform":
            return self._sample_uniform()
        if kind == "log_normal":
            return self._sample_log_normal()
        if kind == "poisson":
            return self._sample_poisson()
        raise ValueError(f"Unsupported stochastic distribution: {kind!r}")

    def _sample_schulz_zimm(self) -> int:
        Mw = float(self._dist.p1) if self._dist.p1 is not None else 0.0
        Mn = float(self._dist.p2) if self._dist.p2 is not None else 0.0
        if Mw <= 0 or Mn <= 0 or Mw <= Mn:
            # Fall back to mean-of-two to avoid runtime errors.
            return max(1, int(round(Mn)) if Mn > 0 else 1)
        pdi = Mw / Mn
        z = 1.0 / (pdi - 1.0)
        shape = z + 1.0
        scale = Mn / (z + 1.0)
        value = max(1, int(round(self._rng.gammavariate(shape, scale))))
        return value

    def _sample_flory_schulz(self) -> int:
        p = float(self._dist.p1) if self._dist.p1 is not None else 0.5
        p = min(max(p, 1e-9), 1 - 1e-9)
        length = 1
        # Geometric distribution with success probability (1 - p).
        while self._rng.random() < p:
            length += 1
        return length

    def _sample_gauss(self) -> int:
        mean = float(self._dist.p1) if self._dist.p1 is not None else 1.0
        stdev = float(self._dist.p2) if self._dist.p2 is not None else max(0.1 * mean, 1.0)
        sample = int(round(self._rng.gauss(mean, stdev)))
        return max(1, sample)

    def _sample_uniform(self) -> int:
        low = float(self._dist.p1) if self._dist.p1 is not None else 1.0
        high = float(self._dist.p2) if self._dist.p2 is not None else low
        if high < low:
            low, high = high, low
        sample = int(round(self._rng.uniform(low, high)))
        return max(1, sample)

    def _sample_log_normal(self) -> int:
        mu = float(self._dist.p1) if self._dist.p1 is not None else 1.0
        sigma = float(self._dist.p2) if self._dist.p2 is not None else 0.5
        mu = max(mu, 1e-6)
        sigma = max(sigma, 1e-6)
        value = self._rng.lognormvariate(math.log(mu), sigma)
        return max(1, int(round(value)))

    def _sample_poisson(self) -> int:
        lam = float(self._dist.p1) if self._dist.p1 is not None else 1.0
        lam = max(lam, 1e-6)
        if lam < 30:
            # Knuth algorithm
            L = math.exp(-lam)
            k = 0
            p = 1.0
            while p > L:
                k += 1
                p *= self._rng.random()
            return max(1, k - 1)
        # Normal approximation for large lambda.
        sample = int(round(self._rng.gauss(lam, math.sqrt(lam))))
        return max(1, sample)


def _build_stochastic_component(
    obj: IRBigStochasticBlock,
    registry: _BigSmilesRegistry,
) -> _StochasticComponent:
    unit_labels: list[str] = []
    unit_weights: list[float] = []
    for unit in obj.units:
        unit_labels.append(registry.register_molecule(unit))
        unit_weights.append(_extract_unit_weight(unit))

    end_group_labels: list[str] = []
    if obj.end_group:
        for end_group in obj.end_group:
            end_group_labels.append(registry.register_molecule(end_group))

    return _StochasticComponent(
        unit_labels=unit_labels,
        unit_weights=unit_weights,
        end_group_labels=end_group_labels,
        distribution=obj.distribution,
        left_descriptor=_descriptor_summary(obj.left_terminal),
        right_descriptor=_descriptor_summary(obj.right_terminal),
    )


def _extract_unit_weight(unit: IRBigSmilesMolecule) -> float:
    values: list[float] = []
    for part in unit.parts:
        if isinstance(part, IRBigConnector) and part.generation and part.generation.values:
            for number in part.generation.values:
                try:
                    values.append(float(number))
                except (TypeError, ValueError):
                    continue
    if not values:
        return 1.0
    return max(1e-6, sum(values) / len(values))


def _descriptor_summary(descriptor: IRBigConnector | None) -> dict[str, object] | None:
    if descriptor is None:
        return None
    summary: dict[str, object] = {
        "mode": descriptor.mode,
        "connector": descriptor.connector,
        "index": descriptor.index,
    }
    if descriptor.generation and descriptor.generation.values:
        summary["generation"] = [float(value) for value in descriptor.generation.values]
    if descriptor.inner:
        summary["inner"] = [_descriptor_summary(inner) for inner in descriptor.inner]
    if descriptor.noncovalent:
        summary["noncovalent"] = {
            "label": descriptor.noncovalent.label,
            "context": _index_context_summary(descriptor.noncovalent.context),
        }
    return summary


def _index_context_summary(context: IRBigIndexContext | None) -> dict[str, object] | None:
    if context is None:
        return None
    return {
        "expr": _index_expr_summary(context.expr),
        "kv": dict(context.kv),
    }


def _index_expr_summary(expr: IRBigIndexExpr | int | None) -> object:
    if expr is None or isinstance(expr, int):
        return expr
    return {
        "left": _index_expr_summary(expr.left),
        "op": expr.op,
        "right": _index_expr_summary(expr.right),
        "unary": expr.unary,
    }


def _normalize_weights(weights: Sequence[float]) -> list[float] | None:
    cleaned = [max(float(weight), 0.0) for weight in weights]
    total = sum(cleaned)
    if total <= 0:
        return None
    return cleaned


def _select_units(
    labels: Sequence[str],
    weights: Sequence[float] | None,
    length: int,
    rng: random.Random,
) -> list[str]:
    if length <= 0:
        return []
    if weights is not None:
        return rng.choices(list(labels), weights=weights, k=length)
    return [rng.choice(list(labels)) for _ in range(length)]


def _generate_polymer_chains(
    unit_labels: Sequence[str],
    weights: Sequence[float] | None,
    sampler: _DistributionSampler,
    rng: random.Random,
    *,
    target_total: int | None = None,
) -> list[PolymerChain]:
    if not unit_labels:
        return []
    weights_list = list(weights) if weights is not None else None
    chains: list[PolymerChain] = []
    total = 0
    iterations = 0

    while True:
        if target_total is not None and total >= target_total:
            break
        iterations += 1
        if iterations > 1_000_000:
            # Safety valve to avoid pathological infinite loops.
            break

        dp = sampler.sample()
        dp = max(1, dp)
        if target_total is not None:
            remaining = target_total - total
            if remaining <= 0:
                break
            if dp > remaining:
                dp = remaining
        labels = _select_units(unit_labels, weights_list, dp, rng)
        chains.append(PolymerChain(labels=labels, degree=dp))
        total += dp
        if target_total is None:
            break

    return chains


def _safe_int(value: float | None) -> int | None:
    if value is None:
        return None
    try:
        as_int = int(round(value))
    except (TypeError, ValueError):
        return None
    return as_int if as_int >= 0 else None


def _to_distribution_spec(distribution: IRBigDistribution | None) -> DistributionSpec | None:
    if distribution is None:
        return None
    p1 = float(distribution.p1) if distribution.p1 is not None else None
    p2 = float(distribution.p2) if distribution.p2 is not None else None
    return DistributionSpec(kind=distribution.kind, parameters=(p1, p2))


def _register_fragments(fragments: Sequence, registry: _BigSmilesRegistry) -> None:
    for fragment in fragments:
        if not fragment.molecules:
            continue
        fragment_builder = _AtomisticBuilder()
        for frag_mol in fragment.molecules:
            fragment_builder.build_molecule(frag_mol)
        registry.register_named(fragment.name, _atomistic_to_monomer(fragment_builder.atomistic))


def ir_to_atomistic(document: IRSmilesDocument | IRBigSmilesDocument | IRGBigSmilesDocument) -> Atomistic:
    """
    Convert a parsed IRDocument into an Atomistic structure.
    """
    builder = _AtomisticBuilder()
    for molecule in document.molecules:
        builder.build_molecule(molecule)
    return builder.atomistic


def ir_to_bigsmiles(
    document: IRBigSmilesDocument,
    *,
    sequence_length: Optional[int] = None,
    seed: int = 0,
) -> tuple[dict[str, Monomer], list[str]]:
    """
    Convert BigSMILES portions of the IR into monomer templates and a mock sequence.
    """
    registry = _BigSmilesRegistry()
    rng = random.Random(seed)

    components: list[_StochasticComponent] = []
    aggregated_labels: list[str] = []
    aggregated_weights: list[float] = []
    distribution_present = False

    for molecule in document.molecules:
        for part in molecule.parts:
            if isinstance(part, IRBigStochasticBlock):
                component = _build_stochastic_component(part, registry)
                components.append(component)
                aggregated_labels.extend(component.unit_labels)
                aggregated_weights.extend(component.unit_weights)
                if component.distribution is not None:
                    distribution_present = True

    # Fragment definitions become named monomers
    _register_fragments(document.fragments, registry)

    if not registry.monomers:
        return {}, []

    if distribution_present:
        sequence: list[str] = []
        for component in components:
            if component.distribution is None or not component.unit_labels:
                continue
            sampler = _DistributionSampler(component.distribution, rng)
            weights = _normalize_weights(component.unit_weights)
            chains = _generate_polymer_chains(
                component.unit_labels,
                weights,
                sampler,
                rng,
                target_total=None,
            )
            for chain in chains:
                sequence.extend(chain.labels)
        return registry.monomers, sequence

    # Otherwise, fallback to old behavior: optional fixed-length random sequence
    if sequence_length is None or sequence_length <= 0:
        return registry.monomers, []

    if aggregated_labels:
        weights = _normalize_weights(aggregated_weights)
        sequence = _select_units(aggregated_labels, weights, sequence_length, rng)
    else:
        pool = list(registry.monomers.keys())
        sequence = _select_units(pool, None, sequence_length, rng)
    return registry.monomers, sequence


def ir_to_gbigsmiles(document: IRGBigSmilesDocument, *, seed: int = 0) -> GBigSmilesResult:
    """
    Convert an IRDocument with G-BigSMILES features into structured stochastic data.
    """
    registry = _BigSmilesRegistry()
    rng = random.Random(seed)

    molecule_results: list[GBigSmilesMolecule] = []

    for molecule in document.molecules:
        system_size = None
        if molecule.generation:
            system_size = _safe_int(molecule.generation.system_size)

        stochastic_total = sum(
            1 for part in molecule.parts if isinstance(part, IRBigStochasticBlock)
        )

        sequences: list[StochasticSequence] = []
        descriptor_summaries: list[dict[str, object] | None] = []

        for part in molecule.parts:
            if isinstance(part, IRBigStochasticBlock):
                component = _build_stochastic_component(part, registry)
                weights = _normalize_weights(component.unit_weights)
                stored_weights = list(weights) if weights is not None else [1.0] * len(component.unit_labels)
                distribution_spec = _to_distribution_spec(component.distribution)
                chains: list[PolymerChain] = []
                target_total = system_size if stochastic_total == 1 else None

                if component.distribution and component.unit_labels:
                    sampler = _DistributionSampler(component.distribution, rng)
                    chains = _generate_polymer_chains(
                        component.unit_labels,
                        weights,
                        sampler,
                        rng,
                        target_total=target_total,
                    )
                elif target_total is not None and component.unit_labels:
                    labels = _select_units(component.unit_labels, weights, target_total, rng)
                    chains = [PolymerChain(labels=labels, degree=target_total)]

                sequences.append(
                    StochasticSequence(
                        unit_labels=component.unit_labels,
                        unit_weights=stored_weights,
                        chains=chains,
                        distribution=distribution_spec,
                        left_terminal=component.left_descriptor,
                        right_terminal=component.right_descriptor,
                        end_group_labels=component.end_group_labels,
                        target_total=target_total,
                    )
                )
            elif isinstance(part, IRBigConnector):
                descriptor_summaries.append(_descriptor_summary(part))

        molecule_results.append(
            GBigSmilesMolecule(
                sequences=sequences,
                system_size=system_size,
                bond_descriptors=descriptor_summaries,
            )
        )

    _register_fragments(document.fragments, registry)
    return GBigSmilesResult(monomers=registry.monomers, molecules=molecule_results)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


class _AtomisticBuilder:
    def __init__(self) -> None:
        self.atomistic = Atomistic()
        self._ring_state: dict[int, tuple[Atom, IRBond | None]] = {}
        self._prev_atom: Atom | None = None

    def build_molecule(self, molecule: object) -> Atomistic:
        self._ring_state.clear()
        self._prev_atom = None

        for part in getattr(molecule, "parts", []):
            if isinstance(part, IRAtom):
                bond_to_prev = part.bond_to_prev
                self._prev_atom = self._emit_atom(part, self._prev_atom, bond_to_prev)
            elif isinstance(part, (IRBigSystemDot, IRSmilesDot)):
                self._prev_atom = None
            # Ignore non-atom parts for atomistic conversion

        return self.atomistic

    def _emit_atom(
        self,
        ir_atom: IRAtom,
        parent_atom: Atom | None,
        bond_to_prev: IRBond | None,
    ) -> Atom:
        atom = self._create_atom(ir_atom)
        # register the atom in the assembly
        self.atomistic.add_entity(atom)
        if parent_atom is not None:
            self._add_bond(parent_atom, atom, bond_to_prev)
        self._handle_rings(ir_atom, atom)
        for branch in ir_atom.branches:
            self._handle_branch(atom, branch.content, branch.bond)
        return atom

    def _handle_branch(
        self,
        origin: Atom,
        content: Sequence[object],
        bond: IRBond | None,
    ) -> None:
        current_parent = origin
        last_atom: Atom | None = None
        for part in content:
            if not isinstance(part, IRAtom):
                continue
            if last_atom is None:
                connection = part.bond_to_prev or bond
                last_atom = self._emit_atom(part, current_parent, connection)
            else:
                connection = part.bond_to_prev
                last_atom = self._emit_atom(part, current_parent, connection)
            current_parent = last_atom

    def _handle_rings(self, ir_atom: IRAtom, atom: Atom) -> None:
        for ring in ir_atom.rings:
            existing = self._ring_state.get(ring.index)
            if existing is None:
                self._ring_state[ring.index] = (atom, ring.bond)
            else:
                other_atom, other_bond = existing
                bond = ring.bond or other_bond
                self._add_bond(other_atom, atom, bond)
                del self._ring_state[ring.index]

    def _create_atom(self, ir_atom: IRAtom) -> Atom:
        symbol = ir_atom.symbol or ("*" if ir_atom.kind == "star" else "X")
        element = _normalize_element(symbol)
        atom = Atom(
            symbol=symbol,
            element=element,
            kind=ir_atom.kind,
            isotope=ir_atom.isotope,
            charge=ir_atom.charge,
            hcount=ir_atom.hcount,
            chiral=ir_atom.chiral,
        )
        if ir_atom.span is not None:
            atom["span"] = {
                "start": ir_atom.span.start,
                "end": ir_atom.span.end,
                "line": ir_atom.span.line,
                "column": ir_atom.span.column,
            }
        return atom

    def _add_bond(self, left: Atom, right: Atom, bond: IRBond | None) -> None:
        kind = bond.kind if bond else "-"
        order = _bond_order(kind)
        # add bond to structure (endpoints will be tracked automatically)
        self.atomistic.add_bond(left, right, order=order, kind=kind)


def _bond_order(kind: str) -> int:
    return {
        "-": 1,
        ":": 1,
        "/": 1,
        "\\": 1,
        ";": 1,
        "=": 2,
        "#": 3,
    }.get(kind, 1)


def _normalize_element(symbol: str) -> str:
    if not symbol or symbol == "*":
        return "X"
    if len(symbol) == 1:
        return symbol.upper()
    return symbol[0].upper() + symbol[1:].lower()


def _molecule_to_atomistic(molecule: IRBigSmilesMolecule) -> Atomistic:
    builder = _AtomisticBuilder()
    builder.build_molecule(molecule)
    return builder.atomistic


def _molecule_to_monomer(molecule: IRBigSmilesMolecule) -> Monomer:
    return _atomistic_to_monomer(_molecule_to_atomistic(molecule))


def _atomistic_to_monomer(structure: Atomistic) -> Monomer:
    """Wrap a copy of the Atomistic structure into a Monomer.

    The Monomer expects an Assembly-like core; Atomistic already satisfies
    this, so we can return a shallow copy to avoid aliasing with the source
    builder state.
    """
    return Monomer(structure.copy())
