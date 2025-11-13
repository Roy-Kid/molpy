import hashlib
import inspect
from abc import ABC, abstractmethod
from functools import lru_cache
from typing import Any, Generic, Iterable, ParamSpec, TypeVar

from molpy.core.logger import get_logger

logger = get_logger(__name__)

P = ParamSpec("P")
T = TypeVar("T")
R = TypeVar("R", covariant=True)


class ComputeMeta:
    """Compute type metadata."""

    def __init__(
        self,
        id: str,
        name: str,
        module: str,
        qualname: str,
        signature: str,
        code_digest: str,
        registered_at: float = 0.0,
    ):
        self.id = id
        self.name = name
        self.module = module
        self.qualname = qualname
        self.signature = signature
        self.code_digest = code_digest
        self.registered_at = registered_at


@lru_cache(maxsize=None)
def _sig_info(func) -> tuple[int, bool, frozenset]:
    """Get the signature information of a function."""
    func = getattr(func, "forward")
    sig = inspect.signature(func)
    params = list(sig.parameters.values())[1:]
    P = inspect.Parameter
    pos_params = [
        p for p in params if p.kind in (P.POSITIONAL_ONLY, P.POSITIONAL_OR_KEYWORD)
    ]
    has_var_pos = any(p.kind == P.VAR_POSITIONAL for p in params)
    kw_names = frozenset(
        p.name for p in params if p.kind in (P.POSITIONAL_OR_KEYWORD, P.KEYWORD_ONLY)
    )
    return len(pos_params), has_var_pos, kw_names


def _route_args_kwargs(
    compute, args: tuple, kwargs: dict[str, Any], *, pos_overflow: str = "ignore"
) -> tuple[tuple, dict[str, Any]]:
    """Crop the actual parameters accepted by compute.forward. on_drop: ignore | warn | error"""
    max_pos, has_var_pos, kw_names = _sig_info(compute.__class__)

    # positional argument processing
    if has_var_pos:
        pass_args = args
    else:
        if len(args) > max_pos:
            extras = args[max_pos:]
            if pos_overflow == "error":
                raise TypeError(
                    f"{compute.__class__.__name__}.forward got {len(args)} positional args, "
                    f"but accepts only {max_pos}. Extras={extras!r}"
                )
            elif pos_overflow == "warn":
                import warnings

                warnings.warn(
                    f"{compute.__class__.__name__}.forward received extra positional args; dropping {extras!r}"
                )
        pass_args = args[:max_pos]

    # kwargs: always only retain explicitly declared keys (even if the target has **kwargs)
    pass_kwargs = {k: v for k, v in kwargs.items() if k in kw_names}

    return pass_args, pass_kwargs


def generate_compute_id(module: str, qualname: str, code_digest: str) -> str:
    """Generate unique compute ID from module, qualname, and code digest."""
    content = f"{module}:{qualname}:{code_digest}"
    return hashlib.sha1(content.encode("utf-8")).hexdigest()


class Compute(ABC, Generic[P, R]):
    """
    Base class for molecular analysis computations that implements the Op protocol.

    This class provides a unified interface for all molecular analysis
    operations, implementing the same protocol as workflow operators.
    """

    def __init__(self, *deps: "Compute") -> None:
        """Initialize the compute with dependencies.

        Args:
            deps: Static (hard) dependencies used for static planning/topology.
        """
        self._deps: tuple[Compute, ...] = tuple(deps)

    @property
    def parents(self) -> Iterable["Compute"]:
        """Return static (hard) dependencies only—used for static planning/topology."""
        return self._deps

    @abstractmethod
    def forward(self, *args: P.args, **kwargs: P.kwargs) -> R:
        """
        Compute this operation's value.
        Subclasses MUST implement this.

        Args:
            *args: Input arguments for the computation.
            **kwargs: Keyword arguments for the computation.
        Returns:
            The computed value of this operation.
        """
        raise NotImplementedError

    def compute(self, *args: P.args, **kwargs: P.kwargs) -> R:
        """
        Alias for forward() method for backward compatibility.

        Args:
            *args: Input arguments for the computation.
            **kwargs: Keyword arguments for the computation.
        Returns:
            The computed value of this operation.
        """
        return self.forward(*args, **kwargs)

    def __repr__(self) -> str:
        """Return string representation of the compute operation."""
        return f"{self.__class__.__name__}()"

    def clone(self) -> "Compute":
        """Clone the compute operation.

        Returns:
            A new instance of the same compute type with the same dependencies.
        """
        return self.__class__(*self._deps)

    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> R:
        """Execute the compute operation with the given arguments.

        This method handles argument routing and can be extended for tracing integration.

        Args:
            *args: Positional arguments for the computation.
            **kwargs: Keyword arguments for the computation.

        Returns:
            The result of the compute operation's forward method.
        """
        # Route arguments to match the compute's signature
        pass_args, pass_kwargs = _route_args_kwargs(self, args, kwargs)

        # Execute the forward method
        return self.forward(*pass_args, **pass_kwargs)

    def get_compute_id(self) -> str:
        """Get unique identifier for this compute type.

        Returns:
            A unique string identifier for this compute type.
        """
        return generate_compute_id(
            module=self.__class__.__module__,
            qualname=self.__class__.__qualname__,
            code_digest=self.get_code_digest(),
        )

    def get_code_digest(self) -> str:
        """Get code digest for this compute type.

        Returns:
            A hash representing the compute's code.
        """
        source = inspect.getsource(self.__class__)
        return hashlib.sha256(source.encode("utf-8")).hexdigest()

    def get_signature(self) -> str:
        """Get string representation of the compute's signature.

        Returns:
            A string representation of the forward method's signature.
        """
        sig = inspect.signature(self.forward)
        return str(sig)

    def get_meta(self) -> ComputeMeta:
        """Get metadata for this compute type.

        Returns:
            ComputeMeta object containing compute metadata.
        """
        return ComputeMeta(
            id=self.get_compute_id(),
            name=self.__class__.__name__,
            module=self.__class__.__module__,
            qualname=self.__class__.__qualname__,
            signature=self.get_signature(),
            code_digest=self.get_code_digest(),
            registered_at=0.0,  # Will be set by the registry
        )


class NullaryCompute(Compute):
    """Nullary compute base class.

    A nullary compute takes no dependencies and operates on input data directly.
    This is the base class for leaf nodes in the computation DAG.
    """

    def __init__(self) -> None:
        """Initialize a nullary compute with no dependencies."""
        super().__init__()


class UnaryCompute(Compute):
    """Unary compute base class.

    A unary compute takes exactly one dependency and applies a transformation
    to its output. This is useful for single-input operations like scaling,
    filtering, or any other unary transformation.
    """

    def __init__(self, dep: Compute) -> None:
        """Initialize the unary compute with a dependency.

        Args:
            dep: The single dependency compute.
        """
        super().__init__(dep)

    @property
    def dep(self) -> Compute:
        """Get the dependency compute.

        Returns:
            The single dependency compute.
        """
        return self._deps[0]


class BinaryCompute(Compute):
    """Binary compute base class.

    A binary compute takes exactly two dependencies and combines their outputs
    using some binary operation. This is useful for operations like addition,
    multiplication, comparison, etc.
    """

    def __init__(self, dep1: Compute, dep2: Compute) -> None:
        """Initialize the binary compute with two dependencies.

        Args:
            dep1: The first dependency compute.
            dep2: The second dependency compute.
        """
        super().__init__(dep1, dep2)

    @property
    def dep1(self) -> Compute:
        """Get the first dependency compute.

        Returns:
            The first dependency compute.
        """
        return self._deps[0]

    @property
    def dep2(self) -> Compute:
        """Get the second dependency compute.

        Returns:
            The second dependency compute.
        """
        return self._deps[1]
