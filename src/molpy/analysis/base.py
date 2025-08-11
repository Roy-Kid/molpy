from abc import ABC, abstractmethod
from typing import Generic, Optional, TypeVar

import numpy as np

from molpy.core.frame import Frame
from molpy.core.logger import get_logger

logger = get_logger(__name__)

T = TypeVar("T", bound="ComputeContext")


class ComputeContext(Generic[T]):
    """
    Context for computing analysis.

    Can wrap another ComputeContext or be the root context with a Frame.
    Provides chain-based access to traverse context history.
    """

    def __init__(self, context: Optional["ComputeContext"] = None):
        self._wrapped_context = context
        self._frame: Frame | None = None
        self.__post_init__()

    def __post_init__(self):
        """
        Initialize the ComputeContext with an empty result dictionary.
        This method can be overridden in subclasses for custom initialization.
        """
        self.result = dict()

    @classmethod
    def attach_frame(cls, frame: Frame) -> "ComputeContext":
        """Create a root ComputeContext attached to a Frame. Avoid mutating Frame directly."""
        context = cls()
        context._frame = frame
        return context

    @property
    def frame(self) -> Frame:
        """Get the Frame, traversing up the context chain if necessary."""
        if self._frame is not None:
            return self._frame
        elif self._wrapped_context is not None:
            return self._wrapped_context.frame
        else:
            raise ValueError("No Frame found in context chain")

    def unwrap(self) -> Optional["ComputeContext"]:
        """Get the wrapped context."""
        return self._wrapped_context

    def get_head(self) -> "ComputeContext":
        """Get the head (topmost) context in the chain."""
        current = self
        while current.unwrap() is not None:
            unwrapped = current.unwrap()
            if unwrapped is not None:
                current = unwrapped
            else:
                break
        return current

    def pop(self) -> "ComputeContext":
        """
        Pop the current context and return the previous context.

        Returns:
            The previous ComputeContext in the chain

        Raises:
            ValueError: If no previous context found
        """
        previous = self.unwrap()
        if previous is None:
            raise ValueError("No previous context found")
        return previous

    def get_depth(self) -> int:
        """Get the depth of the context chain."""
        depth = 0
        current = self
        while current.unwrap() is not None:
            depth += 1
            unwrapped = current.unwrap()
            if unwrapped is not None:
                current = unwrapped
            else:
                break
        return depth

    def get_stack(self) -> list["ComputeContext"]:
        """Get all contexts in the chain from current to root."""
        contexts = []
        current = self
        while current is not None:
            contexts.append(current)
            current = current.unwrap()
        return contexts


class Compute(ABC):

    def __init__(self, name: str):
        self.name = name

    @abstractmethod
    def compute(self, context: ComputeContext) -> ComputeContext: ...

    def __call__(self, context: ComputeContext | Frame) -> ComputeContext:

        context = self._ensure_context(context)

        self._check_inputs(context)

        result = self.compute(context)

        self._check_outputs(result)

        return result

    def map(
        self,
        context: ComputeContext,
        entity_ids: Optional[np.ndarray] = None,
        selection_field: str = "molecules",
    ) -> ComputeContext:
        """
        Apply this computation to each entity in a map-reduce fashion.

        Args:
            context: Input context
            entity_ids: Optional array of entity IDs to map over. If None, will use all entities from selection_field
            selection_field: Field containing the entities to map over (used when entity_ids is None)

        Returns:
            Context with mapped results in result['{field}_results']
        """
        # Get entity IDs to map over
        if entity_ids is None:
            entity_ids = self._get_all_entities(context, selection_field)

        if entity_ids is None or len(entity_ids) == 0:
            raise ValueError(f"No entities found for mapping")

        # Apply computation to each entity
        mapped_results = {}
        for entity_id in entity_ids:
            # Create entity-specific context
            entity_context = self._create_entity_context(
                context, entity_id, selection_field
            )

            # Apply computation
            entity_result = self.compute(entity_context)

            # Store result
            mapped_results[entity_id] = entity_result.result

        # Create result context
        result_context = ComputeContext(context)
        result_context.result[f"{selection_field}_results"] = mapped_results
        result_context.result[f"mapped_{selection_field}"] = entity_ids

        return result_context

    def _get_all_entities(
        self, context: ComputeContext, selection_field: str
    ) -> Optional[np.ndarray]:
        """Get all entities from the specified field in frame."""
        if selection_field in context.frame:
            return context.frame[selection_field]["id"]
        return None

    def _create_entity_context(
        self, base_context: ComputeContext, entity_id: int, selection_field: str
    ) -> ComputeContext:
        """Create a context for a single entity."""
        entity_context = ComputeContext(base_context)
        entity_context.result["entity_id"] = entity_id
        entity_context.result["entity_field"] = selection_field
        return entity_context

    def _ensure_context(self, context: ComputeContext | Frame) -> ComputeContext:
        if isinstance(context, Frame):
            return ComputeContext.attach_frame(context)
        return context

    def _check_inputs(self, _: ComputeContext) -> None:
        logger.debug("Checking inputs for %s", self)

    def _check_outputs(self, _: ComputeContext) -> None:
        logger.debug("Checking outputs for %s", self)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} {self.name}>"
