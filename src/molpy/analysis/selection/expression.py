from typing import Callable, Optional

import numpy as np

from ..base import Compute, ComputeContext


class SelectResult(ComputeContext):
    """
    Result of a selection operation with chain-based access.

    Each SelectResult wraps the previous context, forming a chain.
    The actual selection data is stored in the selected attribute.
    """

    selected: Optional[np.ndarray]

    def __post_init__(self):
        """
        Initialize the SelectResult with an empty result dictionary.
        """
        super().__post_init__()
        self.selected = None


class ExpressionSelection(Compute):
    """
    Selection operation based on a boolean expression.
    """

    def __init__(
        self,
        name: str,
        expression: Callable[[ComputeContext], np.ndarray],
        select_field: str = "atoms",
    ):
        super().__init__(name)
        self.expression = expression
        self.select_field = select_field

    def compute(self, context: ComputeContext) -> ComputeContext:
        """
        Apply selection expression and return SelectResult.

        Args:
            context: Input compute context

        Returns:
            SelectResult wrapping the context with selection applied
        """
        # Get boolean mask from expression
        mask = self.expression(context)

        selected_ids = context.frame[self.select_field]["id"][mask]

        # Create SelectResult with the new selection
        result = SelectResult(context)
        result.selected = selected_ids
        return result
