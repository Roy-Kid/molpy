from typing import Callable
from ..base import Compute, ComputeContext

class ExpressionSelection(Compute):

    def __init__(self, name: str, expression: Callable[ComputeContext, np.ndarray]):
        super().__init__(name)
        self.expression = expression

    def compute(self, context: ComputeContext) -> ComputeContext:

        mask = self.expression(context)
        
