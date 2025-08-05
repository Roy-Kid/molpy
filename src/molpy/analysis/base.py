from abc import ABC, abstractmethod
from typing import Any, Self
import molpy as mp

from molpy.core.wrapper import Wrapper
from molpy.core.frame import Frame
from molpy.core.logger import get_logger

logger = get_logger(__name__)

class ComputeContext(Wrapper[Frame]):
    """
    Context for computing analysis.
    """

    def __init__(self, frame: Frame):
        super().__init__(frame)
        self.result = dict()

    @property
    def frame(self) -> Frame:
        return self.unwrap()


class Compute(ABC):

    def __init__(self, name: str):
        self.name = name

    @abstractmethod
    def compute(self, context: ComputeContext) -> ComputeContext:
        ...
    
    def __call__(self, context: ComputeContext | Frame) -> ComputeContext:

        context = self._ensure_context(context)

        self._check_inputs(context)

        self.compute(context)

        self._check_outputs(context)

        return context

    def _ensure_context(self, context: ComputeContext | Frame) -> ComputeContext:
        if isinstance(context, Frame):
            return ComputeContext(context)
        return context

    def _check_inputs(self, _: ComputeContext) -> Exception:
        logger.debug(f"Checking inputs for {self}")

    def _check_outputs(self, _: ComputeContext) -> Exception:
        logger.debug(f"Checking outputs for {self}")

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} {self.name}>"
