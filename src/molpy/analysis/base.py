from abc import ABC, abstractmethod

from molpy.core.logger import get_logger

logger = get_logger(__name__)


class Compute(ABC): ...
