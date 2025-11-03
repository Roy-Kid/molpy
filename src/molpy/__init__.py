"""Lightweight molpy package initializer for refactored core tests.

Avoid eager imports to prevent legacy dependencies from loading during
unit tests that target the new core architecture.
"""

from .version import version  # noqa: F401

__all__: list[str] = ["version"]
