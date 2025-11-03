"""Lightweight core package initializer for the refactored architecture.

This module intentionally avoids eager imports to prevent import-time side
effects across legacy modules. New components should be imported directly from
their submodules, e.g.:

	from molpy.core.entity import Entity
	from molpy.core.link import Link
	from molpy.core.assembly import Assembly

Legacy wildcard exports have been removed as part of the refactor.
"""

__all__: list[str] = []
