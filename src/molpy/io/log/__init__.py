"""Simulation log readers: the molrs structured LAMMPS log.

``read_lammps_log(path)`` returns a :class:`LammpsLog` — ``header``, one
:class:`LammpsRun` per ``run`` (``thermo["Step"]`` is a float64 column),
``warnings`` — parsed in Rust; ``to_dict()`` gives the JSON form.
"""

from molrs.io import (
    LammpsCpuUse,
    LammpsLoadBalance,
    LammpsLog,
    LammpsLogHeader,
    LammpsLoopTime,
    LammpsMemoryUsage,
    LammpsNeighborStatistics,
    LammpsPerformance,
    LammpsRun,
    LammpsThermo,
    LammpsTimingBreakdown,
    LammpsTimingRow,
    LammpsWarning,
    parse_lammps_log_text,
    read_lammps_log,
)

__all__ = [
    "LammpsCpuUse",
    "LammpsLoadBalance",
    "LammpsLog",
    "LammpsLogHeader",
    "LammpsLoopTime",
    "LammpsMemoryUsage",
    "LammpsNeighborStatistics",
    "LammpsPerformance",
    "LammpsRun",
    "LammpsThermo",
    "LammpsTimingBreakdown",
    "LammpsTimingRow",
    "LammpsWarning",
    "parse_lammps_log_text",
    "read_lammps_log",
]
