from functools import partial, wraps
from typing import Any, Callable, Tuple

import numpy as np

import molpy as mp

from ..base import Potential


class PairPotential(Potential):
    """
    Base class for pair potentials.

    Provides a decorator to allow methods that accept precomputed pair displacements
    (dr, dr_norm, pair_types) to also accept an mp.Frame directly.
    """

    @classmethod
    def or_frame(cls, *extra_fields: str) -> Callable:
        """
        Decorator that enables a method to accept either (dr, dr_norm, pair_types[, extra...])
        or a single mp.Frame argument followed by optional args.

        extra_fields: column names in frame to extract before calling the method.
        """

        def decorator(func: Callable) -> Callable:
            @wraps(func)
            def wrapper(self, *args: Any, **kwargs: Any):
                # Case: first arg is an mp.Frame
                if args and isinstance(args[0], mp.Frame):
                    frame = args[0]
                    # extract coordinates
                    coords = frame["atoms", "xyz"]
                    # external neighborlist & periodic corrections assumed done
                    # extract pair info
                    pairs = frame["pairs"]
                    dr = pairs["dr"]
                    dr_norm = np.linalg.norm(dr, axis=1, keepdims=True)
                    pair_types = pairs["type"]
                    # extract extras
                    extra_args = [frame[field] for field in extra_fields]
                    return func(self, dr, dr_norm, pair_types, *extra_args, **kwargs)
                # otherwise, call as regular method
                return func(self, *args, **kwargs)

            return wrapper

        # allow use as @or_frame or @or_frame('field1', 'field2')
        if len(extra_fields) == 1 and callable(extra_fields[0]):  # used without args
            func = extra_fields[0]
            extra_fields = ()
            return decorator(func)
        return decorator
