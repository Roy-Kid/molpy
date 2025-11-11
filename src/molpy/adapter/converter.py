from __future__ import annotations

# ---- Registry -----------------------------------------------------------------
class ConverterRegistry:
    """
    A lightweight registry that stores converter callables between types.

    - Exact matches are stored in `_by_pair`: (src_type, dst_type) -> converter.
    - Fallbacks are stored in `_preds`: (predicate, dst_type, converter). These
      are evaluated in registration order when no exact MRO match is found.

    Resolution strategy:
      1) If the caller passes a `type` for `obj_or_src`, perform exact lookup.
      2) Otherwise, walk the MRO of `type(obj)` and return the first exact match.
      3) If still unresolved, try predicate-based converters whose `dst` matches.
    """

    __slots__ = ("_by_pair", "_preds")

    def __init__(self):
        # dict preserves insertion order since Python 3.7
        self._by_pair = {}   # (src_type, dst_type) -> converter
        self._preds = []     # list of (predicate, dst_type, converter)

    # --- Registration APIs (public) -------------------------------------------

    def register(self, src, dst, fn):
        """
        Register a converter for an exact (src -> dst) pair.
        Last registration wins if the same (src, dst) is re-registered.
        """
        self._by_pair[(src, dst)] = fn

    def register_predicate(self, pred, dst, fn):
        """
        Register a predicate-based converter for `dst`.
        Predicates are evaluated in registration order when no exact MRO hit exists.
        """
        self._preds.append((pred, dst, fn))

    # --- Resolution APIs (public) --------------------------------------------

    def resolve(self, obj_or_src, dst):
        """
        Resolve a converter function to `dst`.

        - If `obj_or_src` is a type, try an exact (src, dst) match.
        - If it is an instance, walk MRO for exact (mro_class, dst) matches.
        - If no exact match, try predicate-based registrations (ordered).
        """
        # Case 1: direct type provided
        if isinstance(obj_or_src, type):
            return self._by_pair.get((obj_or_src, dst))

        # Case 2: instance -> scan MRO for an exact match
        obj = obj_or_src
        for src in type(obj).mro():
            fn = self._by_pair.get((src, dst))
            if fn is not None:
                return fn

        # Case 3: try predicates in registration order
        for pred, d, fn in self._preds:
            if d is dst:
                try:
                    if pred(obj):
                        return fn
                except Exception:
                    # Swallow predicate errors to avoid blocking other fallbacks
                    pass

        return None

    def get_converter(self, src, dst):
        """
        Get an exact (src, dst) converter without MRO or predicate fallback.
        """
        return self._by_pair.get((src, dst))


# Global registry instance (intended singleton for the module)
REG = ConverterRegistry()


# ---- Unified entry: callable and subscriptable facade ------------------------
class _Convert:
    """
    Facade over REG that provides two ergonomic entry points:

      - Call style:       convert(obj, ToType)
      - Subscript style:  convert[ToType](obj)

    Both delegate to `REG.resolve`. Errors are raised only when invoked,
    never during subscription (so `convert[To]` itself is safe to obtain).
    """

    __slots__ = ()

    def __call__(self, obj, to):
        fn = REG.resolve(obj, to)
        if fn is None:
            # Build a helpful message while preserving original API/semantics
            pairs_preview = ", ".join(
                f"{s.__name__}->{d.__name__}" for (s, d) in REG._by_pair.keys()
            )
            raise TypeError(
                f"No converter registered for {type(obj).__name__} -> {to.__name__}. "
                f"Available pairs: [{pairs_preview}]"
            )
        return fn(obj)

    def __getitem__(self, dst):
        # Returns a thin wrapper that defers resolution until invocation time.
        def _f(obj):
            return self(obj, dst)
        return _f


# Public singleton (API must remain the same)
convert = _Convert()


# ---- Decorator sugar ---------------------------------------------------------
def register(src, dst):
    """
    Decorator: register a concrete (src -> dst) converter.

        @register(Source, Target)
        def to_target(x: Source) -> Target: ...
    """
    def _wrap(fn):
        REG.register(src, dst, fn)
        return fn
    return _wrap


def register_if(pred, dst):
    """
    Decorator: register a predicate-based converter for `dst`.

        @register_if(lambda x: condition(x), Target)
        def to_target(x) -> Target: ...
    """
    def _wrap(fn):
        REG.register_predicate(pred, dst, fn)
        return fn
    return _wrap
