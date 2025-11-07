try:
    import rdkit
    from .rdkit import to_rdkit, from_rdkit
except ImportError:
    pass    