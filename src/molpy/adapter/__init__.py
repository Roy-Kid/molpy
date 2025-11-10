try:
    import rdkit
    from .rdkit import to_rdkit, generate_3d_coords, draw_molecule
except ImportError:
    pass    