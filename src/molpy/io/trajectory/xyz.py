from .base import TrajectoryWriter
from molpy.core import Frame
import xarray as xr

class XYZTrajectoryWriter(TrajectoryWriter):

    def __init__(self, fpath: str):
        
        self.fpath = fpath
        self.fobj = open(fpath, "w")

    def __del__(self):
        if not self.fobj.closed:
            self.fobj.close()

    def write_frame(self, frame: Frame):

        atoms = frame["atoms"]
        if frame.box is None:
            box = frame.box
        n_atoms = len(atoms)

        self.fobj.write(f"{n_atoms}\n")
        self.fobj.write(f"Step={frame.get('step')} Lattic\"{box.matrix.tolist()}\"\n")

        n = atoms.sizes.get("index", 0)
        elem = atoms.get("element", xr.DataArray(["X"] * n))
        for i in range(n):
            x = atoms["x"].values[i]
            y = atoms["y"].values[i]
            z = atoms["z"].values[i]
            e = elem.values[i]
            self.fobj.write(f"{e} {x} {y} {z}\n")

    def write_traj(self, trajectory):

        for frame in trajectory:
            self.write_frame(frame)

    def close(self):

        self.fobj.close()