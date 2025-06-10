import molpy as mp
from typing import Literal
import random
from .target import Target
from .packer import get_packer

class Session:

    def __init__(self, packer: Literal["packmol", "nlopt"] = "packmol"):
        
        self.targets = []
        self.packer = get_packer(packer)

    def add_target(self, frame:mp.Frame, number:int, constrain:mp.Region):
        self.targets.append(Target(frame, number, constrain))

    def optimize(self, max_steps: int = 1000, seed: int|None = None):
        if seed is None:
            seed = random.randint(1, 10000)
        return self.packer.pack(self.targets, max_steps=max_steps, seed=seed)
