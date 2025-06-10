from pathlib import Path

import numpy as np
import molpy as mp
from molpy.core.arraydict import ArrayDict

class AmberInpcrdReader:

    def __init__(self, file: str | Path, ):
        self._file = Path(file)

    @staticmethod
    def sanitizer(line: str) -> str:
        return line
    
    def read(self, frame: mp.Frame):

        with open(self._file, "r") as f:

            lines = filter(
                lambda line: line,
                map(AmberInpcrdReader.sanitizer, f),
            )

            title = next(lines).strip()

            num_atoms = int(next(lines))

            coordinates = []
            for line in lines:
                values = [float(line[i-12:i].strip()) for i in range(12, len(line), 12)]
                coordinates.extend(values)

            x_coords = coordinates[0::3]
            y_coords = coordinates[1::3]
            z_coords = coordinates[2::3]

            table = ArrayDict({
                'id': np.arange(num_atoms) + 1,
                'x': x_coords,
                'y': y_coords,
                'z': z_coords
            })

        frame['props']['name'] = title
        return frame