from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter
from pymatgen import *
import numpy as np
from matplotlib import pyplot as plt
import random
from DosData import *


TiCoO3 = 'mp-19424'


if __name__ == "__main__":

	dos = DosData(TiCoO3)
	
	bounds = np.array(dos.get_all_bounds())
	
	max_width = bounds[np.argmax(np.abs(bounds[:2]))]

	min_gap = bounds[np.argmin(np.abs(bounds[2:])) + 2]

	print(f'Width: {max_width}, Bound: {min_gap}')
