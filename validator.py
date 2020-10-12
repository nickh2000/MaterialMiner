
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter
from pymatgen import *
import numpy as np
from matplotlib import pyplot as plt
import random
import get_candidates


if __name__ == "__main__":
	with open('candidates.txt', 'r') as c:
		with MPRester('UgRqoHkuZyJEVX2d') as m:
			for ID in c.readlines():
				ID = ID.rstrip()


				name = m.query(criteria={'task_id':ID}, properties=['pretty_formula'])[0]['pretty_formula']
				dos = get_candidates.center_dos(get_candidates.get_dos_array(ID))
				plt.plot(list(dos.values()), list(dos.keys()))
				plt.ylabel(name)
				plt.show()