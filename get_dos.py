import json
import requests

from pymatgen import *

import numpy as np

from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter
import numpy as np

TiCoO3 = 'mp-19424'

with MPRester(api_key="UgRqoHkuZyJEVX2d") as m:

		dos = m.get_dos_by_material_id(TiCoO3)

		energies = np.array(dos.as_dict()['energies'])
		e_fermi = dos.as_dict()['efermi']
		densities = np.array(dos.get_densities())

		for density in zip(energies - e_fermi, densities):
				print(density)


		'''dos_plotter = DosPlotter()
		dos_plotter.add_dos_dict(dos.get_spd_dos())
		dos_plotter.show()'''