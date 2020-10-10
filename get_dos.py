import json
import requests

from pymatgen import *

import numpy as np

from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter
import numpy as np





def get_dos_tuples(material_id):
	with MPRester(api_key="UgRqoHkuZyJEVX2d") as m:

		dos = m.get_dos_by_material_id(TiCoO3)
		energies = np.array(dos.as_dict()['energies'])
		e_fermi = dos.as_dict()['efermi']
		densities = np.array(dos.get_densities())

	return zip(energies - e_fermi, densities):