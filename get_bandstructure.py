import json
import requests

from pymatgen import *

import numpy as np

from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter


TiCoO3 = 'mp-19424'
ID = 'mp-2352'
with MPRester(api_key="UgRqoHkuZyJEVX2d") as m:

		structure = m.get_bandstructure_by_material_id(TiCoO3)

		bands = structure.bands[Spin.up]

		for i, band in enumerate(bands):
			print(f'Band: {i}\n Kpoints: ')
			for kpoint in band:
				print(kpoint - structure.efermi, end = " ")
			print("\n")
