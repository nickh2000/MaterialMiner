from pymatgen import *

import numpy as np


with MPRester(api_key="UgRqoHkuZyJEVX2d") as m:

	
	
	with open("candidates.txt", "w") as c:

		materials = m.query(criteria={"has_bandstructure": True}, properties=["task_id"])
		


		for material in materials:
			ID = material['task_id']
			c.write(f'{ID}\n')