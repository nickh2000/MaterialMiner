from pymatgen import *
import csv
import pandas
import json 
with MPRester(api_key="UgRqoHkuZyJEVX2d") as m:
	
	
	with open('database.txt', 'w') as c:


		materials = m.query(criteria={"has_bandstructure": True}, properties=["task_id"])
		

		for material in materials:
			ID = material['task_id']
			c.write(f'{ID}\n')
		