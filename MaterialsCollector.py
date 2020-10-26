from pymatgen import *
import random
from DosData import *
import requests
import json
import pandas as pd
from sqlalchemy import create_engine
import sqlite3

test_upper_gap = 1.3
test_lower_gap = -1.1
test_width = .35

API_KEY = 'UgRqoHkuZyJEVX2d'


def store_all_material_ids_in_database():
	with MPRester(API_KEY) as m:	
		with open('database.txt', 'w') as c:

			#retreve and store all known materials with bandstructures by task ID
			materials = m.query(criteria={"has_bandstructure": True}, properties=["task_id"])
			
			for material in materials:
				ID = material['task_id']
				c.write(f'{ID}\n')

def find_candidates_by_material_id(num):
	#go through materials with known DOS's
	with open('database.txt', 'r')as d:
		#where the candidate materials will be stored
		with open('candidates.txt', 'w')  as c:

			#get first 1000 materials
			materials = [next(d) for x in range(10000)]
			random.shuffle(materials)

			for ID in materials[:num]:
				ID = ID.rstrip()
				DOS = DosData(ID)
				#store the material if it is valid
				if DOS.is_valid(test_upper_gap, test_upper_gap, test_width):
					c.write(f'{ID}\n')

def analyze_candidates():
	with open('candidates.txt', 'r') as c:
		with MPRester(API_KEY) as m:

			materials = c.readlines()

			for ID in materials:
				ID = ID.rstrip()
				DOS = DosData(ID)
				
				centered_array  = DOS.center_dos(10, -10)
				plt.plot(list(centered_array.values()), list(centered_array.keys()))
				plt.ylabel(DOS.get_formula())
				plt.show()

def compile_dos():
	with open('database.txt', 'r') as d:
		materials = d.readlines()
		chunk = []
		for i, ID in enumerate(materials):
			ID = ID.rstrip()
			r = requests.get(f'https://www.materialsproject.org/rest/v2/materials/{ID}/vasp/dos?API_KEY={API_KEY}')
			r = r.json()['response'][0]['dos']
			energies =list(r['energies'])
			densities = list(r['densities']['1'])
			efermi = r['efermi']
			chunk.append([ID, energies, densities, efermi])
			if not i % 1000:
				df = pd.DataFrame(chunk, columns = ['ID', 'Energies', 'Densities', 'efermi'])
				df.set_index("ID")
				print(f'Computed entries: {i}')
				chunk = []
				df.to_csv('dos.csv')

if __name__ == '__main__':
	compile_dos()
	