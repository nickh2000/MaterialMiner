from pymatgen import *
import random
from DosData import *
import requests
import pandas as pd
import itertools
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

#plot the DOS of each material in candidates.txt
def analyze_candidates():
	with open('candidates.txt', 'r') as c:
		materials = c.readlines()

		for ID in materials:
			ID = ID.rstrip()
			DOS = DosData(ID)
			
			centered_array  = DOS.center_dos(10, -10)
			plt.plot(list(centered_array.values()), list(centered_array.keys()))
			plt.ylabel(DOS.get_formula())
			plt.show()



#function that collects and stores the DOS of all materials in a CSV file
def compile_dos():
	with open('database.txt', 'r') as d:
		#get all 70k materials from the database
		materials = d.readlines()

		#chunk of pandas DataFrames
		chunk = []

		#iterate through each material ID
		for i, ID in enumerate(materials[41001:]):
			ID = ID.rstrip()

			#Get the DOS data for that material
			r = requests.get(f'https://www.materialsproject.org/rest/v2/materials/{ID}/vasp/dos?API_KEY={API_KEY}')
			r = r.json()['response'][0]['dos']

			#get different properties of the DOS array
			energies =list(r['energies'])
			densities = list(r['densities']['1'])
			efermi = r['efermi']

			#add the DOS information to the current chunk
			chunk.append([ID, energies, densities, efermi])

			#every thousand entries, store the chunk
			if not i % 1000 and i != 0:

				#convert array of entries into a Dataframe
				df = pd.DataFrame(chunk, columns = ['ID', 'Energies', 'Densities', 'efermi'])
				#Print progress
				print(f'Computed entries: {i}')
				#Clear the current chunk
				chunk = []
				#store the DataFrame
				df.to_csv('dos.csv', mode='a', header=False)


def plot_database(lower = 0, upper = 5):

	#values representing energies rounded to the nearest precision value
	energies = np.arange(0, max_bound, precision)
	num_energies = len(energies)

	data = np.zeros((num_energies, num_energies))

	#width and gap energies respectively
	X, Y = np.meshgrid(energies, energies)

	#size of chunks gathered from CSV file
	chunksize = 5000

	#open the CSV file into chunks
	reader = pd.read_csv('dos.csv', chunksize=chunksize)

	#process each chunk
	for c, chunk in enumerate(itertools.islice(reader, 10)):
		#process each data entry in a chunk (a pd.Series row object)
		for r, (index, row) in enumerate(chunk.iterrows()):
			#output every 1000 entries
			if not r%1000:
				print(f'{c*chunksize + r} entries processed')

			#create a DOS object from the row
			dos = DosData(row)

			#find the maximum width and the minimum gap of the DOS
			width, gap = dos.get_parameters()

			#increment the number of entries with the nearest width and gap
			data[round(width / precision), round(gap / precision)] += 1

	#start a new 3D figure
	fig = plt.figure()
	ax = Axes3D(fig)
	ax.set_xlabel('Gap')
	ax.set_ylabel('Spaghetti band width')

	#start energy index 
	start = round(lower/precision)

	#set the axis bounds
	ax.axes.set_xlim(right = 5, left = lower)
	ax.axes.set_ylim(bottom = 5, top = 0)

	#cut off data below the lower bound
	ax.plot_surface(X[:, start:], Y[:, start:], data[:, start:])
	plt.show()

def store_candidates():
	#where the candidate materials will be stored
	with open('candidates.txt', 'w')  as f:
		#size of chunks gathered from CSV file
		chunksize = 5000

		#open the CSV file into chunks
		reader = pd.read_csv('dos.csv', chunksize=chunksize)

		#process each chunk
		for c, chunk in enumerate(itertools.islice(reader, 10)):
			#process each data entry in a chunk (a pd.Series row object)
			for r, (index, row) in enumerate(chunk.iterrows()):
				dos = DosData(row)

				if dos.is_valid(test_upper_gap, test_lower_gap, test_width):
					f.write(f'{dos.material_id}\n')

				if not r%1000:
					print(f'{c*chunksize + r} entries processed')

def create_training_data():

	#setup headers
	df = pd.DataFrame(columns = ['ID', 'Formula', 'Is_Candidate'])
	df.to_csv('training_data.csv')

	#size of chunks gathered from CSV file
	chunksize = 5000

	#open the CSV file into chunks
	reader = pd.read_csv('dos.csv', chunksize=chunksize)

	#process each chunk
	for c, chunk in enumerate(itertools.islice(reader, 10)):

		formulas = MPRester(API_KEY).query(criteria={'task_id': {'$in': list(chunk['ID'])}}, properties=['pretty_formula', 'task_id'])
		currentChunk = []

		for i, entry in enumerate(formulas):
			row = chunk.loc[chunk['ID'] == entry['task_id']]
			newEntry = {}
			newEntry['ID'] = entry['task_id']
			newEntry['Formula'] = entry['pretty_formula']
			newEntry['Is_Candidate'] = DosData(row.squeeze()).is_valid(test_upper_gap, test_lower_gap, test_width)
			currentChunk.append(newEntry)

			if not i%100:
				print(f'{c*chunksize + i} Computeted entries')

		df = pd.DataFrame(currentChunk, columns = ['ID', 'Formula', 'Is_Candidate'])
		df.to_csv('training_data.csv', mode='a', header=False)



if __name__ == '__main__':
	create_training_data()