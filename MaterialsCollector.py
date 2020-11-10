from pymatgen import *
import random
from DosData import *
import requests
import pandas as pd
import itertools
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import periodictable 
import json
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.sampledata.periodic_table import elements
from bokeh.transform import dodge, linear_cmap, log_cmap
import bokeh.palettes
from bokeh.models import ColorBar
import math

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
		for i, ID in enumerate(materials[69001:]):
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
		for c, chunk in enumerate(reader):
			#process each data entry in a chunk (a pd.Series row object)
			for r, (index, row) in enumerate(chunk.iterrows()):
				dos = DosData(row)

				for offset in np.arange(0, -5, -.1):
					#test for thin-width at various energy levels
					if dos.has_isolated_band(test_upper_gap, test_lower_gap, test_width, offset):
						f.write(f'{dos.material_id}\n')
						break

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



def elements_in_candidates():
	elementDict = {}
	elementDict = {str(element): 0 for element in periodictable.elements}
	elementDict['NumCandidates'] = 0
	with open("candidates.txt") as f:
		materials = [material.rstrip() for material in f.readlines()]
		response = MPRester(API_KEY).query(criteria={'task_id': {'$in': list(materials)}}, properties=['elements'])

		for material in response:
			elements = material['elements']
			for element in elements:
				elementDict[element] += 1
			elementDict['NumCandidates'] += 1

	json.dump(elementDict, open('elements.json', 'w+'), indent = 4)


def plot_periodic_table():
	output_file("periodic.html")

	occurences = json.load(open('elements.json', 'r'))

	periods = ["I", "II", "III", "IV", "V", "VI", "VII"]
	groups = [str(x) for x in range(1, 19)]

	df = elements.copy()
	df["atomic mass"] = df["atomic mass"].astype(str)
	df["group"] = df["group"].astype(str)
	df["period"] = [periods[x-1] for x in df.period]
	df = df[df.group != "-"]
	df = df[df.symbol != "Lr"]
	df = df[df.symbol != "Lu"]
	df['occurences'] = [f'{round(occurences[x]/occurences["NumCandidates"] * 100, 2)}%' for x in df.symbol]
	df['colors'] = [math.log(occurences[x] / 2 + 1) for x in df.symbol]
	
	p = figure(title="Periodic Table (omitting LA and AC Series)", plot_width=1000, plot_height=450,
	           x_range=groups, y_range=list(reversed(periods)))

	#Use the field name of the column source
	mapper = linear_cmap(field_name='colors', palette=bokeh.palettes.Spectral6, low = 1, high = max(df.colors))

	r = p.rect("group", "period", 0.95, 0.95, source=df, fill_alpha=0.6,
	           color=mapper)

	color_bar = ColorBar(color_mapper=mapper['transform'], width=8,  location=(0,0))

	p.add_layout(color_bar, 'right')


	text_props = {"source": df, "text_align": "left", "text_baseline": "middle"}

	x = dodge("group", -0.4, range=p.x_range)

	p.text(x=x, y="period", text="symbol", text_font_style="bold", **text_props)

	p.text(x=x, y=dodge("period", 0.3, range=p.y_range), text="atomic number",
	       text_font_size="11px", **text_props)

	p.text(x=x, y=dodge("period", -0.35, range=p.y_range), text="name",
	       text_font_size="7px", **text_props)

	p.text(x=x, y=dodge("period", -0.2, range=p.y_range), text="occurences",
	       text_font_size="7px", **text_props)

	p.text(x=["3", "3"], y=["VI", "VII"], text=["LA", "AC"], text_align="center", text_baseline="middle")

	p.outline_line_color = None
	p.grid.grid_line_color = None
	p.axis.axis_line_color = None
	p.axis.major_tick_line_color = None
	p.axis.major_label_standoff = 0
	p.legend.orientation = "horizontal"
	p.legend.location ="top_center"
	p.hover.renderers = [r] # only hover element boxes

	show(p)

def clean_dos():
	with open('database.txt') as d:

		reader = pd.read_csv('dos.csv', chunksize = 5000)

		materials = [m.rstrip() for m in d.readlines()]
		missing = {m: True for m in materials}
		

		for chunk in reader:
			for i, material in enumerate(materials):

				if missing[material]:
					missing[material] = not (material in list(chunk['ID']))
		
		missingMaterials = []
		for k, v in missing.items():
			if v:

				r = requests.get(f'https://www.materialsproject.org/rest/v2/materials/{material}/vasp/dos?API_KEY={API_KEY}')
				r = r.json()['response'][0]['dos']

				#get different properties of the DOS array
				energies =list(r['energies'])
				densities = list(r['densities']['1'])
				efermi = r['efermi']

				missingMaterials.append([material, energies, densities, efermi])


		df = pd.DataFrame(missingMaterials, columns = ['ID', 'Energies', 'Densities', 'efermi'])
		df.to_csv('dos.csv', mode='a', header=False)


if __name__ == '__main__':
	plot_periodic_table()
