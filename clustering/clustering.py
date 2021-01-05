import periodictable 
from pymatgen import *
import pandas as pd
import numpy as np
from ast import literal_eval
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from bokeh.sampledata import periodic_table
import sklearn
import matplotlib.pyplot as plt # for data visualization
import seaborn as sns # for statistical data visualization



API_KEY = 'UgRqoHkuZyJEVX2d'

def parse_formula(formula):
	array = []
	currentSpecies = ''
	delimiters = 0
	for i, char in enumerate(formula):
		if char == '(':
			delimiters += 1
			if delimiters == 1:
				if currentSpecies:
					array.append(currentSpecies)
				currentSpecies = char
			else:
				currentSpecies += char
		elif char.isupper() and not '(' in currentSpecies:
			if len(currentSpecies) > 0:
				array.append(currentSpecies) 
			currentSpecies = char
		elif char == ')':
			delimiters -= 1
			if delimiters == 0:
				array.append(currentSpecies[1:])
				currentSpecies = ''
			else:
				currentSpecies += char
		elif char.isdigit() and not '(' in currentSpecies:
			if currentSpecies and not currentSpecies[-1].isdigit():
				array.append(currentSpecies)
				currentSpecies = ''
			currentSpecies += char
		else:
			currentSpecies += char
	if currentSpecies:
		array.append(currentSpecies)
	return array
def get_ratio(formula):
	parsed_formula = parse_formula(formula)
	ratio =[]
	for i, elem in enumerate(parsed_formula):
		if not elem.isdigit() and not parsed_formula[i-1].isdigit():
			ratio.append(1)
		elif elem.isdigit():
			ratio.append(int(elem))

	ratio = sorted(ratio)
	ratio = ':'.join([str(num) for num in ratio])
	return ratio

def cosine_similarity(x, y):

	ratio_similarity = x['Ratios'] == y['Ratios']
	
	x = np.concatenate(([x['Heaviest Element']/118.0], [int(ratio_similarity)], x[5:].values))

	#we give matching elemental ratios a weight of .3
	y = np.concatenate(([y['Heaviest Element']/118.0], [.01], y[5:].values))

	x = np.array(x)
	y = np.array(y)

	x.reshape(1, -1)
	y.reshape(1, -1)
	
	return 1 - np.dot(x, y) / (np.sqrt(x.dot(x) * y.dot(y)))

def curate_data(): 
	columns = ['ID', 'Formula', "Heaviest Element", 'Ratios', 'Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5', 'Group 6', 'Group 7', 'Group 8', 'Transition Metals', 'Lanthanides', 'Actinides']
	table_data = periodic_table.elements.copy()

	#setup headers
	df = pd.DataFrame(columns = columns)
	df.to_csv('candidates.csv')

	candidates = [m.rstrip() for m in open('candidates.txt', 'r').readlines()]

	formulas = MPRester(API_KEY).query(criteria={'task_id': {'$in': candidates}}, properties=['pretty_formula', 'task_id', 'elements', 'nelements'])
	data = []

	for i, entry in enumerate(formulas):
		elements = entry['elements']
		#create new entry
		newEntry = {column: 0 for column in columns}
		newEntry['ID'] = entry['task_id']
		newEntry['Formula'] = entry['pretty_formula']

		#define atomic ratios
		newEntry['Ratios'] = get_ratio(entry['pretty_formula'])

		#define heaviest element
		max_elem = max([table_data[table_data['symbol']==element]['atomic number'].values[0] for element in elements])
		newEntry['Heaviest Element'] = max_elem / 118.0

		#Append quantity of elements in each group
		for elem in elements:

			group = table_data[table_data['symbol']==elem]['group'].values[0]
			period = int(table_data[table_data['symbol']==elem]['period'].values[0])
			
			if group == '-':
				if period == 6:
					newEntry['Lanthanides'] += 1
				elif period == 7:
					newEntry['Actinides'] += 1
				continue

			group = int(group)
			
			if group == 1:
				newEntry['Group 1'] += 1
			elif group == 2:
				newEntry['Group 2'] += 1
			elif group == 13:
				newEntry['Group 3'] += 1
			elif group == 14:
				newEntry['Group 4'] += 1
			elif group == 15:
				newEntry['Group 5'] += 1
			elif group == 16:
				newEntry['Group 6'] += 1
			elif group == 17:
				newEntry['Group 7'] += 1
			elif group == 18:
				newEntry['Group 8'] += 1
			elif group > 2 and group < 13:
				newEntry['Transition Metals'] += 1

		data.append(newEntry)

		if not i%100:
			print(f'{i} Computed entries')
	
	df = pd.DataFrame(data)
	df.to_csv('candidates.csv', mode='a', header=False)

def compute_similarities(): 
	df = pd.read_csv('candidates.csv')

	similarity_matrix = []
	for i, row in df.iterrows(): 
		elem_similarity = []
		for j, row2 in df.iterrows():
			elem_similarity.append(cosine_similarity(row, row2))
		similarity_matrix.append(elem_similarity)
		if not i % 10:
			print(f'{i} entries computed')

	df = pd.DataFrame(similarity_matrix)
	df.to_csv('similarity.csv')


def find_elbow_point():
	df = pd.read_csv('candidates.csv', index_col = 0)
	df.drop(['ID', 'Formula', 'Ratios', 'Cluster'], axis=1, inplace=True)
	X = df
	X['Heaviest Element'] /= 118.0
	
	cs = []
	for i in range(1, 20):
	    kmeans = KMeans(n_clusters = i, init = 'k-means++', max_iter = 300, n_init = 10, random_state = 0)
	    kmeans.fit(X)
	    cs.append(kmeans.inertia_)
	plt.plot(range(1, 20), cs)
	plt.title('The Elbow Method')
	plt.xlabel('Number of clusters')
	plt.ylabel('CS')
	plt.show()

def kmeans_cluster():

	df = pd.read_csv('candidates.csv', index_col = 0)
	df.drop(['ID', 'Formula', 'Ratios', 'Cluster'], axis=1, inplace=True)

	
	df['Heaviest Element'] /= 118

	X = df.values	
	
	kmeans = KMeans(n_clusters=3,random_state=0)

	kmeans.fit(X)

	columns = ["Heaviest Element", 'Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5', 'Group 6', 'Group 7', 'Group 8', 'Transition Metals', 'Lanthanides', 'Actinides']
	df = pd.DataFrame(kmeans.cluster_centers_, columns = columns)
	df.to_csv('kmeans_clusters.csv', index_label = False)

def cosine_cluster():
	df = pd.read_csv('similarity.csv', index_col = 0)
	mat = df.values

	spectral = SpectralClustering(4, affinity = 'precomputed')
	cluster = spectral.fit_predict(mat)
	
	#add cluster labels to csv file
	df = pd.read_csv('candidates.csv', index_col = 0)
	df['Cluster'] = cluster
	df.to_csv('candidates.csv', index_label = False)



if __name__ == '__main__':
	
	cosine_cluster()