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
import chemparse
import parser
import re
from math import gcd
from functools import reduce

API_KEY = 'UgRqoHkuZyJEVX2d'



def parse_formula(formula):
	elements  = chemparse.parse_formula(formula)
	for element in list(elements.keys()):
		if len(element) > 2:
			correct_elements = parse_segment(element)
			elements.pop(element)
			for elem in correct_elements:
				elements[elem] = correct_elements[elem]

	return elements

def parse_segment(segment):


	elems = {}
	currentElem = ''
	for i, char in enumerate(segment):
		if char.isdigit():
			elems[currentElem] = int(char)
			currentElem = ''
		elif char.isupper():
			if i != 0:
				elems[currentElem] = 1
			currentElem = char
		else: 
			currentElem += char
	if len(currentElem) > 0:
		elems[currentElem] = 1
	return elems


def get_ratio(parsed_formula):
	nums = np.array(list(parsed_formula.values())).astype(int)
	factor  = reduce(gcd, nums)
	parsed_formula = {elem: num/factor for elem, num in parsed_formula.items()}

	ratio = ''
	for num in sorted(parsed_formula.values()):
		ratio += str(int(num)) + ':'		
	ratio = ratio[:-1]

	return ratio

def cosine_similarity(x, y):

	#1 if ratios are the same, zero otherwise
	ratio_similarity = x['Ratios'] == y['Ratios']
	
	#get the normalized heaviest element, similarity, and group values
	x = np.concatenate(([x['Heaviest Element']/118.0], [int(ratio_similarity)], x[5:].values))

	#do the same with the other candidate, ratio similarity is given a weight of zero
	y = np.concatenate(([y['Heaviest Element']/118.0], [1], y[5:].values))

	x = np.array(x)
	y = np.array(y)

	x.reshape(1, -1)
	y.reshape(1, -1)
	
	return 1 - np.dot(x, y) / (np.sqrt(x.dot(x) * y.dot(y)))

def curate_data(): 
	os.chdir('..')
	columns = ['ID', 'Formula', "Heaviest Element", 'Ratios', 'Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5', 'Group 6', 'Group 7', 'Group 8', 'Transition Metals', 'Lanthanides', 'Actinides']
	table_data = periodic_table.elements.copy()

	#setup headers
	df = pd.DataFrame(columns = columns)
	df.to_csv('clustering/candidates.csv')

	candidates = [m.rstrip() for m in open('candidates.txt', 'r').readlines()]

	formulas = MPRester(API_KEY).query(criteria={'task_id': {'$in': candidates}}, properties=['pretty_formula', 'task_id'])
	data = []

	for i, entry in enumerate(formulas):
		elements = parse_formula(entry['pretty_formula'])
		#create new entry
		newEntry = {column: 0 for column in columns}
		newEntry['ID'] = entry['task_id']
		newEntry['Formula'] = entry['pretty_formula']

		#define heaviest element
		max_elem = max([table_data[table_data['symbol']==element]['atomic number'].values[0] for element in list(elements.keys())])
			
		newEntry['Heaviest Element'] = max_elem

		#Append quantity of elements in each group
		for elem, num in elements.items():
			num = int(num)
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
				newEntry['Group 1'] += num
			elif group == 2:
				newEntry['Group 2'] += num
			elif group == 13:
				newEntry['Group 3'] += num
			elif group == 14:
				newEntry['Group 4'] += num
			elif group == 15:
				newEntry['Group 5'] += num
			elif group == 16:
				newEntry['Group 6'] += num
			elif group == 17:
				newEntry['Group 7'] += num
			elif group == 18:
				newEntry['Group 8'] += num
			elif group > 2 and group < 13:
				newEntry['Transition Metals'] += num


		#define atomic ratios
		newEntry['Ratios'] = get_ratio(elements)

		data.append(newEntry)

		if not i%100:
			print(f'{i} Computed entries')
	
	df = pd.DataFrame(data)
	df.to_csv('clustering/candidates.csv', mode='a', header=False)

def compute_similarities(): 
	os.chdir('clustering')
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
	
	curate_data()
	compute_similarities()
	cosine_cluster()
	