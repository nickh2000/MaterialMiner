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


import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from DosData import * 

API_KEY = 'UgRqoHkuZyJEVX2d'


def candidate_dos():

	df = pd.read_csv('candidates.csv')
	os.chdir('..')
	reader = pd.read_csv('dos.csv', chunksize = 5000)	

	with open('candidates.txt') as f:
		materials = [candidate.rstrip() for candidate in f.readlines()]

	for c, chunk in enumerate(reader):
		for material in materials:
			search = chunk[chunk['ID']==material]
			if len(search) > 0:
				df.loc[material, 'Energies'] = chunk[chunk['ID']==material].iloc[0]['Energies']
				df.loc[material, 'Densities'] = chunk[chunk['ID']==material].iloc[0]['Densities']
				df.loc[material, 'efermi'] = chunk[chunk['ID']==material].iloc[0]['efermi']


		print(f'{5000*c} entries processed')

	df.to_csv('clustering/candidates.csv')



def curve_distance(x, y):

	x = DosData(x)
	y = DosData(y)
	x_curve = x.center_dos(0, -2)
	y_curve = y.center_dos(0, -2)

	similarity = 0
	sample_points = [0, -.5, -1, -1.5, -2]
	#iterate through different energies in the first curve
	for energy in sample_points:
		#find the nearest energy in the second curve
		y_index = find_nearest_index(list(y_curve.keys()), energy)
		x_index = find_nearest_index(list(x_curve.keys()), energy)
		#difference between densities at these approximately equal energies
		difference = list(y_curve.values())[y_index] - list(x_curve.values())[x_index]
		similarity += difference**2
	similarity /= len(sample_points)

	return similarity


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


#find the distance between two candidates
def candidate_distance(x, y):
	weights = [1/118.0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .01, 1, 1]

	#1 if ratios are the same, zero otherwise
	ratio_distance = int(x['Ratios'] != y['Ratios'])
	ratio_weight = 0
	
	curve_weight = .3
	x = np.concatenate(([x['Heaviest Element']], x[5:19].values))

	#do the same with the other candidate, ratio similarity is given a weight of zero
	y = np.concatenate(([y['Heaviest Element']], y[5:19].values))

	
	x = np.multiply(np.array(x), weights)
	y = np.multiply(np.array(y), weights)

	x.reshape(1, -1)
	y.reshape(1, -1)
	
	return 1 - np.dot(x, y) / (np.sqrt(x.dot(x) * y.dot(y))) + ratio_weight*ratio_distance

def curate_data(): 
	os.chdir('..')
	chiral_spacegroups = np.concatenate(([1], np.arange(3, 6), np.arange(16, 25), np.arange(75, 81), np.arange(89, 99), np.arange(143, 147), np.arange(149, 156), np.arange(168, 174), np.arange(177, 183), np.arange(195, 200), np.arange(207, 215)))	#[77, 79, 91, 95, 92, 96, 152, 154, 151, 153, 169, 170, 178, 179, 172, 173, 180, 181, 212, 213]
	inversion_symmetry_spacegroups = [1] + list(range(10, 15)) + list(range(47, 74)) + list(range(83, 88)) + list(range(123, 142)) + [147] + [148] + list(range(162, 167)) + [175] + [176] + list(range(191, 194)) + list(range(200, 206)) + list(range(221,230))
	
	columns = ['ID', 'Formula', "Heaviest Element", 'Ratios', 'Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5', 'Group 6', 'Group 7', 'Group 8', 'Transition Metals', 'Lanthanides', 'Actinides', 'Spacegroup', 'Is_Chiral', 'Has_Inversion_Symmetry', 'Energies', 'Densities', 'efermi']
	table_data = periodic_table.elements.copy()

	#setup headers
	df = pd.DataFrame(columns = columns)
	df.to_csv('clustering/candidates.csv')

	candidates = [m.rstrip() for m in open('candidates.txt', 'r').readlines()]

	formulas = MPRester(API_KEY).query(criteria={'task_id': {'$in': candidates}}, properties=['pretty_formula', 'task_id', 'spacegroup'])
	data = []

	for i, entry in enumerate(formulas):
		elements = parse_formula(entry['pretty_formula'])
		#create new entry
		newEntry = {column: 0 for column in columns}
		newEntry['ID'] = entry['task_id']
		newEntry['Formula'] = entry['pretty_formula']
		newEntry['Is_Chiral'] = int(entry['spacegroup']['number'] in chiral_spacegroups)
		newEntry['Has_Inversion_Symmetry'] = int(entry['spacegroup']['number'] in inversion_symmetry_spacegroups)
		newEntry['Spacegroup'] = entry['spacegroup']['number']
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
	materials= df['ID'].values
	df.set_index('ID', inplace=True, drop=False)
	reader = pd.read_csv('dos.csv', chunksize = 5000)
	
	for c, chunk in enumerate(reader):
		for r, row in chunk.iterrows():
			material = row['ID']
			if material in materials:
				df.loc[material, 'Energies'] = row['Energies']
				df.loc[material, 'Densities'] = row['Densities']
				df.loc[material, 'efermi'] = row['efermi']
		print(c)

	df.to_csv('clustering/candidates.csv', mode='a', header=False)

def compute_distances():

	df = pd.read_csv('candidates.csv')

	curve_weight = .3
	num_candidates = len(df.values)
	
	similarity_matrix = np.zeros((num_candidates, num_candidates))

	for i, row in df.iterrows(): 
		elem_similarity = []
		
		for j in np.arange(i, num_candidates):
			#Get DOS objects for two candidates

			x = row
			y = df.iloc[j]

			elem_similarity.append(candidate_distance(x, y) + curve_distance(x, y)*curve_weight)
		
		elem_similarity = np.concatenate((np.zeros(num_candidates - len(elem_similarity)), elem_similarity))
		similarity_matrix[i] += elem_similarity
		
		print(f'{i+1} entries computed')

	reflection = similarity_matrix.T 
	similarity_matrix += reflection

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

def spectral_cluster():
	df = pd.read_csv('similarity.csv', index_col = 0)
	mat = df.values

	spectral = SpectralClustering(4, affinity = 'precomputed')
	cluster = spectral.fit_predict(mat)
	
	#add cluster labels to csv file
	df = pd.read_csv('candidates.csv', index_col = 0)
	df['Cluster'] = cluster
	df.to_csv('candidates.csv', index_label = False)



if __name__ == '__main__':
	compute_distances()
	
