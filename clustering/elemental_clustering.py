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
import os




API_KEY = 'UgRqoHkuZyJEVX2d'




def elements_in_candidates():
	os.chdir('..')
	elements = [str(element) for element in periodictable.elements][1:]
	columns = ['ID'] + ['Formula'] + elements + ['Is_Candidate']
	df = pd.DataFrame(columns = columns)
	df.to_csv('clustering/elemental_analysis.csv', index_label = False)

	data = []
	candidates = [candidate.rstrip() for candidate in open('candidates.txt').readlines()]

	with open("database.txt") as f:

		materials = [material.rstrip() for material in f.readlines()]
		response = MPRester(API_KEY).query(criteria={'task_id': {'$in': list(materials)}}, properties=['elements', 'pretty_formula', 'task_id'])

		for i, material in enumerate(response):
			if not i % 100: 
				print(i)


			newEntry = {column: 0 for column in columns}

			newEntry['ID'] = material['task_id']
			newEntry['Formula'] = material['pretty_formula']

			elements = material['elements']
			for element in elements:
				newEntry[element] += 1

			newEntry['Is_Candidate'] = int(material['task_id'] in candidates)

			data.append(newEntry)

	df = pd.DataFrame(data)	
	df.to_csv('clustering/elemental_analysis.csv', mode = 'a', header=False, index_label = False)



if __name__ == '__main__':
	elements_in_candidates()