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
import chemparse



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
		response = MPRester(API_KEY).query(criteria={'task_id': {'$in': list(materials)}}, properties=['pretty_formula', 'task_id'])

		for i, material in enumerate(response):
			if not i % 100: 
				print(i)


			newEntry = {column: 0 for column in columns}

			newEntry['ID'] = material['task_id']
			newEntry['Formula'] = material['pretty_formula']

			elements = parse_formula(material['pretty_formula'])
			for element, num in elements.items():
				newEntry[element] += num

			newEntry['Is_Candidate'] = int(material['task_id'] in candidates)

			data.append(newEntry)

	df = pd.DataFrame(data)	
	df.to_csv('clustering/elemental_analysis.csv', mode = 'a', header=False, index_label = False)



if __name__ == '__main__':
	elements_in_candidates()