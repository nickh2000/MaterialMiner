import pandas as pd
import itertools
import periodictable
import numpy as np
import chemparse
import tensorflow as tf
import random



'''This program is for training a machine learning model using the mined materials data '''




#converts a chemical formula to an array
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

def parse_formula2(formula):
	delimiters = 0
	for char in formula:
		if char == '(':
			delimiters += 1
	if delimiters > 2:
		print(formula)

#parse the pretty formulas fom training_data.csv and tokenize them such that they can be read by a neural network
def preprocess_data():

	#open the CSV file into chunks
	reader = pd.read_csv('training_data.csv')

	#shuffle data and seperate labels and candidates
	formulas = list(reader['Formula'].to_numpy())
	labels = list(reader['Is_Candidate'].to_numpy().astype(int))
	data = list(zip(formulas, labels))

	random.shuffle(data)

	formulas, labels = zip(*data)
	formulas = list(formulas)
	labels = list(labels)

	parsed_formulas = []
	elements = [str(element.symbol) for element in periodictable.elements]
	elements = elements[1:]

	for formula in formulas:
		parsed_formula = []
		if isinstance(formula, str):
			parsed_formula = {k:0 for k in elements}
			for k, v in chemparse.parse_formula(formula).items():
				parsed_formula[k] = v
			parsed_formulas.append(list(parsed_formula.values()))

	print(parsed_formulas[0])


	training_size = 50000

	training_formulas = parsed_formulas[:training_size]
	testing_formulas = parsed_formulas[:training_size]

	training_labels = labels[:training_size]
	testing_labels = labels[:training_size]
	
	training_formulas = np.array(training_formulas)
	training_labels = np.array(training_labels)
	testing_formulas = np.array(testing_formulas)
	testing_labels = np.array(testing_labels)

	print(len(training_formulas[0]))
	model = tf.keras.Sequential([
    tf.keras.layers.InputLayer(input_shape=(len(training_formulas[0]))),
    tf.keras.layers.Dense(24, activation='relu'),
    tf.keras.layers.Dense(1, activation='sigmoid')
])

	model.compile(loss='binary_crossentropy',optimizer='adam',metrics=['accuracy'])

	model.summary()
	
	num_epochs = 30

	history = model.fit(training_formulas, training_labels, epochs=num_epochs, validation_data=(testing_formulas, testing_labels), verbose=2)

if __name__ == '__main__':
	preprocess_data()