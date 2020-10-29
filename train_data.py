import tensorflow as tf
import pandas as pd
import itertools
from periodictable import *
import numpy as np
import chemparse
import tensorflow as tf 
import random

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



#parse the pretty formulas fom training_data.csv and tokenize them such that they can be read by a neural network
def preprocess_data():

	#open the CSV file into chunks
	reader = pd.read_csv('training_data.csv')

	
	formulas = list(reader['Formula'].to_numpy())
	labels = list(reader['Is_Candidate'].to_numpy().astype(int))
	data = list(zip(formulas, labels))

	random.shuffle(data)

	formulas, labels = zip(*data)
	formulas = list(formulas)
	labels = list(labels)



	parsed_formulas = []
	for formula in formulas:
		parsed_formulas.append(parse_formula(str(formula)))

	training_size = 30000

	training_formulas = parsed_formulas[:training_size]
	testing_formulas = parsed_formulas[:training_size]

	training_labels = labels[:training_size]
	testing_labels = labels[:training_size]


	token = tf.keras.preprocessing.text.Tokenizer(
    num_words=None, lower=False,
    oov_token=None
)

	token.fit_on_texts(training_formulas)

	training_sequences = token.texts_to_sequences(training_formulas)
	training_sequences = tf.keras.preprocessing.sequence.pad_sequences(training_sequences)

	testing_sequences = token.texts_to_sequences(testing_formulas)
	testing_sequences = tf.keras.preprocessing.sequence.pad_sequences(testing_sequences)


	model = tf.keras.Sequential([
    tf.keras.layers.Embedding(len(token.word_index), 16, input_length = len(testing_sequences[0])),
    tf.keras.layers.GlobalAveragePooling1D(),
    tf.keras.layers.Dense(24, activation='relu'),
    tf.keras.layers.Dense(1, activation='sigmoid')
])

	model.compile(loss='binary_crossentropy',optimizer='adam',metrics=['accuracy'])

	model.summary()
	
	num_epochs = 30
	history = model.fit(training_sequences, training_labels, epochs=num_epochs, validation_data=(testing_sequences, testing_labels), verbose=2)

if __name__ == '__main__':
	preprocess_data()