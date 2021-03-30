#!/usr/bin/env python3

"""
parses the content of the given file into clusters based upon a given range
"""

from numpy import concatenate
from numpy import array
from sgt import SGT
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pprint
import utils
import random
import csv


def main():

	print("Enter the number of clusters")
	degree = int(input())

	# define dataset
	print('____________________________________')
	print()
	print('Experiment 3')

	transformed_pca = array(utils.parse_values("./protein_embedding.csv", [ 1, 2]))

	data_frame = pd.DataFrame(data=transformed_pca, columns=['x1', 'x2'])


	data_frame['x1'] = data_frame['x1'].map(float)
	data_frame['x2'] = data_frame['x2'].map(float)

	X = array(data_frame)
	y = array(utils.parse_is_interacted("./interaction_data.csv"))
	print(X)
	kfold = KFold(n_splits=5, random_state=1, shuffle=True)
	scores = []
	for train, test in kfold.split(X):
		X_train, X_test = X[train], X[test]
		y_train, y_test = y[train], y[test]

		norm = MinMaxScaler().fit(X_train)
		X_train_norm = norm.transform(X_train)
		X_test_norm = norm.transform(X_test)

		model = LabelPropagation(max_iter = 1000, n_jobs = -1)
		# fit model on training dataset
		model.fit(X_train_norm, y_train)
		# make predictions on hold out test set
		yhat = model.predict(X_test_norm)
		# calculate score for test set
		score = accuracy_score(y_test, yhat)
		scores.append(score*100)
		# summarize score
		print('Accuracy: %.3f' % (score*100))
	print('Mean Accuracy: %.3f (SD: %.3f)' % (np.mean(scores), np.std(scores)))
	print('____________________________________')


	kmeans = KMeans(n_clusters=degree, max_iter=1000)
	kmeans.fit(data_frame)

	labels = kmeans.predict(data_frame)


	data_frame['label'] = labels

	plt.figure(figsize=(5, 5))
	colmap = generate_colors(degree)
	colors = list(map(lambda x: colmap[x+1], labels))
	plt.scatter(data_frame['x1'], data_frame['x2'], color=colors, alpha=0.5)
	plt.show()

def generate_colors(number_of_colors):
	"""
	generates the differant colors for the scatter plot
	"""
	colors = list(mcolors.CSS4_COLORS)
	random.shuffle(colors)

	colors_to_use = {}
	index = 0

	while index < number_of_colors:
		colors_to_use[index+1] = colors[index]
		index = index + 1
	return colors_to_use


if __name__ == '__main__':
	main()
