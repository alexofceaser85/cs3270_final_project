#!/usr/bin/env python3

"""
parses the content of the given file into clusters based upon a given range
"""
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import pandas as pd
from sgt import SGT
import matplotlib.pyplot as plt

__author__ = "Alex DeCesare"
__version__ = "20-March-2021"

DEGREE = 33
SIGMA = 1.0
KAPPA = 10

def split_clusters(file_path):

    """
    clusters the data into different groups
    """

    print("Enter the number of clusters")
    degree = int(input())

    corpus = pd.read_csv(file_path)
    corpus = corpus.loc[:, ['protein1', 'similarity']]
    corpus.columns = ['sequence', 'id']
    corpus['sequence'] = corpus['sequence'].map(list)
    sgt_ = SGT(kappa=1, lengthsensitive=False, mode='multiprocessing')
    sgtembedding_df = sgt_.fit_transform(corpus)
    sgtembedding_df.set_index('id')

    pca = PCA(n_components=2)
    pca.fit(sgtembedding_df)
    transformed_pca = pca.transform(sgtembedding_df)

    print(np.sum(pca.explained_variance_ratio_))

    data_frame = pd.DataFrame(data=transformed_pca, columns=['x1', 'x2'])
    data_frame.head()

    kmeans = KMeans(n_clusters=degree, max_iter=300)
    kmeans.fit(data_frame)

    labels = kmeans.predict(data_frame)

    plt.figure(figsize=(5, 5))
    colmap = generate_colors(degree)
    print(colmap)
    colors = list(map(lambda x: colmap[x+1], labels))
    plt.scatter(data_frame['x1'], data_frame['x2'], color=colors, alpha=0.5)
    plt.show()
    print(corpus)

def generate_colors(number_of_colors):
    """
    generates the differant colors for the scatter plot
    """
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    colors_to_use = {}
    index = 0

    while index < number_of_colors:
        colors_to_use[index+1] = colors[index]
        index = index + 1
    return colors_to_use

def split(word):
    """
    splits a word into letters
    """
    return [char for char in word]

if __name__ == '__main__':
    split_clusters('./interaction_data.csv')