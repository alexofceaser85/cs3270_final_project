#!/usr/bin/env python3

"""
parses the content of the given file into clusters based upon a given range
"""
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
from sgt import SGT
import csv

__author__ = "Alex DeCesare"
__version__ = "20-March-2021"

alphabets = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
             'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
             'W', 'X', 'Y', 'U', 'O']

def sequence_embedding(file_path):
    """
    clusters the data into different groups
    """

    corpus = pd.read_csv(file_path)
    corpus = corpus.loc[:, ['protein1', 'sequence_1']]
    
    corpus.columns = ['id', 'sequence']
    corpus['sequence'] = corpus['sequence'].map(list)
    sgt_ = SGT(kappa=1, lengthsensitive=False, mode='multiprocessing', alphabets=alphabets)
    sgtembedding_df = sgt_.fit_transform(corpus)
    sgtembedding_df = sgtembedding_df.set_index('id')

    pca = PCA(n_components=2)
    pca.fit(sgtembedding_df)
    transformed_pca = pca.transform(sgtembedding_df)
    print(np.sum(pca.explained_variance_ratio_))
    with open("protein_embedding.csv", mode='w', newline='') as csv_file:
        fieldnames = ['protein_id', 'x1', 'x2']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        index = 0
        for protein_id, vector in sgtembedding_df.iterrows():
            data = transformed_pca[index]
            writer.writerow(
                {'protein_id': protein_id, 'x1': data[0], 'x2': data[1]})
            index += 1


if __name__ == '__main__':
    sequence_embedding('./interaction_data.csv')
