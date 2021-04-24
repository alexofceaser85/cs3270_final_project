#!/usr/bin/env python3

"""
Experiment 3 that utilizes the Euclidean 
sequence embedding and Label Propagation 
to compare the PPIs.
"""

from numpy import array
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
from sklearn.preprocessing import MinMaxScaler
from progress.bar import IncrementalBar
import numpy as np
import pandas as pd
import math


def calculate_distancs(embedding_frame, interaction_labels_orig):
    """
    Calculates the distances between protein interactions.
    """
    bar = IncrementalBar('Processing', max=len(embedding_frame)**2)
    distances = {}
    interaction_labels = {}
    for _, vectorA in embedding_frame.iterrows():
        for _, vectorB in embedding_frame.iterrows():
            index = (vectorA['protein_id'], vectorB['protein_id'])
            if index in interaction_labels_orig:
                data = interaction_labels_orig.get(index)
                distance = math.sqrt(
                    (vectorA['x1'] - vectorB['x1']) ** 2 + (vectorA['x2'] - vectorB['x2']) ** 2)
                interaction_labels[index] = data['interaction_label']
                distances[index] = [distance]
            bar.next()
    bar.finish()
    return distances, interaction_labels


def parse_sequence_embedding():
    """
    Parses the embedding given a CSV of sequences.
    """
    embedding_frame = pd.read_csv('./protein_embedding.csv')
    embedding_frame.columns = ['protein_id', 'x1', 'x2']
    embedding_frame['x1'] = embedding_frame['x1'].map(float)
    embedding_frame['x2'] = embedding_frame['x2'].map(float)
    embedding_frame = embedding_frame.drop_duplicates(subset=['protein_id'])

    X_frame = embedding_frame.loc[:, ['x1', 'x2']]
    X_frame.columns = ['x1', 'x2']

    interaction_frame = pd.read_csv('./interaction_data2.csv')
    interaction_frame.columns = ['protein1', 'protein2', 'interaction_score',
                                 'similarity', 'interaction_label', 'jaccard', 'l3_score', 'sequence_1', 'sequence_2']
    interaction_frame['interaction_label'] = interaction_frame['interaction_label'].map(
        int)
    interaction_frame.set_index(['protein1', 'protein2'])

    interaction_labels_orig = {}
    for _, vector in interaction_frame.iterrows():
        interaction_labels_orig[vector['protein1'],
                                vector['protein2']] = vector

    return calculate_distancs(embedding_frame, interaction_labels_orig)


def main():
    """
    Main execution point for the experiment.
    """
    # define dataset
    print('____________________________________')
    print()
    print('Experiment 3')
    distances, interaction_labels = parse_sequence_embedding()

    X = array(list(distances.values()))
    y = array(list(interaction_labels.values()))
    kfold = KFold(n_splits=5, random_state=1, shuffle=True)
    scores = []
    for train, test in kfold.split(X):
        X_train, X_test = X[train], X[test]
        y_train, y_test = y[train], y[test]

        norm = MinMaxScaler().fit(X_train)
        X_train_norm = norm.transform(X_train)
        X_test_norm = norm.transform(X_test)

        model = LabelPropagation(max_iter=1000, n_jobs=-1)
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


if __name__ == '__main__':
    main()
