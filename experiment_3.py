#!/usr/bin/env python3

"""
Experiment 3 that utilizes the Euclidean
sequence embedding and Label Propagation
to compare the PPIs.
"""

import math
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.semi_supervised import LabelPropagation
from sklearn.preprocessing import MinMaxScaler
from progress.bar import IncrementalBar
import numpy as np
import pandas as pd


def calculate_distancs(embedding_frame, interaction_labels_orig):
    """
    Calculates the distances between protein interactions.
    """
    progress_bar = IncrementalBar('Processing', max=len(embedding_frame)**2)
    distances = {}
    interaction_labels = {}
    for _, vector_a in embedding_frame.iterrows():
        for _, vector_b in embedding_frame.iterrows():
            index = (vector_a['protein_id'], vector_b['protein_id'])
            if index in interaction_labels_orig:
                data = interaction_labels_orig.get(index)
                distance = math.sqrt(
                    (vector_a['x1'] - vector_b['x1']) ** 2 + (vector_b['x2'] - vector_b['x2']) ** 2)
                interaction_labels[index] = data['interaction_label']
                distances[index] = [distance]
            progress_bar.next()
    progress_bar.finish()
    return distances, interaction_labels


def parse_sequence_embedding():
    """
    Parses the embedding given a CSV of sequences.
    """
    embedding_frame = pd.read_csv('./data/protein_embedding.csv')
    embedding_frame.columns = ['protein_id', 'x1', 'x2']
    embedding_frame['x1'] = embedding_frame['x1'].map(float)
    embedding_frame['x2'] = embedding_frame['x2'].map(float)
    embedding_frame = embedding_frame.drop_duplicates(subset=['protein_id'])

    x_frame = embedding_frame.loc[:, ['x1', 'x2']]
    x_frame.columns = ['x1', 'x2']

    interaction_frame = pd.read_csv('./data/interaction_data2.csv')
    interaction_frame.columns = ['protein1', 'protein2', 'interaction_score',
                                 'similarity', 'interaction_label', 'jaccard',
                                 'l3_score', 'sequence_1', 'sequence_2']
    interaction_frame['interaction_label'] = interaction_frame['interaction_label'].map(
        int)
    interaction_frame.set_index(['protein1', 'protein2'])

    interaction_labels_orig = {}
    for _, vector in interaction_frame.iterrows():
        interaction_labels_orig[vector['protein1'],
                                vector['protein2']] = vector

    return calculate_distancs(embedding_frame, interaction_labels_orig)

def calculate_scores(x_value, y_value):
    """
    calculates the accuracy scores
    """
    kfold = KFold(n_splits=5, random_state=1, shuffle=True)
    scores = []
    for train, test in kfold.split(x_value):
        x_train, x_test = x_value[train], x_value[test]
        y_train, y_test = y_value[train], y_value[test]

        norm = MinMaxScaler().fit(x_train)
        x_train_norm = norm.transform(x_train)
        x_test_norm = norm.transform(x_test)

        yhat = create_label_propagation_model(x_train_norm, y_train, x_test_norm)
        # calculate score for test set
        score = accuracy_score(y_test, yhat)
        scores.append(score*100)
        # summarize score
    return scores

def create_label_propagation_model(x_train_norm, y_train, x_test_norm):
    """
    creates a new label propagation model
    """
    model = LabelPropagation(max_iter=1000, n_jobs=-1)
    # fit model on training dataset
    model.fit(x_train_norm, y_train)
    # make predictions on hold out test set
    return model.predict(x_test_norm)

def main():
    """
    Main execution point for the experiment.
    """
    # define dataset
    print('____________________________________')
    print()
    print('Experiment 3')
    distances, interaction_labels = parse_sequence_embedding()

    x_value = np.array(list(distances.values()))
    y_value = np.array(list(interaction_labels.values()))

    scores = calculate_scores(x_value, y_value)
    for score in scores:
        print('Accuracy: %.3f' % (score))

    print('Mean Accuracy: %.3f (SD: %.3f)' % (np.mean(scores), np.std(scores)))
    print('____________________________________')


if __name__ == '__main__':
    main()
