#!/usr/bin/env python3

"""
Experiment 3 that utilizes the Euclidean
sequence embedding and Label Propagation
to compare the PPIs.
"""

import math
from progress.bar import IncrementalBar
import numpy as np
import pandas as pd
import base_experiment

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

def main():
    """
    Main execution point for the experiment.
    """
    # define dataset
    distances, interaction_labels = parse_sequence_embedding()
    x_values = np.array(list(distances.values()))
    y_values = np.array(list(interaction_labels.values()))

    base_experiment.execute(x_values, y_values, "3")



if __name__ == '__main__':
    main()
