#!/usr/bin/env python3

import csv

"""
This is the utils class for the project
"""

__author__ = "Alex DeCesare"
__version__ = "10-March-2021"

def parse_values(file_location, columns):

    interaction_score_and_similarity = []
    is_on_header = True

    with open(file_location, newline='') as csvfile:
            file_content = csv.reader(csvfile, delimiter = ' ', quotechar = '|')

            for row in file_content:
                if (is_on_header == False):
                    interaction_file_data = row[0].split(',')
                   
                    interaction_score_and_similarity.append(
                        [float(interaction_file_data[i]) for i in columns])
                   
                else:
                    is_on_header = False

    return interaction_score_and_similarity


def parse_is_interacted(file_location):

    is_interacted = []
    is_on_header = True

    with open(file_location, newline='') as csvfile:
            file_content = csv.reader(csvfile, delimiter = ' ', quotechar = '|')

            for row in file_content:
                if (is_on_header == False):
                    interaction_file_data = row[0].split(',')
                    is_interacted.append(int(interaction_file_data[4]))
                else:
                    is_on_header = False

    return is_interacted

def parse_csv(file_location):

    """
    Parses a csv into a dictionary and returns that dictionary
    """

    aq1_interactions = {}

    with open(file_location, newline='') as csvfile:
        file_content = csv.reader(csvfile, delimiter = ' ', quotechar = '|')

        for row in file_content:
            interaction_file_data = row[0].split(',')
            aq1_interactions[(interaction_file_data[0], interaction_file_data[1])] = int(interaction_file_data[2])

    return aq1_interactions

def construct_lists_of_features_and_labels(filename):
    """
    Gets all features and labels from given csv file and
    returns them as a tuple of lists.
    """
    interactions = parse_csv(filename)
    proteins = []
    interacts = []

    for key in interactions:
        proteins.append(key)
    for key in interactions:
        interacts.append(interactions[key])

    return (proteins, interacts)
