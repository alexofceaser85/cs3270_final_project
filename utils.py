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
        file_content = csv.reader(csvfile, delimiter=' ', quotechar='|')

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
        file_content = csv.reader(csvfile, delimiter=' ', quotechar='|')

        for row in file_content:
            if (is_on_header == False):
                interaction_file_data = row[0].split(',')
                is_interacted.append(int(interaction_file_data[4]))
            else:
                is_on_header = False

    return is_interacted
