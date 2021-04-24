#!/usr/bin/env python3

"""
This is the utils class for the project
"""

import csv

__author__ = "Alex DeCesare"
__version__ = "10-March-2021"


def parse_values(file_location, columns):
    """
    Parses a file into a list of data with a given number of columns.
    """
    score_and_similarity = []
    is_on_header = True

    with open(file_location, newline='') as csvfile:
        file_content = csv.reader(csvfile, delimiter=' ', quotechar='|')

        for row in file_content:
            if not is_on_header:
                interaction_file_data = row[0].split(',')

                score_and_similarity.append(
                    [float(interaction_file_data[i]) for i in columns])

            else:
                is_on_header = False

    return score_and_similarity


def parse_is_interacted(file_location):
    """
    Parses a file into a list of proteins that do interact.
    """
    is_interacted = []
    is_on_header = True

    with open(file_location, newline='') as csvfile:
        file_content = csv.reader(csvfile, delimiter=' ', quotechar='|')

        for row in file_content:
            if not is_on_header:
                interaction_file_data = row[0].split(',')
                is_interacted.append(int(interaction_file_data[4]))
            else:
                is_on_header = False

    return is_interacted
