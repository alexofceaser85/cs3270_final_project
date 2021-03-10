#!/usr/bin/env python3

import csv

"""
This is the utils class for the project
"""

__author__ = "Alex DeCesare"
__version__ = "10-March-2021"

def parse_csv(file_location):

    """
    parses a csv into a dictionary and returns that dictionary
    """

    aq1_interactions = {}

    with open(file_location, newline='') as csvfile:
        file_content = csv.reader(csvfile, delimiter = ' ', quotechar = '|')
        counter = 0

        for row in file_content:
            interaction_file_data = row[0].split(',')
            aq1_interactions[(interaction_file_data[0], interaction_file_data[1])] = interaction_file_data[2]

    return aq1_interactions
