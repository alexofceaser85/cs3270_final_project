#!/usr/bin/env python3
'''
Parses out information from UNIPROT into various formats.clscls
'''
import urllib.parse
import urllib.request
import csv
import random
from progress.bar import IncrementalBar
import uniprot_parser
import sequence_similarity

__author__ = "Alexander Ayers"
__version__ = "2-17-2021"

DEGREE = 33
SIGMA = 1.0
KAPPA = 0.37


def random_combination(iterable, choices):
    """
    Makes a random combination of an iterator.
    """
    pool = tuple(iterable)
    length = len(pool)
    indices = sorted(random.sample(range(length), choices))
    return tuple(pool[i] for i in indices)


def read_interaction_file(loaded_filename):
    """
    Read the interaction file into protein and interaction collections
    """
    full_proteins = []
    interactions_string = {}
    with open(loaded_filename, newline="") as interactions_file:
        interaction_reader = csv.DictReader(interactions_file)

        for interaction in interaction_reader:
            protein1 = interaction["node1_string_id"]
            protein2 = interaction["node2_string_id"]
            interactions_string[protein1,
                                protein2] = interaction["combined_score"]
            if protein1 not in full_proteins:
                full_proteins.append(protein1)
            if protein2 not in full_proteins:
                full_proteins.append(protein2)
    return full_proteins, interactions_string


def send_and_recieve_interactions(full_proteins):
    """
    Send request to unpiprot and read response data
    """
    url = 'https://www.uniprot.org/uploadlists/'
    protein_line = " ".join(full_proteins)
    params = {
        'from': 'STRING_ID',
        'to': 'ACC',
        'format': 'txt',
        'query': protein_line
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as url_response:
        response = url_response.read()

    protein_data, id_conversion = uniprot_parser.parse_txt_file(
        response.decode("utf-8"))
    return protein_data, id_conversion


def calculate_similarity_and_format_output(total_interactions, interactions_string,
                                           id_conversion, protein_data):
    """
    Calculate similarities and format for output
    """
    progress_bar = IncrementalBar('Processing', max=total_interactions)
    output_data = []
    for proteins in interactions_string.items():
        if proteins[0][0] in id_conversion and proteins[0][1] in id_conversion:
            uniprot1 = id_conversion[proteins[0][0]]
            uniprot2 = id_conversion[proteins[0][1]]

            uniprot1_sequeunce = protein_data[uniprot1]["sequence"]
            uniprot2_sequeunce = protein_data[uniprot1]["sequence"]

            similarity = sequence_similarity.compute_similarity(
                uniprot1_sequeunce, uniprot2_sequeunce, DEGREE, SIGMA, KAPPA)

            interaction_score = float(
                interactions_string[proteins[0][0], proteins[0][1]])
            interaction_label = 1 if interaction_score > 0.5 else 0

            output_data.append(
                {'protein1': uniprot1,
                 'protein2': uniprot2,
                 'interaction_score': interaction_score,
                 'similarity': similarity,
                 'interaction_label': interaction_label,
                 'sequence_1': uniprot1_sequeunce,
                 'sequence_2': uniprot2_sequeunce})
            progress_bar.next()

    progress_bar.finish()
    return output_data


def convert_id_and_calculate_similarity(loaded_filename, saved_filename):
    """
    Converts a Uniprot report into a CSV of Sequences and Other Data.
    """

    full_proteins, interactions_string = read_interaction_file(loaded_filename)
    protein_data, id_conversion = send_and_recieve_interactions(full_proteins)

    all_interactions = []
    total_interactions = 0
    for proteins in interactions_string.items():
        if proteins[0][0] in id_conversion and proteins[0][1] in id_conversion:
            total_interactions = total_interactions + 1
            all_interactions.append(proteins[0][0])
            all_interactions.append(proteins[0][1])

    for _ in range(1000):
        combination = random_combination(all_interactions, 2)
        if combination not in interactions_string:
            interactions_string[combination] = 0
            total_interactions = total_interactions + 1

    output_data = calculate_similarity_and_format_output(
        total_interactions, interactions_string, id_conversion, protein_data)

    with open(saved_filename, mode='w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=['protein1',
                                                      'protein2',
                                                      'interaction_score',
                                                      'similarity',
                                                      'interaction_label',
                                                      'sequence_1',
                                                      'sequence_2'])

        writer.writeheader()
        for data in output_data:
            writer.writerow(data)


if __name__ == "__main__":
    convert_id_and_calculate_similarity(
        "./data/interactions.csv", "./data/interaction_data.csv")
