#!/usr/bin/env python3
"""
Methods for calculating Jacard and L3 scores.
"""

import csv
import math


def parse_data(file_location):
    """
    Parses data into a format for the calculation.
    """
    is_on_header = True
    interactions = {}
    raw_data = {}
    with open(file_location, newline='') as csvfile:
        file_content = csv.reader(csvfile, delimiter=' ', quotechar='|')

        for row in file_content:
            if not is_on_header:
                file_data = row[0].split(',')
                data = interactions.get(file_data[0], [])
                if float(file_data[2]) > 0.1:
                    data.append(file_data[1])
                    interactions[file_data[0]] = data
                raw_data[(file_data[0], file_data[1])] = file_data[2:]
            else:
                is_on_header = False

    return interactions, raw_data


def l3_build_paths(nodes, source):
    """
    Builds a list utilizing Length 3 for exploration.
    """
    paths = set()
    for layer1 in nodes[source]:
        if layer1 not in nodes:
            continue
        for layer2 in nodes[layer1]:
            if layer2 not in nodes:
                continue
            for layer3 in nodes[layer2]:
                if layer3 not in nodes:
                    continue
                path = (source, layer1, layer2, layer3)
                if path not in paths:
                    paths.add(path)

    return list(paths)


def calculate_jaccard(list1, list2):
    """
    Calculates the Jaccard score between two lists
    of PPIs.
    """
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union


def calculate_l3_score(nodes, path):
    """
    Calculates the L3 score between two lists
    of PPIs.
    """
    second_node_degree = len(nodes[path[1]])
    third_node_degree = len(nodes[path[2]])
    score = 1 / math.sqrt(second_node_degree*third_node_degree)
    return score


def build_jaccard_and_l3_values(interactions):
    """
    Builds jaccard and l3 score values
    """
    tcp_jaccard = {}
    l3_scores = {}
    nodes = interactions.items()

    for protein_a in nodes:
        for protein_b in nodes:
            if protein_a == protein_b:
                continue
            score = calculate_jaccard(protein_a[1], protein_b[1])
            tcp_jaccard[(protein_a[0], protein_b[0])] = score
            l3_scores[(protein_a[0], protein_b[0])] = 0

        protein_a_l3_paths = l3_build_paths(interactions, protein_a[0])
        for path in protein_a_l3_paths:
            curr_score = l3_scores.get((path[0], path[3]), 0)
            l3_scores[(path[0], path[3])] = curr_score + \
                calculate_l3_score(interactions, path)
    return tcp_jaccard, l3_scores


def format_output_data(raw_data, tcp_jaccard, l3_scores):
    """
    Formats raw_data, jaccard and l3 values for writing to file
    """
    output_data = {}
    for proteins, l3_score in l3_scores.items():
        data = raw_data.get(proteins, raw_data.get(
            (proteins[1], proteins[0]), None))
        if data is None:
            continue
        entry = {'protein1': proteins[0],
                 'protein2': proteins[1],
                 'interaction_score': data[0],
                 'similarity': data[1],
                 'interaction_label': data[2],
                 'jaccard': tcp_jaccard[proteins],
                 'l3_score': l3_score,
                 'sequence_1': data[3],
                 'sequence_2': data[4]}
        if proteins in output_data or (proteins[1], proteins[0]) in output_data:
            continue
        output_data[proteins] = entry
    return output_data


def main():
    """
    Main execution point for the calculation.
    """
    interactions, raw_data = parse_data("./data/interaction_data.csv")
    tcp_jaccard, l3_scores = build_jaccard_and_l3_values(interactions)
    output_data = format_output_data(raw_data, tcp_jaccard, l3_scores)

    # print(entry)
    with open("./data/interaction_data2.csv", mode='w', newline='') as csv_file:
        fieldnames = ['protein1', 'protein2', 'interaction_score',
                      'similarity', 'interaction_label',
                      'jaccard', 'l3_score', 'sequence_1', 'sequence_2']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()
        for data in output_data.values():
            writer.writerow(data)


if __name__ == '__main__':
    main()
