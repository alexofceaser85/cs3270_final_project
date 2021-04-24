#!/usr/bin/env python3
'''
Parses out information from UNIPROT into various formats.clscls
'''
import urllib.parse
import urllib.request
import csv
import uniprot_parser
import sequence_similarity
import pprint
import random
from progress.bar import IncrementalBar

__author__ = "Alexander Ayers"
__version__ = "2-17-2021"

degree = 33
sigma = 1.0
kappa = 0.37


def random_combination(iterable, r):
	pool = tuple(iterable)
	n = len(pool)
	indices = sorted(random.sample(range(n), r))
	return tuple(pool[i] for i in indices)

def convert_id_and_calculate_similarity(loaded_filename, saved_filename):
	url = 'https://www.uniprot.org/uploadlists/'

	full_proteins = []
	interactions_string = {}
	with open(loaded_filename, newline="") as interactions_file:
		interaction_reader = csv.DictReader(interactions_file)

		for interaction in interaction_reader:
			protein1 = interaction["node1_string_id"]
			protein2 = interaction["node2_string_id"]
			interactions_string[protein1, protein2] = interaction["combined_score"]
			if protein1 not in full_proteins:
				full_proteins.append(protein1)
			if protein2 not in full_proteins:
				full_proteins.append(protein2)
	
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
	with urllib.request.urlopen(req) as f:
		response = f.read()

	interactions_uniprot = {}
	protein_data, id_conversion = uniprot_parser.parse_txt_file(response.decode("utf-8"))


	all_interactions = []
	total_interactions = 0
	for protein1, protein2 in interactions_string:
		if protein1 in id_conversion and protein2 in id_conversion:
			total_interactions = total_interactions + 1
			all_interactions.append(protein1)
			all_interactions.append(protein2)

	for _ in range(1000):
		combination = random_combination(all_interactions, 2)
		if combination not in interactions_string:
			interactions_string[combination] = 0
			total_interactions = total_interactions + 1

	bar = IncrementalBar('Processing', max=total_interactions)
	output_data = []
	for protein1,protein2 in interactions_string:
		if protein1 in id_conversion and protein2 in id_conversion:
			uniprot1 = id_conversion[protein1]
			uniprot2 = id_conversion[protein2]
			
			uniprot1_data = protein_data[uniprot1]
			uniprot2_data = protein_data[uniprot2]

			uniprot1_sequeunce = uniprot1_data["sequence"]
			uniprot2_sequeunce = uniprot2_data["sequence"]

			similarity = sequence_similarity.compute_similarity(
				uniprot1_sequeunce, uniprot2_sequeunce, degree, sigma, kappa)

			interaction_score = float(interactions_string[protein1, protein2])
			interaction_label = (interaction_score > 0.5 and 1) or 0
			interaction_features = [interaction_score, similarity, interaction_label]

			interactions_uniprot[uniprot1, uniprot2] = interaction_features
			interactions_uniprot[uniprot2, uniprot1] = interaction_features
			output_data.append(
				{'protein1': uniprot1, 'protein2': uniprot2, 'interaction_score': interaction_score, 'similarity': similarity, 'interaction_label': interaction_label, 'sequence_1': uniprot1_sequeunce, 'sequence_2': uniprot2_sequeunce})
			bar.next()

	bar.finish()
	pprint.pprint(interactions_uniprot["P01019", "P29972"])

	with open(saved_filename, mode='w', newline='') as csv_file:
		fieldnames = ['protein1', 'protein2', 'interaction_score', 'similarity', 'interaction_label', 'sequence_1', 'sequence_2']
		writer = csv.DictWriter(csv_file, fieldnames = fieldnames)

		writer.writeheader()
		for data in output_data:
			writer.writerow(data)

if __name__ == "__main__":
    convert_id_and_calculate_similarity("interactions.csv", "interaction_data.csv")
