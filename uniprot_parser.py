#!/usr/bin/env python3

'''
Parses out information from UNIPROT into various formats.
'''

__author__ = "Alexander Ayers"
__version__ = "2-17-2021"


def parse_descriptions(metadata_by_seqid):
    """
    Parses description
    """
    for entry in metadata_by_seqid.values():
        descriptions = entry['descriptions']
        for i in reversed(range(len(descriptions))):
            description = descriptions[i]
            if 'Short' in description or 'Full' in description:
                j = description.find('=')
                descriptions[i] = description[j+1:].replace(';', '')
            else:
                del descriptions[i]
        entry['description'] = '; '.join(descriptions)


def parse_id(tag, words, index, entry):
    """
    Parses id
    """
    if tag == "ID":
        uniprot_id = words[0]
        is_reviewed = words[1].startswith('Reviewed')
        length = int(words[2])
        entry = {
            'index': index,
            'id': uniprot_id,
            'is_reviewed': is_reviewed,
            'length': length,
            'sequence': '',
            'accs': [],
        }
        index = index + 1
    return entry, index


def parse_sequence(tag, entry, words):
    """
    Parses sequence
    """
    if tag == "SQ":
        if words[0] != "SEQUENCE":
            entry['sequence'] += ''.join(words)


def parse_accs(tag, entry, words, metadata_by_seqid, accs):
    """
    Parses ac
    """
    if tag == "AC":
        accs = [w.replace(";", "") for w in words]
        entry['accs'].extend(accs)
        for acc in accs:
            metadata_by_seqid[acc] = entry
    return accs


def parse_dr(tag, entry, words, protein_conversion, accs):
    """
    Parses dr
    """
    if tag == "DR":
        if 'STRING' in words[0]:
            if 'string_id' not in entry:
                entry['string_id'] = []
            ids = [w[:-1] for w in words[1:]]
            ids = filter(lambda w: len(w) > 1, ids)
            entry['string_id'].extend(ids)
            string_id = entry['string_id'][0]
            protein_conversion[accs[0]] = string_id
            protein_conversion[string_id] = accs[0]


def parse_gn( tag, entry, words):
    """
    Parses gn
    """
    if tag == "GN":
        if 'gene' not in entry and len(words) > 0:
            pieces = words[0].split("=")
            if len(pieces) > 1 and 'name' in pieces[0].lower():
                entry['gene'] = pieces[1].replace(';', '').replace(',', '')


def parse_os(line, tag, entry):
    """
    Parses os
    """
    if tag == "OS":
        if 'organism' not in entry:
            entry['organism'] = ""
        entry['organism'] += line


def parse_de(line, tag, entry):
    """
    Parses de
    """
    if tag == "DE":
        if 'descriptions' not in entry:
            entry['descriptions'] = []
        entry['descriptions'].append(line)


def parse_cc(line, tag, entry):
    """
    Parses cc
    """
    if tag == "CC":
        if 'comment' not in entry:
            entry['comment'] = ''
        entry['comment'] += line + '\n'


def parse_txt_file(cache_txt):
    """
    Parses the text of metadata retrieved from uniprot.org.
    Only a few fields have been parsed, but this provides a
    template for the other fields.
    A single description is generated from joining alternative
    descriptions.
    Returns a dictionary with the main UNIPROT ACC as keys.
    """

    tag = None
    metadata_by_seqid = {}
    protein_conversion = {}
    accs = []
    entry = {}
    index = 0
    for line in cache_txt.splitlines():
        test_tag = line[:5].strip()
        if test_tag and test_tag != tag:
            tag = test_tag
        line = line[5:].strip()
        words = line.split()
        entry, index = parse_id(tag, words, index, entry)
        parse_sequence(tag, entry, words)
        accs = parse_accs(tag, entry, words, metadata_by_seqid, accs)
        parse_dr(tag, entry, words, protein_conversion, accs)
        parse_gn(tag, entry, words)
        parse_os(line, tag, entry)
        parse_de(line, tag, entry)
        parse_cc(line, tag, entry)

    parse_descriptions(metadata_by_seqid)

    return metadata_by_seqid, protein_conversion
