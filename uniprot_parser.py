#!/usr/bin/env python3

'''
Parses out information from UNIPROT into various formats.
'''

__author__ = "Alexander Ayers"
__version__ = "2-17-2021"


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
  uniprot_id = None
  metadata_by_seqid = {}
  protein_conversion = {}
  index = 0
  for l in cache_txt.splitlines():
    test_tag = l[:5].strip()
    if test_tag and test_tag != tag:
      tag = test_tag
    line = l[5:].strip()
    words = line.split()
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
    if tag == "SQ":
      if words[0] != "SEQUENCE":
        entry['sequence'] += ''.join(words)
    if tag == "AC":
      accs = [w.replace(";", "") for w in words]
      entry['accs'].extend(accs)
      for acc in accs:
        metadata_by_seqid[acc] = entry
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
    if tag == "GN":
      if 'gene' not in entry and len(words) > 0:
        pieces = words[0].split("=")
        if len(pieces) > 1 and 'name' in pieces[0].lower():
          entry['gene'] = pieces[1].replace(';', '').replace(',', '')
    if tag == "OS":
      if 'organism' not in entry:
        entry['organism'] = ""
      entry['organism'] += line
    if tag == "DE":
      if 'descriptions' not in entry:
        entry['descriptions'] = []
      entry['descriptions'].append(line)
    if tag == "CC":
      if 'comment' not in entry:
        entry['comment'] = ''
      entry['comment'] += line + '\n'
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

  return metadata_by_seqid, protein_conversion
