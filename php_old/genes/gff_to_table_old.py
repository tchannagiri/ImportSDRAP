import pandas as pd
# from Bio.Seq import Seq
# import Bio.SeqIO as SeqIO
# import gffutils
import os

GENE_DIR = 'genes'
OXYTRI_MAC_2020_GENE_GFF = os.path.join(GENE_DIR, 'O_trifallax_2020-upd3.gff3')
OXYTRI_MAC_2020_GENE_CSV = os.path.join(GENE_DIR, 'O_trifallax_2020-upd3.csv')
OXYTRI_MAC_2012_GENE_GFF = os.path.join(GENE_DIR, 'Oxytricha_trifallax_022112.gff3')
OXYTRI_MAC_2012_GENE_CSV = os.path.join(GENE_DIR, 'Oxytricha_trifallax_022112.csv')

# RENAME_ATTR_KEY = {
#   'id': 'attr_id',
#   'name': 'attr_name',
#   'parent': 'attr_parent',
#   'parent_gene': 'attr_parent_gene',
#   'note': 'attr_note',
# }

# def get_db_file(gff_file):
#   return gff_file + '.db'

# def get_csv_file(gff_file):
#   return gff_file + '.csv'

def parse_gff_attrs(attr_str):
  attr_list = attr_str.split(';')

  # Use fuzzy logic to check if ';' ocurred in a field, and combine with the previous field
  attr_list_2 = []
  for attr in attr_list:
    if '=' in attr:
      attr_list_2.append(attr)
    else:
      if len(attr_list_2) == 0:
        raise Exception('Can''t combine the attribute, no previous one')
      attr_list_2[-1] += ';' + attr
  attr_list = attr_list_2

  attr_dict = {
    'id': None,
    'parent': None,
    'name': None,
    'note': None,
  }
  for attr in attr_list:
    if len(attr) == 0:
      continue
    key, value = attr.split('=')
    key = key.strip().lower()
    value = value.strip().strip(';')
    if key in attr_dict:
      attr_dict[key] = value
  return attr_dict


# Problems to fix:
# 1. The 'attr_parent' attributes are missing. These are extracted from the IDs.
# 2. Duplicate exon/intron IDs. We add a number to the ID
#    denoting the ordinal of the exon/intron.
# 3. We only need features of type 'gene', 'mRNA', 'exon', 'intron'.
def fix_oxytri_mac2020_gene_data(data):
  data = data.to_dict('records')
  data = [x for x in data if x['type'] in ['gene', 'mRNA', 'exon', 'intron']]

  feature_counts = {}

  for row in data:
    if row['type'] == 'mRNA':
      row['attr_parent'] = row['attr_id'].split('.')[0];
    elif row['type'] in ['exon', 'intron']:
      id = row['attr_id'];

      # IDs look like cds.t13.g45
      id_prefix, parent = id.split('.', 1)
      
      feature_counts.setdefault(parent, {'exon': 0, 'intron': 0})
      feature_counts[parent][row['type']] += 1

      id_num = feature_counts[parent][row['type']]
      

      id = f'{id_prefix}{id_num}.{parent}'

      row['attr_id'] = id
      row['attr_parent'] = parent
  
  return pd.DataFrame.from_records(data)

# Problems to fix:
# 1. Duplicate exon/intron IDs. We add a number to the ID
#    denoting the ordinal of the exon/intron.
# 2. We only need features of type 'gene', 'mRNA', 'exon', 'intron'.
def fix_oxytri_mac2012_gene_data(data):
  data = data.to_dict('records')
  data = [x for x in data if x['type'] in ['gene', 'mRNA', 'exon', 'intron']]

  feature_counts = {}

  for row in data:
    if row['type'] in ['exon', 'intron']:
      id = row['attr_id']
      parent = row['attr_parent']

      # IDs look like exon_of_ContigXXX.g12.t1
      id_prefix, id_suffix = id.split('_', 1)
      
      feature_counts.setdefault(parent, {'exon': 0, 'intron': 0})
      feature_counts[parent][row['type']] += 1

      id_num = feature_counts[parent][row['type']]

      id = f'{id_prefix}{id_num}_{id_suffix}'

      row['attr_id'] = id
  
  return pd.DataFrame.from_records(data)

def parse_gff(gff_file, nucleus):
  data = pd.read_csv(gff_file, sep='\t', header=None, comment='#')

  if data.shape[1] != 9:
    raise Exception('Unexpected # of columns: ' + str(data.shape[1]))

  data.columns = [
    'contig_name',
    'source',
    'type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'attrs',
  ]

  attr_data = data['attrs'].apply(parse_gff_attrs)
  data = data.drop('attrs', axis='columns')
  attr_data = pd.DataFrame.from_records(attr_data)
  attr_data.columns = 'attr_' + attr_data.columns
  data = pd.concat([data, attr_data], axis='columns')

  data['score'] = data['score'].replace('.', float('nan')).astype(float)
  data['length'] = data['end'] - data['start'] + 1
  data['attr_parent_gene'] = pd.Series([None] * data.shape[0], dtype=str)
  data['contig_nucleus'] = nucleus

  data = data[[
    'contig_name',
    'contig_nucleus',
    'source',
    'type',
    'start',
    'end',
    'length',
    'score',
    'strand',
    'phase',
    'attr_id',
    'attr_parent',
    'attr_name',
    'attr_note',
  ]]
  return data

# Since the GFF format allows features to be hierarchically nested
# usually like gene -> mRNA -> exon/intron, we want each feature
# (ie. exon/intron) to have a handle to it's parent gene, not just
# parent mRNA.
def get_attr_parent_gene(data):
  data = data.to_dict('records')
  
  id_to_parent_gene_id = {}

  # First pass do the genes and mRNA
  for row in data:
    if row['type'] == 'gene':
      id_to_parent_gene_id[row['attr_id']] = row['attr_id']
    elif row['type'] == 'mRNA':
      id_to_parent_gene_id[row['attr_id']] = row['attr_parent']

  # 2nd pass do the exons and introns
  for row in data:
    if row['type'] in ['exon', 'intron']:
      id_to_parent_gene_id[row['attr_id']] = id_to_parent_gene_id[row['attr_parent']]

  # Assign the final parent gene ids
  for row in data:
    row['attr_parent_gene'] = id_to_parent_gene_id[row['attr_id']]

  data = pd.DataFrame.from_records(data)
  return data

# def gff_to_table(gff_file):
#   try:
#     gff = gffutils.create_db(gff_file, get_db_file(gff_file))
#   except:
#     gff = gffutils.FeatureDB(get_db_file(gff_file))
#   data = {
#     # "gene_id": [],
#     # "contig_id": [],
#     "contig_name": [],
#     # "contig_nucleus": [],
#     "source": [],
#     "type": [],
#     "start": [],
#     "end": [],
#     "length": [],
#     "score": [],
#     "strand": [],
#     "phase": [],
#     "attr_id": [],
#     "attr_parent": [],
#     "attr_name": [],
#     "attr_note": [],
#   }
#   for feature in gff.all_features():
#     data['contig_name'].append(feature.seqid)
#     data['source'].append(feature.source)
#     data['type'].append(feature.featuretype)
#     data['start'].append(feature.start)
#     data['end'].append(feature.end)
#     data['length'].append(feature.end - feature.start + 1)
#     data['score'].append(feature.score)
#     data['strand'].append(feature.strand)
#     data['phase'].append(feature.frame)
    
#     attr_keys = list(feature.attributes.keys())
#     attr_values = list(feature.attributes.values())

#     attr_keys = [x.lower() for x in attr_keys]
    
#     attr_dict = dict(zip(attr_keys, attr_values))

#     for key, value in attr_dict.items():
#       if (key != 'note') and (len(value) > 1):
#         raise Exception(f'Duplicate attribute values: key={key} value={value})')
#       else:
#         attr_dict[key] = ','.join(value)
      
#     data['attr_id'].append(attr_dict['id'])
#     data['attr_parent'].append(attr_dict.get('parent', None))
#     data['attr_name'].append(attr_dict.get('name', None))
#     data['attr_note'].append(attr_dict.get('note', None))

#   data = pd.DataFrame(data)
#   data.to_csv(get_csv_file(gff_file), index=False)

# SEE IF CAN USE SQL HERE!
def make_oxytri_mac2020_gene_table():
  data = parse_gff(OXYTRI_MAC_2020_GENE_GFF, 'mac')
  data = fix_oxytri_mac2020_gene_data(data)
  data = get_attr_parent_gene(data)
  data.to_csv(OXYTRI_MAC_2020_GENE_CSV, index=False)

def make_oxytri_mac2012_gene_table():
  data = parse_gff(OXYTRI_MAC_2012_GENE_GFF, 'mac')
  data = fix_oxytri_mac2012_gene_data(data)
  data = get_attr_parent_gene(data)
  data.to_csv(OXYTRI_MAC_2012_GENE_CSV, index=False)
  
  
# make_oxytri_mac2020_gene_table()
make_oxytri_mac2012_gene_table()

# gff_to_table('genes/O_trifallax_2020-upd3.fix.gff3')
# data  = pd.read_csv(get_csv_file('genes/O_trifallax_2020-upd3.fix.gff3'))
# data = data.loc[data['type'] == 'gene']
# data['num_genes'] = data.groupby('contig_name')['contig_name'].transform('count')
# data = data.loc[data['num_genes'] == 1]
# print((data['strand'] == '-').sum())
# print((data['strand'] == '+').sum())
# print(data.loc[data['strand'] == '-'].head())
# exit()

# def get_attrs_dict(attrs_str):
#   attrs_list = attrs_str.split(';')
#   attrs_dict = {}
#   for attr in attrs_list:
#     key, value = attr.split('=')
#     key = key.lower()
#     attrs_dict[RENAME_ATTR_KEY[key]] = value
#   return attrs_dict


# gff_file_name = 'O_trifallax_2020-upd3.fix.gff3'
# fasta_file_name = 'O_trifallax_2020.assembly.fasta'

# gff = pd.read_csv(gff_file_name, header=None, sep='\t')
# gff.columns = ['contig_name', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attrs']

# gff_records = gff.to_dict('records')

# for record in gff_records:
#   record['start'] = int(record['start']) - 1
#   record['end'] = int(record['end']) - 1
#   record['length'] = record['end'] - record['start'] + 1
#   attrs_dict = get_attrs_dict(record['attrs'])
#   del record['attrs']
#   record.update(attrs_dict)
#   record['protein_sequence'] = None

# fasta = SeqIO.index(fasta_file_name, 'fasta')

# mRNA_children = {}
# for record in gff_records:
#   if record['type'] == 'mRNA':
#     mRNA_children[record['attr_id']] = []
#   elif record['type'] == 'exon':
#     mRNA_children[record['attr_parent']].append(record)
#     dna_seq = fasta[record['contig_name']][record['start']:(record['end']+1)]
#     prot_seq = dna_seq.translate()
#     record['protein_sequence'] = prot_seq
