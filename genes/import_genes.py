import pandas as pd
import numpy as np
import mysql.connector
import time
import os

GENE_DIR = 'genes'
OXYTRI_MAC_2020_GENE_GFF = os.path.join(GENE_DIR, 'O_trifallax_2020-upd3.gff3')
OXYTRI_MAC_2020_GENE_CSV = os.path.join(GENE_DIR, 'O_trifallax_2020-upd3.csv')
OXYTRI_MAC_2012_GENE_GFF = os.path.join(GENE_DIR, 'Oxytricha_trifallax_022112.gff3')
OXYTRI_MAC_2012_GENE_CSV = os.path.join(GENE_DIR, 'Oxytricha_trifallax_022112.csv')
OXYTRI_MIC_2014_GENE_GFF = os.path.join(GENE_DIR, 'oxy_tri_jrb310_mic_2014.gff')
OXYTRI_MIC_2014_GENE_CSV = os.path.join(GENE_DIR, 'oxy_tri_jrb310_mic_2014.csv')

def log(str):
    print(time.strftime("%H:%M:%S", time.localtime()) + ': ' + str)

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

# Problems to fix:
# 1. Change 'transcript' type to 'mRNA'.
def fix_oxytri_mic2014_gene_data(data):
  return data.assign(type=data['type'].replace('transcript', 'mRNA'))

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

def make_oxytri_mic2014_gene_csv():
  data = parse_gff(OXYTRI_MIC_2014_GENE_GFF, 'mic')
  data = fix_oxytri_mic2014_gene_data(data)
  data = get_attr_parent_gene(data)
  data.to_csv(OXYTRI_MIC_2014_GENE_CSV, index=False)

GENE_SQL_COLUMN_SPEC = {
  'gene_id': "int(11) NOT NULL AUTO_INCREMENT COMMENT 'Primary index of the `gene` table'",
  'contig_id': "int(10) NOT NULL COMMENT 'ID of the parent contig in the `contig` table'",
  'contig_name': "varchar(50) NOT NULL COMMENT 'Name of the parent contig'",
  'contig_nucleus': "ENUM('mac','mic') NOT NULL COMMENT 'Type of nucleus of the parent contig'",
  'source': "varchar(50) NOT NULL COMMENT 'Software that produced this annotation'",
  'type': "varchar(50) NOT NULL COMMENT 'Type of feature'",
  'start': "int(11) NOT NULL COMMENT 'Start base pair'",
  'end': "int(11) NOT NULL COMMENT 'End base pair'",
  'length': "int(11) NOT NULL COMMENT 'Length of feature'",
  'score': "float COMMENT 'Score'",
  'strand': "enum('+','-','.') NOT NULL COMMENT 'Orientation of how this feature is read'",
  'phase': "enum('.','0','1','2') NOT NULL COMMENT 'Where the feature begins with reference to the reading frame'",
  'attr_id': "varchar(50) NOT NULL COMMENT 'Unique ID of the feature'",
  'attr_parent': "varchar(50) COMMENT 'ID of the parent feature or NULL if none'",
  'attr_name': "varchar(50) COMMENT 'Human readable name of the feature'",
  'attr_note': "varchar(300) COMMENT 'Additional details about the feature'",
  'attr_parent_gene': "varchar(300) COMMENT 'ID of the gene this feature is on'",
}
GENE_SQL_KEY_SPEC = [
  ('PRIMARY KEY', ('gene_id',)),
  ('KEY', ('attr_id',)),
  ('KEY', ('contig_id',))
]

def get_sql_connection():
  # Note: must have a ssh tunnel setup on your computer to forward port 8888 to the server
  # ssh -N -L 8888:localhost:3306 <username>@knot.math.usf.edu
  return mysql.connector.connect(
    host = 'localhost',
    port = '8888', # To run on the server comment out this line
    user = 'web',
    passwd = 'RZhRwsau6HZrMUXf',
  )

def make_sql_quoted_list(objs, quote):
  return ', '.join(
    'NULL'
    if pd.isna(x) else
    quote + str(x).replace("'", "''") + quote
    for x in objs
  )

def make_sql_values(data, columns):
  data = data[columns]
  # data = data.applymap(
  #   lambda x: "'" + str(x).replace("'", "''") + "'", # escape single quotes
  #   na_action = 'ignore',
  # )
  # data = data.fillna('NULL')
  data = [
    '(' + make_sql_quoted_list(row, "'") + ')'
    for row in data.itertuples(index=False)
  ]
  return data

def make_sql_drop_table(database, table_name):
  return f'DROP TABLE IF EXISTS `{database}`.`{table_name}`;'

def make_sql_create_table(database_name, table_name, column_spec, key_spec, comment):
  column_spec = [f'`{name}` {spec}' for name, spec in column_spec.items()]
  key_spec = [
    f'{key_type} ({make_sql_quoted_list(key_columms, "`")})'
    for key_type, key_columms in key_spec
  ]
  return (
    f"""CREATE TABLE `{database_name}`.`{table_name}` (
      {', '.join(column_spec + key_spec)}
    ) COMMENT='{comment}';"""
  )

def make_sql_column_names(column_names, database_name=None):
  database_prefix = f'`{database_name}`.' if database_name else ''
  return ', '.join(f'{database_prefix}`{name}`' for name in column_names)

def sql_chunk_upload(
  connection,
  cursor,
  database_to,
  table_to,
  column_names,
  sql_values,
  chunk_rows = 10000,
):
  for i in range(0, len(sql_values), chunk_rows):
    log(f'Uploading chunk {i} - {i + chunk_rows - 1}')
    sql_values_chunk = sql_values[i : i + chunk_rows]
    cursor.execute(
      f"""INSERT INTO `{database_to}`.`{table_to}`
      ({make_sql_column_names(column_names)})
      VALUES
      {", ".join(sql_values_chunk)};"""
    )
    connection.commit()

def create_gene_sql_table(connection, cursor, database):
  log(f'Creating `{database}`.`gene`')
  cursor.execute(make_sql_drop_table(database, 'gene'))
  cursor.execute(make_sql_create_table(
    database,
    'gene',
    GENE_SQL_COLUMN_SPEC,
    GENE_SQL_KEY_SPEC,
    'Table of genes and gene features',
  ))

def upload_gene_sql_table_impl(connection, cursor, database, data):
  cursor.execute(make_sql_drop_table(database, 'temp_gene'))
  temp_gene_column_spec = {
    'contig_name': "varchar(50) NOT NULL",
    'contig_nucleus': "ENUM('mac','mic') NOT NULL",
    'source': "varchar(50) NOT NULL",
    'type': "varchar(50) NOT NULL",
    'start': "int(11) NOT NULL",
    'end': "int(11) NOT NULL",
    'length': "int(11) NOT NULL",
    'score': "float",
    'strand': "enum('+','-','.') NOT NULL",
    'phase': "enum('.','0','1','2') NOT NULL",
    'attr_id': "varchar(50) NOT NULL",
    'attr_parent': "varchar(50)",
    'attr_name': "varchar(50)",
    'attr_note': "varchar(300)",
    'attr_parent_gene': "varchar(300)",
  }
  temp_gene_key_spec = [
    ('KEY', ('contig_name',)),
  ]
  temp_gene_columns = list(temp_gene_column_spec)

  cursor.execute(make_sql_create_table(
    database,
    'temp_gene',
    temp_gene_column_spec,
    temp_gene_key_spec,
    '',
  ))
  sql_values = make_sql_values(data, temp_gene_columns)
  sql_chunk_upload(
    connection,
    cursor,
    database,
    'temp_gene',
    temp_gene_columns,
    sql_values,
  )

  cursor.execute(
    f"""INSERT INTO `{database}`.`gene`
    (`contig_id`, {make_sql_column_names(temp_gene_columns)})
    SELECT
      `C`.`contig_id`,
      {make_sql_column_names(temp_gene_columns, "TG")}
    FROM `{database}`.`temp_gene` AS `TG`
    INNER JOIN `{database}`.`contig` AS `C`
    ON `C`.`name` = `TG`.`contig_name`;"""
  )
  connection.commit()

def upload_gene_sql_table(connection, cursor, database, file_name):
  log(f'Uploading `{database}`.`gene`' <- "{file_name}")
  data = pd.read_csv(file_name)
  upload_gene_sql_table_impl(connection, cursor, database, data)

def make_gene_sql_table():
  pass

def make_gene_csvs():
  # make_oxytri_mac2020_gene_table()
  # make_oxytri_mac2012_gene_table()
  make_oxytri_mic2014_gene_csv()
  pass

# make_gene_csvs()

# This is not joining with the reset of the data because of using the OXYTRI_MIC_NAMES
# so... re run SDRAP dammit.
conn = get_sql_connection()
cursor = conn.cursor()
# data = pd.read_csv(OXYTRI_MAC_2020_GENE_CSV)
data = pd.read_csv(OXYTRI_MIC_2014_GENE_CSV)
create_gene_sql_table(conn, cursor, 'mds_ies_db_data_1')
upload_gene_sql_table_impl(conn, cursor, 'mds_ies_db_data_1', data)
cursor.close()
conn.close()