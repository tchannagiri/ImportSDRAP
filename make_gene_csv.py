import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))) # allow importing the utils dir

import pandas as pd
import numpy as np

import common_utils
import constants

def parse_gff_attrs(attr_str):
  attr_list = attr_str.split(";")

  # Use fuzzy logic to check if ";" ocurred in a field, and combine with the previous field
  attr_list_2 = []
  for attr in attr_list:
    if "=" in attr:
      attr_list_2.append(attr)
    else:
      if len(attr_list_2) == 0:
        raise Exception("Can't combine the attribute, no previous one")
      attr_list_2[-1] += ";" + attr
  attr_list = attr_list_2

  attr_dict = {
    "id": None,
    "parent": None,
    "name": None,
    "note": None,
  }
  for attr in attr_list:
    if len(attr) == 0:
      continue
    key, value = attr.split("=")
    key = key.strip().lower()
    value = value.strip().strip(";")
    if key in attr_dict:
      attr_dict[key] = value
  return attr_dict


# Problems to fix:
# 1. The "attr_parent" attributes are missing. These are extracted from the IDs.
# 2. Duplicate exon/intron IDs. We add a number to the ID
#    denoting the ordinal of the exon/intron.
# 3. We only need features of type "gene", "mRNA", "exon", "intron".
def fix_oxytri_mac2020_gene_data(data):
  data = data.to_dict("records")
  data = [x for x in data if x["type"] in ["gene", "mRNA", "exon", "intron"]]

  feature_counts = {}

  for row in data:
    if row["type"] == "mRNA":
      row["attr_parent"] = row["attr_id"].split(".")[0];
    elif row["type"] in ["exon", "intron"]:
      id = row["attr_id"]

      # IDs look like cds.t13.g45
      id_prefix, parent = id.split(".", 1)
      
      feature_counts.setdefault(parent, {"exon": 0, "intron": 0})
      feature_counts[parent][row["type"]] += 1

      id_num = feature_counts[parent][row["type"]]
      

      id = f"{id_prefix}{id_num}.{parent}"

      row["attr_id"] = id
      row["attr_parent"] = parent
  
  return pd.DataFrame.from_records(data)

# Problems to fix:
# 1. Duplicate exon/intron IDs. We add a number to the ID
#    denoting the ordinal of the exon/intron.
# 2. We only need features of type "gene", "mRNA", "exon", "intron".
def fix_oxytri_mac2012_gene_data(data):
  data = data.to_dict("records")
  data = [x for x in data if x["type"] in ["gene", "mRNA", "exon", "intron"]]

  feature_counts = {}

  for row in data:
    if row["type"] in ["exon", "intron"]:
      id = row["attr_id"]
      parent = row["attr_parent"]

      # IDs look like exon_of_ContigXXX.g12.t1
      id_prefix, id_suffix = id.split("_", 1)
      
      feature_counts.setdefault(parent, {"exon": 0, "intron": 0})
      feature_counts[parent][row["type"]] += 1

      id_num = feature_counts[parent][row["type"]]

      id = f"{id_prefix}{id_num}_{id_suffix}"

      row["attr_id"] = id
  
  return pd.DataFrame.from_records(data)

# Problems to fix:
# 1. Change "transcript" type to "mRNA".
def fix_oxytri_mic2014_gene_data(data):
  return data.assign(type=data["type"].replace("transcript", "mRNA"))

def parse_gff(gff_file, nucleus):
  data = pd.read_csv(gff_file, sep="\t", header=None, comment="#")

  if data.shape[1] != 9:
    raise Exception("Unexpected # of columns: " + str(data.shape[1]))

  data.columns = [
    "contig_name",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attrs",
  ]

  attr_data = data["attrs"].apply(parse_gff_attrs)
  data = data.drop("attrs", axis="columns")
  attr_data = pd.DataFrame.from_records(attr_data)
  attr_data.columns = "attr_" + attr_data.columns
  data = pd.concat([data, attr_data], axis="columns")

  data["score"] = data["score"].replace(".", float("nan")).astype(float)
  data["length"] = data["end"] - data["start"] + 1
  data["attr_parent_gene"] = pd.Series([None] * data.shape[0], dtype=str)
  data["contig_nucleus"] = nucleus

  data = data[[
    "contig_name",
    "contig_nucleus",
    "source",
    "type",
    "start",
    "end",
    "length",
    "score",
    "strand",
    "phase",
    "attr_id",
    "attr_parent",
    "attr_name",
    "attr_note",
  ]]
  return data

# Since the GFF format allows features to be hierarchically nested
# usually like gene -> mRNA -> exon/intron, we want each feature
# (ie. exon/intron) to have a handle to it"s parent gene, not just
# parent mRNA.
def get_attr_parent_gene(data):
  data = data.to_dict("records")
  
  id_to_parent_gene_id = {}

  # First pass do the genes and mRNA
  for row in data:
    if row["type"] == "gene":
      id_to_parent_gene_id[row["attr_id"]] = row["attr_id"]
    elif row["type"] == "mRNA":
      id_to_parent_gene_id[row["attr_id"]] = row["attr_parent"]

  # 2nd pass do the exons and introns
  for row in data:
    if row["type"] in ["exon", "intron"]:
      id_to_parent_gene_id[row["attr_id"]] = id_to_parent_gene_id[row["attr_parent"]]

  # Assign the final parent gene ids
  for row in data:
    row["attr_parent_gene"] = id_to_parent_gene_id[row["attr_id"]]

  data = pd.DataFrame.from_records(data)
  return data

def make_oxytri_mac2020_gene_tsv():
  data = parse_gff(constants.OXYTRI_MAC_2020_GENE_GFF, "mac")
  data = fix_oxytri_mac2020_gene_data(data)
  data = get_attr_parent_gene(data)
  file_out = constants.OXYTRI_MAC_2020_GENE_TSV
  common_utils.log(file_out)
  common_utils.write_tsv(data, file_out)

def make_oxytri_mac2012_gene_tsv():
  data = parse_gff(constants.OXYTRI_MAC_2012_GENE_GFF, "mac")
  data = fix_oxytri_mac2012_gene_data(data)
  data = get_attr_parent_gene(data)

  file_out = constants.OXYTRI_MAC_2012_GENE_TSV
  common_utils.log(file_out)
  common_utils.write_tsv(data, file_out)

def make_oxytri_mic2014_gene_tsv():
  data = parse_gff(constants.OXYTRI_MIC_2014_GENE_GFF, "mic")
  data = fix_oxytri_mic2014_gene_data(data)
  data = get_attr_parent_gene(data)
  file_out = constants.OXYTRI_MIC_2014_GENE_TSV
  common_utils.log(file_out)
  common_utils.write_tsv(data, file_out)

def main():
  make_oxytri_mac2020_gene_tsv()
  make_oxytri_mac2012_gene_tsv()
  make_oxytri_mic2014_gene_tsv()

if __name__ == "__main__":
  main()