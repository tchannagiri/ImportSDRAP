import os

SQL_BATCH_UPLOAD_ROWS = 1000
SQL_BATCH_DOWNLOAD_ROWS = 10000

DATA_DIR = "python/data"
OXYTRI_MAC_2020_GENE_GFF = os.path.join(DATA_DIR, "O_trifallax_2020-upd3.gff3")
OXYTRI_MAC_2020_GENE_TSV = os.path.join(DATA_DIR, "O_trifallax_2020-upd3_gene.tsv")
OXYTRI_MAC_2012_GENE_GFF = os.path.join(DATA_DIR, "Oxytricha_trifallax_022112.gff3")
OXYTRI_MAC_2012_GENE_TSV = os.path.join(DATA_DIR, "Oxytricha_trifallax_022112_gene.tsv")
OXYTRI_MIC_2014_GENE_GFF = os.path.join(DATA_DIR, "oxy_tri_jrb310_mic_2014.gff")
OXYTRI_MIC_2014_GENE_TSV = os.path.join(DATA_DIR, "oxy_tri_jrb310_mic_2014_gene.tsv")
OXYTRI_MIC_2014_GENE_HIGHTXN_GFF = os.path.join(DATA_DIR, "tglo_final_hightxn.gff")
OXYTRI_MIC_2014_GENE_LOWXN_GFF = os.path.join(DATA_DIR, "tglo_final_lowtxn.gff")

OXYTRI_MAC_2012_ALIAS = os.path.join(DATA_DIR, "oxy_tri_jrb310_mac_2012_alias.tsv")
OXYTRI_MAC_2020_ALIAS = os.path.join(DATA_DIR, "oxy_tri_jrb310_mac_2020_alias.tsv")
OXYTRI_MIC_2014_ALIAS = os.path.join(DATA_DIR, "oxy_tri_jrb310_mic_2014_alias.tsv")

OXYTRI_MAC_2020_VARIANT = os.path.join(DATA_DIR, "oxy_tri_jrb310_mac_2020_variant.tsv")

OXYTRI_MAC_2012_PROTEIN = os.path.join(DATA_DIR, "Oxytricha_trifallax_022112_aa.tsv")
OXYTRI_MAC_2020_PROTEIN = os.path.join(DATA_DIR, "oxytrijrb310pacbio_aa.tsv")


TABLE_DUMP_LIST = [
  "alias",
  "contig",
  "count",
  "coverage",
  "gene",
  "ies_strict",
  "ies_weak",
  "match",
  "parameter",
  "pointer",
  "properties",
  "protein",
  "stats",
  "variant",
]

PRESETS = {
  "none": {
    "mac_name_regex": ".*", 
    "mic_name_regex": ".*", 
  },
  "oxytri_mac2012_mic2014": {
    "mac_name_regex": "Contig.*", 
    "mic_name_regex": "OXYTRI.*",
    "gene_files": [
      OXYTRI_MAC_2012_GENE_TSV,
      OXYTRI_MIC_2014_GENE_TSV,
    ],
    "alias_files": [
      {"file": OXYTRI_MAC_2012_ALIAS, "table": "contig", "nucleus": "mac"},
      {"file": OXYTRI_MIC_2014_ALIAS, "table": "contig", "nucleus": "mic"},
    ],
    "protein_files": [
      OXYTRI_MAC_2012_PROTEIN
    ],
    "name": "Oxytricha trifallax JRB310 (MAC 2013/MIC 2014)",
    "description": "Oxytricha trifallax JRB310 - SDRAP MDS/IES Annotation - MAC 2013 - MIC 2014",
    "organism": "oxy_tri_jrb310",
    "assembly": "mac_2012,mic_2014",
    "url": "oxy_tri_jrb310_mac_2013_mic_2014",
    "genus": "Oxytricha",
    "species": "Oxytricha trifallax",
    "strain": "JRB310",
    "taxonomy_id": "1172189",
  },
  "oxytri_mac2020_mic2014": {
    "mac_name_regex": "Contig.*", 
    "mic_name_regex": "OXYTRI.*",
    "gene_files": [
      OXYTRI_MAC_2020_GENE_TSV,
      OXYTRI_MIC_2014_GENE_TSV,
    ],
    "alias_files": [
      {"file": OXYTRI_MAC_2020_ALIAS, "table": "contig", "nucleus": "mac"},
      {"file": OXYTRI_MIC_2014_ALIAS, "table": "contig", "nucleus": "mic"},
    ],
    "variant_files": [
      OXYTRI_MAC_2020_VARIANT
    ],
    "protein_files": [
      OXYTRI_MAC_2020_PROTEIN
    ],
    "name": "Oxytricha trifallax JRB310 (MAC 2019/MIC 2014)",
    "description": "Oxytricha trifallax JRB310 - SDRAP MDS/IES Annotation - MAC 2019 - MIC 2014",
    "organism": "oxy_tri_jrb310",
    "assembly": "mac_2020,mic_2014",
    "url": "oxy_tri_jrb310_mac_2019_mic_2014",
    "genus": "Oxytricha",
    "species": "Oxytricha trifallax",
    "strain": "JRB310",
    "taxonomy_id": "1172189",
  },
  "tet": {
    "mac_name_regex": ".*", 
    "mic_name_regex": ".*",
    "gene_files": [],
    "alias_files": [],
    "variant_files": [],
    "protein_files": [],
    "name": "Tetmemena sp. SeJ-2015 (MAC 2015/MIC 2022)",
    "description": "Tetmemena sp. SeJ-2015 - SDRAP MDS/IES Annotation - MAC 2015 - MIC 2022",
    "organism": "tet_sp_sej2015",
    "assembly": "mac_2015,mic_2022",
    "url": "tet_sp_sej2015_mac_2015_mic_2022",
    "genus": "Tetmemena",
    "species": "Tetmemena sp.",
    "strain": "SeJ-2015",
    "taxonomy_id": "200606",
  },
  "ewoo": {
    "mac_name_regex": ".*", 
    "mic_name_regex": ".*",
    "gene_files": [],
    "alias_files": [],
    "variant_files": [],
    "protein_files": [],
    "name": "Euplotes woodruffi Iz01 (MAC 2022/MIC 2022)",
    "description": "Euplotes woodruffi Iz01 - SDRAP MDS/IES Annotation - MAC 2022 - MIC 2022",
    "organism": "eup_woo_iz01",
    "assembly": "mac_2022,mic_2022",
    "url": "eup_woo_iz01_mac_2022_mic_2022",
    "genus": "Euplotes",
    "species": "Euplotes woodruffi",
    "strain": "Iz01",
    "taxonomy_id": "5942",
  },
}

