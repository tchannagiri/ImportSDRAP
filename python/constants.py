import os

SQL_BATCH_UPLOAD_ROWS = 1000
SQL_BATCH_DOWNLOAD_ROWS = 10000

DATA_DIR = "python/data"
OXYTRI_MAC_2020_GENE_GFF = os.path.join(DATA_DIR, "O_trifallax_2020-upd3.gff3")
OXYTRI_MAC_2020_GENE_TSV = os.path.join(DATA_DIR, "O_trifallax_2020-upd3.tsv")
OXYTRI_MAC_2012_GENE_GFF = os.path.join(DATA_DIR, "Oxytricha_trifallax_022112.gff3")
OXYTRI_MAC_2012_GENE_TSV = os.path.join(DATA_DIR, "Oxytricha_trifallax_022112.tsv")
OXYTRI_MIC_2014_GENE_GFF = os.path.join(DATA_DIR, "oxy_tri_jrb310_mic_2014.gff")
OXYTRI_MIC_2014_GENE_TSV = os.path.join(DATA_DIR, "oxy_tri_jrb310_mic_2014.tsv")
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