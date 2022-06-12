import os

SQL_BATCH_UPLOAD_ROWS = 1000
SQL_BATCH_DOWNLOAD_ROWS = 10000

DATA_DIR = 'data'
OXYTRI_MAC_2020_GENE_GFF = os.path.join(DATA_DIR, "O_trifallax_2020-upd3.gff3")
OXYTRI_MAC_2020_GENE_TSV = os.path.join(DATA_DIR, "O_trifallax_2020-upd3.tsv")
OXYTRI_MAC_2012_GENE_GFF = os.path.join(DATA_DIR, "Oxytricha_trifallax_022112.gff3")
OXYTRI_MAC_2012_GENE_TSV = os.path.join(DATA_DIR, "Oxytricha_trifallax_022112.tsv")
OXYTRI_MIC_2014_GENE_GFF = os.path.join(DATA_DIR, "oxy_tri_jrb310_mic_2014.gff")
OXYTRI_MIC_2014_GENE_TSV = os.path.join(DATA_DIR, "oxy_tri_jrb310_mic_2014.tsv")
