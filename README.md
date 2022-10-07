# Import scripts for `mds_ies_db` data

Contains scripts for importing [SDRAP](https://github.com/JasperBraun/SDRAP) SQL output and gene annotations into a format usable by [`mds_ies_db`](https://knot.math.usf.edu/mds_ies_db/2022/).

Table of contents:

* [Scripts/data](#scriptsdata)
* [Adding new downloads to `mds_ies_db`](#adding-new-downloads-to-mdsiesdb)
* [Adding new organisms to `mds_ies_db`](#adding-new-organisms-to-mdsiesdb)

## Scripts/data

### `create_main.py`

The main script for importing the SQL output of SDRAP.

Arguments:

* `-o`, `--output_db`: The name of the SQL database where the output will be written.

* `-i`, `--input_db`: Then name of the SQL data where the output of SDRAP is located. The following tables are expected to exist: `alias`, `coverage`, `feature`, `match`, `nucleotide`, `organism`, `parameter`, `pointer`, `properties`, `telomere`. The following tables in SDRAP's output are not used: `gap`, `hsp`. The column format of each table must also match what is expected.

* `-p`, `--preset`: Choices: `oxytri_mac2012_mic2014`, `oxytri_mac2020_mic2014`, `tetsp_mac2015_mic2022`, `eupwoo_mac2022_mic2022`, and `none`. A preset which determines the auxiliary files used for gene annotations, contig aliases, and contig variants. If `none` is given, no auxiliary data will be added to the database. If new organisms or assemblies must be imported using this script, new presets must be added to the `PRESET` variable in `constants.py` and the appropriate auxiliary files must be created.

* `-s`, `--stages`: Choices: `all` (default), `create`, `contig`, `match`, `pointer`, `properties`, `parameter`, `coverage`, `gene`, `ies_strict`, `ies_weak`, `count`, `alias`, `variant`, `stats`, `protein`, `directory`, `dump`, and `drop_temp`. See [stages](#stages) for a description of each stage. By default all stages are run, but you may specify any number of specific stages to be run. Note, due to dependencies between stages, they must be run in the correct order to work correctly. It is recommended to only run all stages unless you know what you are doing.

* `--port`: Port to use for connecting to the MYSQL server. If running on the main server should be kept at the default. If running on a local machine, should be the port which is forwarded to the main server. If additional MYSQL settings must be configured the the `get_connection()` function in `mysql_utils.py` should be modified directly.

#### Stages

The stages of `create_main.py` are:

* `create`: Create the SQL database `OUTPUT_DB`. In the database already exists, it is dropped and recreated.

* `contig`: Create the SQL table `OUTPUT_DB.contig` which contains information about each contig in the assembly. This mainly uses the data in the SDRAP output tables `INPUT_DB.nucleotide`, `INPUT_DB.alias`, and `INPUT_DB.telomere`.

* `match`: Create the SQL table `OUTPUT_DB.match` which contains information about each match (AKA MDS) between a precurcor and product contig. Mainly uses the data in the SDRAP output table `INPUT_DB.match`.

* `pointer`: Create the SQL table `OUTPUT_DB.pointer` which contains information about the pointers (overlapping sergments) of the matches in the `INPUT_DB.match` table. Mainly uses the data in the SDRAP output table `pointer`.

* `properties`: Create the SQL table `OUTPUT_DB.pointer` which contains information about the types of rearrangement maps between precursor and product segments. Mainly uses the data in the SDRAP output table `INPUT_DB.properties`.

* `parameter`: Create the SQL table `OUTPUT_DB.parameter` which contains information about the parameters uses while running SDRAP to create `INPUT_DB`. Mainly uses the data in the SDRAP output table `INPUT_DB.parameter`.

* `coverage`: Create the SQL table `OUTPUT_DB.parameter` which contains information about the parameters uses while running SDRAP to create `INPUT_DB`. Mainly uses the data in the SDRAP output table `INPUT_DB.parameter`.

* `gene`: Create the SQL table `OUTPUT_DB.gene` which contains information about the gene annotations for the assembly. For this table to be populated, a `PRESET` must be used and the appropriate gene TSV file must exists in the `data` directory.

* `ies_strict`: Create the SQL table `OUTPUT_DB.ies_strict` which contains information about the strict IESs of the matches determined SDRAP. Note, this stage is computed in this script (not by SDRAP), although the computations depend on the mathes determined by SDRAP. See [`create_ies.py`](#createiespy) for the definition of a strict IES. Note only strict IESs are currently displayed on `mds_ies_db`, though weak IESs are also provided as downlaod for completeness.

* `ies_weak`: Create the SQL table `OUTPUT_DB.ies_weak` which contains information about the weak IESs of the matches determined SDRAP. Note, this stage is computed in this script (not by SDRAP), although the computations depend on the mathes determined by SDRAP. See [`create_ies.py`](#createiespy) for the definition of a weak IES. Note only strict IESs are currently displayed on `mds_ies_db`, though weak IESs are also provided as downlaod for completeness.

* `count`: Create the SQL table `OUTPUT_DB.count` which contains summary statistics on each contig in the assembly, such as number of genes, number of hits (i.e., other contigs that this one has a match with), number of MDSs (matches), etc. Mainly uses the data in tables created in previous steps: `INPUT_DB.contig`, `INPUT_DB.match`, `INPUT_DB.pointer`, `INPUT_DB.gene`, `INPUT_DB.ies_strict`, and `INPUT_DB.properties`.

* `alias`: Create the SQL table `OUTPUT_DB.alias` which contains the aliases of each contig and/or gene in the assembly. This is used for search functionality so a contig/gene may be searched by its aliases. For this table to be populated, a `PRESET` must be used and the appropriate alias file located in `data`.

* `variant`: Create the SQL table `OUTPUT_DB.variant` which contains the variants of each contig in the assembly (if any). This is used to cross-reference contigs with their variants, though it is currently only relevant for the Oxytricha 2019 MAC assembly. For this table to be populated, a `PRESET` must be used and the appropriate variant file located in `data`.

* `stats`: Create the SQL table `OUTPUT_DB.stats` which contains summary statistics for the whole database, such as the total number of contigs, genes, matches, etc. Mainly uses data from `INPUT_DB.contig`, `INPUT_DB.match`, `INPUT_DB.pointer`, `INPUT_DB.gene`, `INPUT_DB.ies_strict`, and `INPUT_DB.properties`.

* `protein`: Create the SQL table `OUTPUT_DB.protein` which contains amino acid sequence for protein coding sequences identified on contigs. For this table to be populated, a `PRESET` must be used and the appropriate protein file located in `data`.

* `directory`: Adds a row to the SQL table `mds_ies_db.db` which contains information about each of the annotation database currently hosted on `mds_ies_db`. To use this, a `PRESET` should be specified, although it will also work without a present and fill all relevant values with the value of `OUTPUT_DB`. The directory entry is used to determine: name of the SQL database (i.e., value of `OUTPUT_DB`); human readable name/decription of the database; ID of the orgranism and assemblies (e.g., to determine proper BLAST files for sequence search); download directory for the download files for this database; URL string to use for website pages related to this database; genus/species/strain/taxonomy-id for the organism for website display.

* `dump`: Dumps all the tables in `OUTPUT_DB` to provide for download. The output is written to `data/OUTPUT_DB`. Of course, all tables from the previous steps shouls have been created for this to run properly.

* `drop_temp`: Drops the temporary tables created in `OUTPUT_DB` in the previous steps. Currently, these are: `OUTPUT_DB.alias_temp`,  `OUTPUT_DB.name_temp`,  `OUTPUT_DB.protein_temp`,  `OUTPUT_DB.variant_temp`, though not all will have been created if the corresponding tables were not populated.

Note, above `OUTPUT_DB` and `INPUT_DB` refer to the arguments of [`create_main.py`](#createmainpy).

#### Examples

The following commands may be used to reproduce the imports of the current databases on `mds_ies_db`:

```
> python create_main.py -o mds_ies_db_data_5 -i sdrap_oxy_mac2012_May_30_2022 -p oxytri_mac2012_mic2014 -s all --port 8888
> python create_main.py -o mds_ies_db_data_6 -i sdrap_oxy_mac2020_Jun_13_2022 -p oxytri_mac2020_mic2014 -s all --port 8888
> python create_main.py -o mds_ies_db_data_7 -i sdrap_ewoo_Oct_1_2022 -p eupwoo_mac2022_mic2022 -s all --port 8888
> python create_main.py -o mds_ies_db_data_8 -i sdrap_tet_Oct_1_2022 -p tetsp_mac2015_mic2022 -s all --port 8888
```

These examples assume that the local port 8888 is forwarded to the server MYSQL port 3306. For example, with:

```
  ssh -NL 8888:localhost:3306 <username>@knot.math.usf.edu
```


### `make_gene_tsv.py`

Contain scripts for fixing raw GFF files and converting them to tabular for for use by the [`create_main.py`](#createmainpy) script.

### `constants.py`

Contains constants for use by other other scripts. Notably, contains the `PRESETS` variable which must be modified to support new assemblies. Please see the source code for examples of new presets should be created.

### `common_utils.py`

Contains utility functions used by the other scripts.

### `create_gene.py`

Support function for importing gene data into a SQL table.

### `create_ies.py`

Support function for computing strict/weak IESs in the `ies_strict` and `ies_weak` [stage](#stages) of [`create_main.py`](#createmainpy). The definitions of strict/weak IES are:

* Strict IES: The region between two consecutive matches on a precursor contig, such that the two matches are with the same product contig, such that no other matches with other product contigs overlap the region.

* Weak IES: The region between two consecutive matches on a precursor contig, such that the two matches are with the same product contig. However, other matches with other product contigs may overlap the region.

### `interval_utils.py`

Support function for manipulating intervals. Mainly used for [computing IESs](#createiespy).

### `data`

Directory containing auxiliary files for gene annotations, contig aliases, and contig variants. Not contained in the Git repositoy, but should be available on the server in `/home/public/www_safe/channagiri/Ciliate_Project/import_sdrap`.

## Adding new downloads to `mds_ies_db`

After running the [`create_main.py`](#createmainpy) script, additional steps must be taken to make downloads available on `mds_ies_db`. The dumped database files in `data/OUTPUT_DB` must be copied to the directory `WEBSITE_ROOT/data/OUTPUT_DB/database`. The SDRAP output BED and TSV files must be copied to the directory `WEBSITE_ROOT/data/OUTPUT_DB/annotation`. The alias and variant files (if any) must be copied to `WEBSITE_ROOT/data/OUTPUT_DB/annotation`. FASTA nucleotide assembly files must be copied to `WEBSITE_ROOT/data/genome/sequence`; FASTA protein files must be copied to `WEBSITE_ROOT/data/genome/protein`; GFF gene annotation files must be copied to `WEBSITE_ROOT/data/genome/gene`. Finally, the appropriate file names must be added to the `$assemblyInfo` variable `WEBSITE_ROOT/constants.php`. `WEBSITE_ROOT` denotes the location of `mds_ies_db`'s web code, and `OUTPUT_DB` denotes the argument of [`create_main.py`](#createmainpy).

## Adding new organisms to `mds_ies_db`

Some of the files, such as those for gene annotations, proteins, variants, and aliases are premade and kept in the `data` directory in TSV format. In order to incorporate new organisms, these TSV files must be created from the raw GFF/FASTA files for the organism and formatted with the appropriate columns. Please see the examples of these files the `data` directory and the script `make_gene_tsv.py` for how these files are produced.. Once the appropriate TSV files are located in `data`, a new entry in the `PRESETS` constant in [`constants.py`](#createmainpy) must be added. For sequence searching to work on `mds_ies_db` the BLAST databases must be generated for the new organism using it's FASTA files and placed in `WEBSITE_ROOT/blast/data` (see `WEBSITE_ROOT/blast/README.txt` for more details). The `$assemblyInfo` variable in `WEBSITE_ROOT/constants.php` must be updated with a new entry with the appropriate names of the download files and BLAST assemblies to use. `WEBSITE_ROOT` denotes the location of `mds_ies_db`'s web code.
