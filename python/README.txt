This is the script for importing the SDRAP SQL output to a format usable for viewing on <mds_ies_db>

Examples of running the script for the currently hosted databases:
> python python/create_main.py -o mds_ies_db_data_5 -i sdrap_oxy_mac2012_May_30_2022 -p oxytri_mac2012_mic2014 -s all --port 8888
> python python/create_main.py -o mds_ies_db_data_6 -i sdrap_oxy_mac2020_Jun_13_2022 -p oxytri_mac2020_mic2014 -s all --port 8888
> python python/create_main.py -o mds_ies_db_data_7 -i sdrap_ewoo_Sep_29_2022 -p eupwoo_mac2022_mic2022 -s all --port 8888
> python python/create_main.py -o mds_ies_db_data_8 -i sdrap_tet_Sep_28_2022 -p tetsp_mac2015_mic2022 -s all --port 8888

These examples assume that the local port 8888 is forwarded to the server MYSQL port 3306.

Download section on <mds_ies_db>
-------------------------------
After running the script additional steps must be taken to make downloads available on <mds_ies_db>. For example, if the first command above is run to create the database mds_ies_db_data_5, the dumped database files in python/data/mds_ies_db_data_5 must be copied to the directory <mds_ies_db root>/data/mds_ies_db_dev_5/database. The SDRAP output BED and TSV files must be copied to the directory <mds_ies_db root>/data/mds_ies_db_dev_5/annotation.

Adding new organisms
--------------------
Note that some of the files, such as those for gene annotations, proteins, variants, and aliases are premade and kept in the python/data directory in TSV format. In order to incorporate new organisms, these files must be created from the raw GFF/FASTA files for the organism and formatted with the appropriate columns. Please see the examples in the python/data directory. For the script o automatically use these file for the new organism, a new entry in the PRESETS constant in python/constants.py must be added. For sequence searching to work on <mds_ies_db> the BLAST databases must be generated for the new organism using it's FASTA files. The constants in <mds_ies_db root>/constants.php must be updated to reflect the appropriate download files and BLAST assemblies to use.

Old PHP code
------------
The php_old directory is the first version of the import script intended to be run on the main server, or on a local server such as WAMP by forwarding MYSQL ports to the main server. It is kept for now incase it needs to be used for reference.

Connection to MYSQL
-------------------
This script is easiert to run directly on the server. However, to run from a local machine SSH port forwading must be used in order to connect with the main MYSQL server. The port used by the script can be configure on the command line. If additional parameters for connected must be change, they should be directly modified in the get_connection() funciton in mysql_utils.py.