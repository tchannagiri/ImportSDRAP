
import pandas as pd
import mysql.connector


import constants
import common_utils
import mysql_utils

def create_gene_table(db: str):
  common_utils.log(f"create_gene_table {db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)

  cursor.execute(f"DROP TABLE IF EXISTS `{db}`.`gene`;")
  cursor.execute(
    f"""
    CREATE TABLE `{db}`.`gene` (
      `gene_id` int(11) NOT NULL AUTO_INCREMENT COMMENT 'Primary index of the `gene` table',
      `contig_id` int(10) NOT NULL COMMENT 'ID of the parent contig in the `contig` table',
      `contig_name` varchar(50) NOT NULL COMMENT 'Name of the parent contig',
      `contig_nucleus` ENUM('mac','mic') NOT NULL COMMENT 'Type of nucleus of the parent contig',
      `source` varchar(50) NOT NULL COMMENT 'Software that produced this annotation',
      `type` varchar(50) NOT NULL COMMENT 'Type of feature',
      `start` int(11) NOT NULL COMMENT 'Start base pair',
      `end` int(11) NOT NULL COMMENT 'End base pair',
      `length` int(11) NOT NULL COMMENT 'Length of feature',
      `score` float COMMENT 'Score',
      `strand` enum('+','-','.') NOT NULL COMMENT 'Orientation of how this feature is read',
      `phase` enum('.','0','1','2') NOT NULL COMMENT 'Where the feature begins with reference to the reading frame',
      `attr_id` varchar(50) NOT NULL COMMENT 'Unique ID of the feature',
      `attr_parent` varchar(50) COMMENT 'ID of the parent feature or NULL if none',
      `attr_name` varchar(50) COMMENT 'Human readable name of the feature',
      `attr_note` varchar(300) COMMENT 'Additional details about the feature',
      `attr_parent_gene` varchar(300) COMMENT 'ID of the gene this feature is on',
      PRIMARY KEY (`gene_id`),
      KEY (`attr_id`),
      KEY `contig_id` (`contig_id`)
    ) COMMENT='Table of genes and gene features';
    """
  )

  cursor.close()
  conn.close()

def insert_gene_file(to_db: str, file: str):
  common_utils.log(f"insert_gene_file {to_db} {file}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)

  data = common_utils.read_tsv(file)

  cursor.execute(
    f"""
    SELECT
      `contig_id`,
      `nucleus` AS `contig_nucleus`,
      `name` AS `contig_name`
    FROM `{to_db}`.`contig`;
    """
  )

  contig_data = pd.DataFrame.from_records(cursor.fetchall())
  
  data = pd.merge(
    data,
    contig_data,
    how = 'left',
    on = ['contig_name', 'contig_nucleus'],
  )

  mysql_utils.upload_in_chunks(
    data,
    [
      "contig_id",
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
      "attr_parent_gene",
    ],
    cursor,
    to_db,
    "gene",
    constants.SQL_BATCH_UPLOAD_ROWS,
  )

  cursor.close()
  conn.close()
