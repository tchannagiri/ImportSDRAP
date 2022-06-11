
import pandas as pd
import mysql.connector


import constants
import common_utils
import mysql_utils

def create_gene_table(to_db: str, file: str):
  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`gene`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`gene` (
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

  gene_data = common_utils.read_tsv(file)

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
  

  gene_data = pd.merge(
    gene_data,
    contig_data,
    how = 'left',
    on = ['contig_name', 'contig_nucleus'],
  )

  mysql_utils.upload_in_chunks(
    gene_data,
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
    f"`{to_db}`.`gene`",
    constants.SQL_BATCH_UPLOAD_ROWS,
  )
  # for i in range(0, data.shape[0], constants.SQL_BATCH_UPLOAD_ROWS):
  #   common_utils.log(f"{i} / {data.shape[0]}")
  #   values = mysql_utils.make_sql_values(
  #     data.loc[
  #       i : (i + constants.SQL_BATCH_UPLOAD_ROWS - 1),
  #       ["id", "name", "alias", "table", "type"]
  #     ]
  #   )
  #   cursor.execute(
  #     f"""
  #     INSERT INTO `gene`
  #     (
  #       `contig_id`,
  #       `contig_name`,
  #       `contig_nucleus`,
  #       `source`,
  #       `type`,
  #       `start`,
  #       `end`,
  #       `length`,
  #       `score`,
  #       `strand`,
  #       `phase`,
  #       `attr_id`,
  #       `attr_parent`,
  #       `attr_name`,
  #       `attr_note`,
  #       `attr_parent_gene`
  #     )
  #     VALUES {values}
  #     """
  #   )

def create_gene_table_2012_mac(to_db: str):
  create_gene_table(to_db, constants.OXYTRI_MAC_2012_GENE_TSV)

def create_gene_table_2020_mac(to_db: str):
  create_gene_table(to_db, constants.OXYTRI_MAC_2020_GENE_TSV)

def create_gene_table_2014_mic(to_db: str):
  create_gene_table(to_db, constants.OXYTRI_MIC_2014_GENE_TSV)


create_gene_table_2012_mac("hello_world")