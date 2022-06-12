import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))) # allow importing the utils dir

import mysql.connector
import pandas as pd
import numpy as np
import collections

import constants
import common_utils
import interval_utils
import mysql_utils

SQL_BATCH_UPLOAD_ROWS = 1000

def create_name_temp_table(
  to_db: str,
  from_db: str,
  mac_name_regex: str,
  mic_name_regex: str,
):
  common_utils.log(
    f"create_name_temp_table {to_db} {from_db}" +
    f" {mac_name_regex} {mic_name_regex}"
  )

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()
  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`name_temp`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`name_temp` (
      `contig_id` INT NOT NULL,
      `name` VARCHAR(50) NOT NULL,
      PRIMARY KEY (`contig_id`),
      KEY `name` (`name`)
    );
    """
  )

  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`name_temp`
    (
      `contig_id`,
      `name`
    )
    SELECT
      `nuc_id`,
      `alias`
    FROM `{from_db}`.`alias`
    WHERE `alias` REGEXP '{mac_name_regex}'
    OR `alias` REGEXP '{mic_name_regex}';"
    """
  )

  cursor.close()
  conn.close()

def create_contig_table(to_db: str, from_db: str):
  common_utils.log(f"make_contig_table {to_db} {from_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()
  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`contig`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`contig` (
      `contig_id` INT NOT NULL COMMENT 'Index of the `contig` table',
      `nucleus` ENUM('mac','mic') NOT NULL COMMENT 'Type of nucleus',
      `name` VARCHAR(50) NOT NULL COMMENT 'Contig name',
      `length` INT NOT NULL COMMENT 'Length of the sequence',
      `non_tel_length` INT NOT NULL COMMENT 'Length of the sequence between telomeres',
      `telomere_five_start` INT NOT NULL COMMENT 'Start of the 5'' telomere or 0 if there is none',
      `telomere_five_end` INT NOT NULL COMMENT 'End of the 5'' telomere or 0 if there is none',
      `telomere_three_start` INT NOT NULL COMMENT 'Start of the 3'' telomere or 0 if there is none',
      `telomere_three_end` INT NOT NULL COMMENT 'End of the 3'' telomere or 0 if there is none',
      `sequence` longtext NOT NULL COMMENT 'Sequence of nucleotides',
      PRIMARY KEY (`contig_id`),
      UNIQUE KEY `name` (`name`)
    );
    """
  )
  
  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`contig`
    (
      `contig_id`,
      `nucleus`,
      `name`,
      `length`,
      `non_tel_length`,
      `telomere_five_start`,
      `telomere_five_end`,
      `telomere_three_start`,
      `telomere_three_end`,
      `sequence`
    )
    SELECT
      `N`.`nuc_id` AS `contig_id`,
      CASE
        WHEN `F`.`type` = 'prec' THEN 'mic'
        WHEN `F`.`type` = 'prod' THEN 'mac'
      END AS `nucleus`,
      `A`.`name` AS `name`,
      LENGTH(`N`.`sequence`) AS `length`,
      (
        IF(`T`.`three_length` = 0 OR `T`.`three_length` IS NULL, LENGTH(`N`.`sequence`), `T`.`three_start` - 1) -
        IF(`T`.`five_length` = 0 OR `T`.`five_length` IS NULL, 1, `T`.`five_start` + `T`.`five_length`) +
        1
      ) AS `non_tel_length`,
      IF(`T`.`five_length` = 0 OR `T`.`five_length` IS NULL, 0, `T`.`five_start`) AS `telomere_five_start`,
      IF(`T`.`five_length` = 0 OR `T`.`five_length` IS NULL, 0, `T`.`five_start` + `T`.`five_length` - 1) AS `telomere_five_end`,
      IF(`T`.`three_length` = 0 OR `T`.`three_length` IS NULL, 0, `T`.`three_start`) AS `telomere_three_start`,
      IF(`T`.`three_length` = 0 OR `T`.`three_length` IS NULL, 0, `T`.`three_start` + `T`.`three_length` - 1) AS `telomere_three_end`,
      UPPER(`N`.`sequence`) AS `sequence`
    FROM `{from_db}`.`nucleotide` AS `N`
    LEFT JOIN `{to_db}`.`name_temp` AS `A`
      ON `N`.`nuc_id` = `A`.`contig_id`
    LEFT JOIN `{from_db}`.`feature` AS `F`
      ON `N`.`feat_id` = `F`.`feat_id`
    LEFT JOIN `{from_db}`.`telomere` AS `T`
      ON `N`.`nuc_id` = `T`.`nuc_id`;
    """
  )
  cursor.close()
  conn.close()

def create_alias_table(
  to_db: str,
  from_db: str,
  alias_file_list: list[str],
):
  common_utils.log(f"create_alias {to_db} {from_db} {alias_file_list}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)
  cursor.execute(f"SELECT * FROM `{from_db}`.`alias`;")
  alias_from = cursor.fetchall()
  alias_from = pd.DataFrame.from_records(alias_from)

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`alias`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`alias` (
      `alias_id` INT NOT NULL AUTO_INCREMENT COMMENT 'primary key for the table',
      `id` INT NOT NULL COMMENT 'primary key of the table for the corresponding sequence',
      `name` VARCHAR(50) NOT NULL COMMENT 'the primary name of the sequence',
      `alias` VARCHAR(50) NOT NULL COMMENT 'alias of the sequence',
      `table` VARCHAR(50) NOT NULL COMMENT 'the table that `id` refers to',
      `type` VARCHAR(50) DEFAULT NULL COMMENT 'the type of alias',
      PRIMARY KEY (`alias_id`),
      KEY `name` (`name`),
      KEY `alias` (`alias`)
    ) COMMENT='table of sequence aliases';
    """
  )

  for alias_file in alias_file_list:
    common_utils.log(alias_file)
    alias_to = pd.read_csv(alias_file)
    alias_to = alias_to.merge(
      alias_from[["nuc_id", "alias"]],
      left_on = "primary",
      right_on = "alias",
      how = "inner",
    )
    alias_to = alias_to.drop("alias", axis="columns")
    alias_to = alias_to.melt(
      id_vars = ["nuc_id", "primary"],
      var_name = "type",
      value_name = "alias",
    )
    alias_to = alias_to.dropna(axis='index')
    alias_to = alias_to.rename({"nuc_id": "id", "primary": "name"}, axis="columns")
    alias_to["table"] = "contig"
    mysql_utils.upload_in_chunks(
      alias_to,
      ["id", "name", "alias", "table", "type"],
      cursor,
      to_db,
      "alias",
      constants.SQL_BATCH_UPLOAD_ROWS,
    )
    # for i in range(0, alias_to.shape[0], SQL_BATCH_UPLOAD_ROWS):
    #   common_utils.log(f"{i} / {alias_to.shape[0]}")
    #   values = mysql_utils.make_sql_values(
    #     alias_to.loc[
    #       i : (i + SQL_BATCH_UPLOAD_ROWS - 1),
    #       ["id", "name", "alias", "table", "type"]
    #     ])
    #   cursor.execute(
    #     f"""
    #     INSERT INTO `{to_db}`.`alias`
    #     (
    #       `id`,
    #       `name`,
    #       `alias`,
    #       `table`,
    #       `type`
    #     )
    #     VALUES {values};
    #     """
    #   )
  cursor.close()
  conn.close()

def create_match_table(to_db: str, from_db: str):
  common_utils.log(f"create_match_table {to_db} {from_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()
  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`match`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`match` (
      `match_id` INT NOT NULL COMMENT 'primary key for the table',
      `mac_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the corresponding precursor sequence',
      `mac_name` VARCHAR(50) NOT NULL COMMENT 'MAC contig name',
      `mic_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the corresponding product sequence',
      `mic_name` VARCHAR(50) NOT NULL COMMENT 'MIC contig name',
      `mic_start` INT NOT NULL COMMENT 'position of the first base pair of the precursor segment of the match on the corresponding precursor sequence',
      `mic_end` INT NOT NULL COMMENT 'position of the last base pair of the precursor segment of the match on the corresponding precursor sequence',
      `mac_start` INT NOT NULL COMMENT 'position of first base pair of the product segment of the match on the corresponding product sequence',
      `mac_end` INT NOT NULL COMMENT 'position of last base pair of the product segment of the match on the corresponding product sequence',
      `orientation` ENUM('+','-') CHARACTER SET latin1 NOT NULL COMMENT 'orientation of the match',
      `length` INT NOT NULL COMMENT 'length of the match',
      `index` INT NOT NULL COMMENT 'index of the match or fragment in the arrangement of the corresponding product on the corresponding precursor',
      `pre_cov` DECIMAL(8,2) NOT NULL COMMENT 'percentage of product segment of preliminary match covered by product segment of additional match',
      `add_cov` DECIMAL(8,2) NOT NULL COMMENT 'percentage of product segment of additional match covered by product segment of preliminary match',
      `is_preliminary` BOOLEAN NOT NULL COMMENT 'indicates whether the match is a preliminary match (1 means yes/true, 0 means no/false)',
      `is_additional` BOOLEAN NOT NULL COMMENT 'indicates whether the match is an additional match (1 means yes/true, 0 means no/false)',
      `is_fragment` BOOLEAN NOT NULL COMMENT 'indicates whether the match is a fragment (1 means yes/true, 0 means no/false)',
      PRIMARY KEY (`match_id`),
      KEY `mac_contig_id` (`mac_contig_id`),
      KEY `mic_contig_id` (`mic_contig_id`)
    );
    """
  )

  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`match`
    (
      `match_id`,
      `mac_contig_id`,
      `mac_name`,
      `mic_contig_id`,
      `mic_name`,
      `mic_start`,
      `mic_end`,
      `mac_start`,
      `mac_end`,
      `orientation`,
      `length`,
      `index`,
      `pre_cov`,
      `add_cov`,
      `is_preliminary`,
      `is_additional`,
      `is_fragment`
    )
    SELECT
      `M`.`match_id`,
      `M`.`prod_nuc_id`,
      `NMAC`.`name`,
      `M`.`prec_nuc_id`,
      `NMIC`.`name`,
      `M`.`prec_start`,
      `M`.`prec_end`,
      `M`.`prod_start`,
      `M`.`prod_end`,
      `M`.`orientation`,
      `M`.`length`,
      ABS(`M`.`index`),
      `M`.`pre_cov`,
      `M`.`add_cov`,
      `M`.`is_preliminary`,
      `M`.`is_additional`,
      `M`.`is_fragment`
    FROM `{from_db}`.`match` AS `M`
    LEFT JOIN `{to_db}`.`contig` AS `NMIC`
      ON `M`.`prec_nuc_id` = `NMIC`.`contig_id`
    LEFT JOIN `{to_db}`.`contig` AS `NMAC`
      ON `M`.`prod_nuc_id` = `NMAC`.`contig_id`;
    """
  )
  cursor.close()
  conn.close()

def create_pointer_table(to_db: str, from_db: str):
  common_utils.log(f"create_pointer_table {to_db} {from_db}")
  
  conn = mysql_utils.get_connection()
  cursor = conn.cursor()
  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`pointer`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`pointer` (
      `ptr_id` INT NOT NULL COMMENT 'primary key for the table',
      `mic_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the corresponding precursor sequence',
      `mic_name` VARCHAR(50) NOT NULL COMMENT 'MIC contig name',
      `mac_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the corresponding product sequence',
      `mac_name` VARCHAR(50) NOT NULL COMMENT 'MAC contig name',
      `left_match_id` INT NOT NULL COMMENT 'primary key of the `match` table for the corresponding match flanking the pointer on the left on the corresponding product sequence',
      `right_match_id` INT NOT NULL COMMENT 'primary key of the `match` table for the corresponding match flanking the pointer on the right on the corresponding product sequence',
      `left_index` INT NOT NULL COMMENT 'index of the match on the left of this pointer',
      `right_index` INT NOT NULL COMMENT 'index of the match on the right of this pointer',
      `mac_start` INT NOT NULL COMMENT 'position of the first base pair of the pointer on the corresponding product sequence',
      `mac_end` INT NOT NULL COMMENT 'position of the last base pair of the pointer on the corresponding product sequence',
      `left_mic_start` INT NOT NULL COMMENT 'position of the first base pair of the pointer at the end of the precursor segment of the left match on the corresponding precursor sequence',
      `left_mic_end` INT NOT NULL COMMENT 'position of the last base pair of the pointer at the end of the precursor segment of the left match on the corresponding precursor sequence',
      `left_match_orientation` ENUM('+','-') CHARACTER SET latin1 COMMENT 'orientation of left match',
      `right_mic_start` INT NOT NULL COMMENT 'position of the first base pair of the pointer at the end of the precursor segment of the right match on the corresponding precursor sequence',
      `right_mic_end` INT NOT NULL COMMENT 'position of the last base pair of the pointer at the end of the precursor segment of the right match on the corresponding precursor sequence',
      `right_match_orientation` ENUM('+','-') CHARACTER SET latin1 COMMENT 'orientation of right match',
      `length` INT NOT NULL COMMENT 'length of the pointer',
      `is_preliminary` BOOLEAN COMMENT 'indicates whether the pointer is computed purely from preliminary matches (1 means yes/true, 0 means no/false)',
      PRIMARY KEY (`ptr_id`),
      KEY `mac_contig_id` (`mac_contig_id`),
      KEY `mic_contig_id` (`mic_contig_id`)
    );
    """
  )

  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`pointer`
    (
      `ptr_id`,
      `mic_contig_id`,
      `mic_name`,
      `mac_contig_id`,
      `mac_name`,
      `left_match_id`,
      `right_match_id`,
      `left_index`,
      `right_index`,
      `mac_start`,
      `mac_end`,
      `left_mic_start`,
      `left_mic_end`,
      `left_match_orientation`,
      `right_mic_start`,
      `right_mic_end`,
      `right_match_orientation`,
      `length`,
      `is_preliminary`
    )
    SELECT
      `P`.`ptr_id`,
      `P`.`prec_nuc_id`,
      `NMIC`.`name`,
      `P`.`prod_nuc_id`,
      `NMAC`.`name`,
      `P`.`left_match_id`,
      `P`.`right_match_id`,
      `MLEFT`.`index`,
      `MRIGHT`.`index`,
      `P`.`prod_start`,
      `P`.`prod_end`,
      `P`.`left_prec_start`,
      `P`.`left_prec_end`,
      `P`.`left_match_orientation`,
      `P`.`right_prec_start`,
      `P`.`right_prec_end`,
      `P`.`right_match_orientation`,
      `P`.`length`,
      `P`.`is_preliminary`
    FROM `{from_db}`.`pointer` AS `P`
    LEFT JOIN `{to_db}`.`contig` AS `NMIC`
      ON `P`.`prec_nuc_id` = `NMIC`.`contig_id`
    LEFT JOIN `{to_db}`.`contig` AS `NMAC`
      ON `P`.`prod_nuc_id` = `NMAC`.`contig_id`
    LEFT JOIN `{to_db}`.`match` AS `MLEFT`
      ON `P`.`left_match_id` = `MLEFT`.`match_id`
    LEFT JOIN `{to_db}`.`match` AS `MRIGHT`
      ON `P`.`right_match_id` = `MRIGHT`.`match_id`;
    """
  )

def create_properties_table(
  to_db: str,
  from_db: str,
):
  common_utils.log(f"create_properties_table {to_db} {from_db}")
  conn = mysql_utils.get_connection()
  cursor = conn.cursor()
  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`properties`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`properties` (
      `prop_id` INT NOT NULL COMMENT 'primary key of the table',
      `mic_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the corresponding precursor sequence',
      `mic_name` VARCHAR(50) NOT NULL COMMENT 'MIC contig name',
      `mac_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the corresponding product sequence',
      `mac_name` VARCHAR(50) NOT NULL COMMENT 'MAC contig name',
      `preliminary_match_number` INT NOT NULL COMMENT 'number of preliminary matches of the corresponding product sequence on the corresponding precursor sequence',
      `total_match_number` INT NOT NULL COMMENT 'total number of matches of the corresponding prduct sequence on the corresponding precursor sequence',
      `non_gapped` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement has no gaps (1 means yes/true, 0 means no/false)',
      `non_overlapping` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement has no overlaps (1 means yes/true, 0 means no/false)',
      `non_repeating` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement has no repeats (1 means yes/true, 0 means no/false)',
      `exceeded_clique_limit` BOOLEAN NOT NULL COMMENT 'indicates whether the number of non repeating and non p overlapping subarrangements of the arrangement exceeded the limit found in the `parameter` table in the `property_clique_limit` column (1 means yes/true, 0 means no/false)',
      `weakly_complete` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is weakly complete (1 means yes/true, 0 means no/false)',
      `strongly_complete` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is strongly complete (1 means yes/true, 0 means no/false)',
      `weakly_consecutive` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is weakly consecutive (1 means yes/true, 0 means no/false)',
      `strongly_consecutive` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is strongly consecutive (1 means yes/true, 0 means no/false)',
      `weakly_ordered` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is weakly ordered (1 means yes/true, 0 means no/false)',
      `strongly_ordered` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is strongly ordered (1 means yes/true, 0 means no/false)',
      `weakly_non_scrambled` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is weakly non scrambled (1 means yes/true, 0 means no/false)',
      `strongly_non_scrambled` BOOLEAN NOT NULL COMMENT 'indicates whether the arrangement is strongly non scrambled (1 means yes/true, 0 means no/false)',
      PRIMARY KEY (`prop_id`),
      KEY `mac_contig_id` (`mac_contig_id`),
      KEY `mic_contig_id` (`mic_contig_id`)
    ) COMMENT='table of properties for each arrangement with coverage above the threshold found in the `parameter` table in the `property_min_coverage` column';
    """
  )

  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`properties`
    (
      `prop_id`,
      `mic_contig_id`,
      `mic_name`,
      `mac_contig_id`,
      `mac_name`,
      `preliminary_match_number`,
      `total_match_number`,
      `non_gapped`,
      `non_overlapping`,
      `non_repeating`,
      `exceeded_clique_limit`,
      `weakly_complete`,
      `strongly_complete`,
      `weakly_consecutive`,
      `strongly_consecutive`,
      `weakly_ordered`,
      `strongly_ordered`,
      `weakly_non_scrambled`,
      `strongly_non_scrambled`
    )
    SELECT
      `P`.`prop_id`,
      `P`.`prec_nuc_id`,
      `NMIC`.`name`,
      `P`.`prod_nuc_id`,
      `NMAC`.`name`,
      `P`.`preliminary_match_number`,
      `P`.`total_match_number`,
      `P`.`non_gapped`,
      `P`.`non_overlapping`,
      `P`.`non_repeating`,
      `P`.`exceeded_clique_limit`,
      `P`.`weakly_complete`,
      `P`.`strongly_complete`,
      `P`.`weakly_consecutive`,
      `P`.`strongly_consecutive`,
      `P`.`weakly_ordered`,
      `P`.`strongly_ordered`,
      `P`.`weakly_non_scrambled`,
      `P`.`strongly_non_scrambled`
    FROM `{from_db}`.`properties` AS `P`
    LEFT JOIN `{to_db}`.`contig` AS `NMIC`
      ON `P`.`prec_nuc_id` = `NMIC`.`contig_id`
    LEFT JOIN `{to_db}`.`contig` AS `NMAC`
      ON `P`.`prod_nuc_id` = `NMAC`.`contig_id`;
    """
  )


def create_parameter_table(to_db: str, from_db: str):
  common_utils.log(f"create_parameter_table {to_db} {from_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`parameter`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`parameter` (
      `parameter_id` INT NOT NULL COMMENT 'primary key for the table',
      `date` TEXT COMMENT 'date of running the pipeline with specified parameter values',
      `database` TEXT COMMENT 'name of the mysql database',
      `username` TEXT COMMENT 'mysql username used to run the pipeline',
      `genus` TEXT COMMENT 'name of the genus containing the species whose genomes are contained in this database',
      `species` TEXT COMMENT 'name of the species whose genomes are contained in this database',
      `strain` TEXT COMMENT 'strain of the species whose genomes are contained in this database',
      `taxonomy_id` TEXT COMMENT 'taxonomy id of the species whose genomes are contained in this database',
      `precursor_filename` TEXT COMMENT 'path to the file containing the precursor genome contained in this database',
      `precursor_delimiter` TEXT COMMENT 'delimiter in description lines of precursor genome file',
      `product_filename` TEXT COMMENT 'path to the file containing the product genome contained in this database',
      `product_delimiter` TEXT COMMENT 'delimiter in description lines of product genome file',
      `telo_pattern` TEXT COMMENT 'repeated DNA sequence characterizing the telomere pattern of the species whose genomes are contained in this database',
      `telo_error_limit` DECIMAL(4,3) DEFAULT NULL COMMENT 'maximum number of erroneous basepairs relative to the current length of telomere allowed during telomere expansion',
      `telo_buffer_limit` INT DEFAULT NULL COMMENT 'maximum number of erronerous basepairs in excess of non erroneous basepairs allowed in buffer during telomere expansion',
      `telo_max_length` INT DEFAULT NULL COMMENT 'maximum length allowed for telomeres',
      `telo_max_offset` INT DEFAULT NULL COMMENT 'maximum distance allowed between telomere and end of contig',
      `telo_min_length` INT DEFAULT NULL COMMENT 'minimum length require_onced for telomeres',
      `hsp_min_length` INT DEFAULT NULL COMMENT 'minimum length require_onced for HSPs',
      `pre_match_min_bitscore` INT DEFAULT NULL COMMENT 'minimum biscore require_onced for HSPs used for preliminary arrangement',
      `pre_match_min_pident` DECIMAL(5,2) DEFAULT NULL COMMENT 'minimum percent identity require_onced for HSPs used for preliminary arrangement',
      `pre_match_min_coverage_addition` INT DEFAULT NULL COMMENT 'minimum number of basepairs an HSP must cover and which are not already covered to be considered to become a preliminary match during computation of preliminary arrangement',
      `merge_tolerance` INT DEFAULT NULL COMMENT 'maximum shift allowed between precursor segments and between product segments of two matches to still be considered for merging',
      `merge_max_gap` INT DEFAULT NULL COMMENT 'maximum gap allowed between precursor segments and between product segments of two matches to still be considered for merging',
      `gap_min_length` INT DEFAULT NULL COMMENT 'minimum length of non covered basepairs require_onced for a gap to be annotated',
      `pointer_min_length` INT DEFAULT NULL COMMENT 'minimum length of overlap between product segments of preliminary matches require_onced for a pointer to be annotated',
      `add_match_min_bitscore` INT DEFAULT NULL COMMENT 'minimum bitscore require_onced for HSPs used for additional matches',
      `add_match_min_pident` DECIMAL(5,2) DEFAULT NULL COMMENT 'minimum percent identity require_onced for HSPs used for additional matches',
      `add_match_min_prod_segment_overlap` DECIMAL(4,3) DEFAULT NULL COMMENT 'minimum portion of product segment of preliminary match contained in product segment of additional match require_onced for additional match to be annotated as match with the index of that preliminary match',
      `fragment_min_prod_segment_overlap` DECIMAL(4,3) DEFAULT NULL COMMENT 'minimum portion of product segment of preliminary match contained in product segment of HSP require_onced for HSP to be annotated as fragment with the index of that preliminary match',
      `property_min_coverage` INT DEFAULT NULL COMMENT 'minimum portion of part of product sequence between telomeres covered by product segments of preliminary matches of a precursor sequence require_onced for properties of the arrangement between the two sequences to be computed',
      `property_max_match_overlap` INT DEFAULT NULL COMMENT 'maximum overlap allowed between precursor segments of matches in an arrangement for the two matches to still be considered disjoint (if overlap as long as one of the precursor segments the two matches are not considered disjoint regardless of value of this parameter)',
      `property_clique_limit` INT DEFAULT NULL COMMENT 'maximum number of maximal cliques in match graph of an arrangement for which properties of the corresponding subarrangements are computed',
      `scr_complete` BOOLEAN DEFAULT NULL COMMENT 'indicates whether or not completeness is require_onced for a non repeating and non p-overlapping subarrangement to be considered non scrambled (1 means yes/true, 0 means no/false)',
      `scr_consecutive` BOOLEAN DEFAULT NULL COMMENT 'indicates whether or not consecutivenss is require_onced for a non repeating and non p-overlapping subarrangement to be considered non scrambled (1 means yes/true, 0 means no/false)',
      `scr_ordered` BOOLEAN DEFAULT NULL COMMENT 'indicates whether or not orderedness is require_onced for a non repeating and non p-overlapping subarrangement to be considered non scrambled (1 means yes/true, 0 means no/false)',
      `output_min_coverage` INT DEFAULT NULL COMMENT 'minimum portion of part of product sequence between telomeres covered by product segments of preliminary matches of a precursor sequence require_onced for annotations of precursor segments, product segments, fragments, gaps, and pointers to be included in output files',
      `output_give_summary` BOOLEAN DEFAULT NULL COMMENT 'indicates whether the a summary of interesting values should be included in the output (1 means yes/true, 0 means no/false)',
      `output_use_alias` BOOLEAN DEFAULT NULL COMMENT 'indicates whether the program should reference sequences in its output by the identifier assigned by the program (1 means yes/true, 0 means no/false)',
      PRIMARY KEY (`parameter_id`)
    ) COMMENT='table of user defined parameters used when running the pipeline';
    """
  )
  cursor.execute(f"INSERT INTO `{to_db}`.`parameter` SELECT * from `{from_db}`.`parameter`;")

def create_coverage_table(to_db: str, from_db: str):
  common_utils.log(f"create_coverage_table {to_db} {from_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()
  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`coverage`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`coverage` (
      `cov_id` INT NOT NULL COMMENT 'primary key for the table',
      `mic_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the MIC contig',
      `mic_name` VARCHAR(50) NOT NULL COMMENT 'MIC contig name',
      `mac_contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the MAC contig',
      `mac_name` VARCHAR(50) NOT NULL COMMENT 'MAC contig name',
      `mic_start` INT NOT NULL COMMENT 'Minimum MDS start position on the MIC contig',
      `mic_end` INT NOT NULL COMMENT 'Maximum MDS end position on the MIC contig',
      `mac_start` INT NOT NULL COMMENT 'Minimum MDS start position on the MAC contig',
      `mac_end` INT NOT NULL COMMENT 'Maximum MDS end position on the MAC contig',
      `mac_coverage_percent` DECIMAL(5,2) NOT NULL COMMENT 'percentage of MAC contig between the telomeres covered by MDSs of the two contigs',
      `mac_coverage_length` INT NOT NULL COMMENT 'Total length that the MDSs cover on the MAC contig',
      `mds_num` INT NOT NULL COMMENT 'Number of MDSs of the two contigs',
      `pointer_num` INT NOT NULL COMMENT 'Number of pointers of MDSs of the two contigs',
      `ies_num` INT NOT NULL COMMENT 'Number of strict of IESs of MDSs of the two contigs',
      PRIMARY KEY (`cov_id`),
      KEY `mac_contig_id` (`mac_contig_id`),
      KEY `mic_contig_id` (`mic_contig_id`)
    ) COMMENT='Summary data about the matches between a pair of contig.';
    """
  )

  cursor.execute(f"SELECT MAX(`mac_contig_id`) AS `max` FROM `{to_db}`.`properties`;")
  mac_id_max = cursor.fetchall()[0][0]
  mac_id_batch_size = 100 # How many mac contigs to process at a time
  cov_id = 1 # Current coverage id
  for mac_id_batch_start in range(0, mac_id_max + 1, mac_id_batch_size):
    common_utils.log(f"{mac_id_batch_start} / {mac_id_max}")
    mac_id_batch_end = min(mac_id_batch_start + mac_id_batch_size - 1, mac_id_max)

    # Need to reopen the connect to avoid "Commands out of sync" error
    # This might be happening because the loop interleaves SELECT and INSERT statements
    cursor.close()
    conn.close()
    conn = mysql_utils.get_connection()
    cursor = conn.cursor(dictionary=True)
    cursor.execute(
      f"""
      SELECT
        `contig_id`,
        `non_tel_length`
      FROM
        `{to_db}`.`contig`
      WHERE
        `contig_id` BETWEEN {mac_id_batch_start} AND {mac_id_batch_end};
      """
    )
    mac_length_dict = {
      int(x["contig_id"]): int(x["non_tel_length"])
      for x in cursor.fetchall()
    }

    cursor.execute(
      f"""
      SELECT
        `mac_contig_id`,
        `mic_contig_id`,
        `mac_start`,
        `mac_end`,
        `mic_start`,
        `mic_end`
      FROM
        `{to_db}`.`match`
      WHERE
        `mac_contig_id` BETWEEN {mac_id_batch_start} AND {mac_id_batch_end} 
      AND
        `is_preliminary` = 1;
      """
    )

    mac_mds_start_dict = collections.defaultdict(lambda: np.inf) # Starting position of any MDS 
    mac_mds_end_dict = collections.defaultdict(lambda: -np.inf) # End position of any MDS
    mic_mds_start_dict = collections.defaultdict(lambda: np.inf) # Starting position of any MDS
    mic_mds_end_dict = collections.defaultdict(lambda: -np.inf) # Starting position of any MDS
    mds_dict = collections.defaultdict(list) # List of MDS's for the MAC/MIC pair
    for mds in cursor.fetchall():
      mds["mac_contig_id"] = int(mds["mac_contig_id"])
      mds["mic_contig_id"] = int(mds["mic_contig_id"])
      mds["mac_start"] = int(mds["mac_start"])
      mds["mac_end"] = int(mds["mac_end"])
      mds["mic_start"] = int(mds["mic_start"])
      mds["mic_end"] = int(mds["mic_end"])

      key = (mds['mac_contig_id'], mds['mic_contig_id'])
      mac_mds_start_dict[key] = min(mac_mds_start_dict[key], mds['mac_start'])
      mac_mds_end_dict[key] = max(mac_mds_end_dict[key], mds['mac_end'])
      mic_mds_start_dict[key] = min(mic_mds_start_dict[key], mds['mic_start'])
      mic_mds_end_dict[key] = max(mic_mds_end_dict[key], mds['mic_end'])
      mds_dict[key].append({
        'start': mds['mac_start'],
        'end': mds['mac_end'],
      })

    sql_table_rows = '' # Table rows in SQL syntax
    for (mac_id, mic_id), mds_list in mds_dict.items():
      mds_list = interval_utils.get_union(mds_list)
      coverageLength = sum(x['end'] - x['start'] + 1 for x in mds_list)
      
      coveragePercent = 100 * coverageLength / mac_length_dict[mac_id]
      mac_start = mac_mds_start_dict[mac_id, mic_id]
      mac_end = mac_mds_end_dict[mac_id, mic_id]
      mic_start = mic_mds_start_dict[mac_id, mic_id]
      mic_end = mic_mds_end_dict[mac_id, mic_id]
      if len(sql_table_rows) == 0:
        sql_table_rows = (
          "SELECT " +
            f"{cov_id} AS `cov_id`, " +
            f"{mac_id} AS `mac_contig_id`, " +
            f"{mac_start} AS `mac_start`, " +
            f"{mac_end} AS `mac_end`, " +
            f"{mic_id} AS `mic_contig_id`, " +
            f"{mic_start} AS `mic_start`, " +
            f"{mic_end} AS `mic_end`, " +
            f"{coverageLength} AS `length`, " +
            f"{coveragePercent} AS `coverage` "
        )
      else:
        sql_table_rows += (
          "UNION ALL SELECT "
            f"{cov_id}, " +
            f"{mac_id}, " +
            f"{mac_start}, " +
            f"{mac_end}, " +
            f"{mic_id}, " +
            f"{mic_start}, " +
            f"{mic_end}, " +
            f"{coverageLength}, " +
            f"{coveragePercent} "
        )
      cov_id += 1

    cursor.execute(
      f"""
      INSERT INTO `{to_db}`.`coverage`
      (
        `cov_id`,
        `mic_contig_id`,
        `mic_name`,
        `mic_start`,
        `mic_end`,
        `mac_contig_id`,
        `mac_name`,
        `mac_start`,
        `mac_end`,
        `mac_coverage_percent`,
        `mac_coverage_length`,
        `mds_num`,
        `pointer_num`,
        `ies_num`
      )
      SELECT
        `T`.`cov_id` AS `cov_id`,
        `T`.`mic_contig_id` AS `mic_contig_id`,
        `NMIC`.`name` AS `mic_name`,
        `T`.`mic_start` AS `mic_start`,
        `T`.`mic_end` AS `mic_end`,
        `T`.`mac_contig_id` AS `mac_contig_id`,
        `NMAC`.`name` AS `mac_name`,
        `T`.`mac_start` AS `mac_start`,
        `T`.`mac_end` AS `mac_end`,
        `T`.`coverage` AS `mac_coverage_percent`,
        `T`.`length` AS `mac_coverage_length`,
        0 AS `mds_num`,
        0 AS `pointer_num`,
        0 AS `ies_num`
      FROM
        ({sql_table_rows}) AS `T`
      LEFT JOIN
        `{to_db}`.`contig` AS `NMAC`
      ON
        `T`.`mac_contig_id` = `NMAC`.`contig_id`
      LEFT JOIN
        `{to_db}`.`contig` AS `NMIC`
      ON
        `T`.`mic_contig_id` = `NMIC`.`contig_id`;"
      """
    )
  cursor.close()
  conn.close()

def create_count_table(to_db: str, from_db: str):
  common_utils.log(f"create_count_table {to_db} {from_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`count`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`count` (
      `contig_id` INT NOT NULL COMMENT 'primary key of the `contig` table for the corresponding contig',
      `name` VARCHAR(50) NOT NULL COMMENT 'contig name',
      `gene_num` INT NOT NULL COMMENT 'Number of genes on this contig',
      `hit_num` INT NOT NULL COMMENT 'Number of contigs that have an MDS with this one',
      `mds_num` INT NOT NULL COMMENT 'Number of total MDSs of this contig with all hits',
      `pointer_num` INT NOT NULL COMMENT 'Number of total pointers of this contig with all hits',
      `ies_num` INT NOT NULL COMMENT 'Number of total strict of IESs of this contig with all hits',
      `properties_hit_num` INT NOT NULL COMMENT 'Number of contigs that have an entry in `properties` with this one',
      `non_gapped_num` INT NOT NULL COMMENT 'Number of `non_gapped` hits in the `properties` table for this contig',
      `non_overlapping_num` INT NOT NULL COMMENT 'Number of `non_overlapping` hits in the `properties` table for this contig',
      `non_repeating_num` INT NOT NULL COMMENT 'Number of `non_repeating` hits in the `properties` table for this contig',
      `exceeded_clique_limit_num` INT NOT NULL COMMENT 'Number of `exceeded_clique_limit` hits in the `properties` table for this contig',
      `weakly_complete_num` INT NOT NULL COMMENT 'Number of `weakly_complete` hits in the `properties` table for this contig',
      `strongly_complete_num` INT NOT NULL COMMENT 'Number of `strongly_complete` hits in the `properties` table for this contig',
      `weakly_consecutive_num` INT NOT NULL COMMENT 'Number of `weakly_consecutive` hits in the `properties` table for this contig',
      `strongly_consecutive_num` INT NOT NULL COMMENT 'Number of `strongly_consecutive` hits in the `properties` table for this contig',
      `weakly_ordered_num` INT NOT NULL COMMENT 'Number of `weakly_ordered` hits in the `properties` table for this contig',
      `strongly_ordered_num` INT NOT NULL COMMENT 'Number of `strongly_ordered` hits in the `properties` table for this contig',
      `weakly_non_scrambled_num` INT NOT NULL COMMENT 'Number of `weakly_non_scrambled` hits in the `properties` table for this contig',
      `strongly_non_scrambled_num` INT NOT NULL COMMENT 'Number of `strongly_non_scrambled` hits in the `properties` table for this contig',
      PRIMARY KEY (`contig_id`)
    ) COMMENT='Summary data about the hits on a contig';
    """
  )

  # initialize table to all 0
  cursor.execute(
    f"""
    INSERT INTO
      `{to_db}`.`count`
      (
        `contig_id`,
        `name`,
        `gene_num`,
        `hit_num`,
        `mds_num`,
        `pointer_num`,
        `ies_num`,
        `properties_hit_num`,
        `non_gapped_num`,
        `non_overlapping_num`,
        `non_repeating_num`,
        `exceeded_clique_limit_num`,
        `weakly_complete_num`,
        `strongly_complete_num`,
        `weakly_consecutive_num`,
        `strongly_consecutive_num`,
        `weakly_ordered_num`,
        `strongly_ordered_num`,
        `weakly_non_scrambled_num`,
        `strongly_non_scrambled_num`
      )
    SELECT
      `contig_id`,
      `name`,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0
    FROM
      `{to_db}`.`contig`;
    """
  )

  # get `gene`
  cursor.execute(
    f"""
    UPDATE
      `{to_db}`.`count` AS `C`,
      (
        SELECT
          `contig_id` AS `contig_id`,
          COUNT(*) AS `num`
        FROM `{to_db}`.`gene`
        WHERE `type` = 'gene'
        GROUP BY `contig_id`
      ) AS `G`
    SET `C`.`gene_num` = `G`.`num`
    WHERE `C`.`contig_id` = `G`.`contig_id`;
    """
  )

  for i in [0, 1]:
    # update mac and mic separately
    if i == 0:
      a_nuc = "mac"
      b_nuc = "mic"
    elif i == 1:
      a_nuc = "mic"
      b_nuc = "mic"
    else:
      raise Exception("Impossible.")

    # get `mds`, `ies_strict`, `pointer`
    for j in [0, 1, 2]:
      if j == 0:
        table = "match"
        field = "mds_num"
        where = "`is_preliminary` = 1"
      elif j == 1:
        table = "ies_strict"
        field = "ies_num"
        where = "1"
      elif j == 2:
        table = "pointer"
        field = "pointer_num"
        where = "`is_preliminary` = 1"
      else:
        raise Exception("Impossible.")

      cursor.execute(
        f"""
        UPDATE
          `{to_db}`.`count` AS `C`,
          (
            SELECT
              `{a_nuc}_contig_id` AS `contig_id`,
              COUNT(*) AS `num`
            FROM `{to_db}`.`{table}`
            WHERE {where}
            GROUP BY `{a_nuc}_contig_id`
          ) AS `T`
        SET `C`.`{field}` = `T`.`num`
        WHERE `C`.`contig_id` = `T`.`contig_id`;
        """
      )

    # get `contig_hits`
    cursor.execute(
      f"""
      UPDATE
        `{to_db}`.`count` AS `C`,
        (
          SELECT 
            `{a_nuc}_contig_id` AS `contig_id`,
            COUNT(DISTINCT `{b_nuc}_contig_id`) AS `num`
          FROM
            `{to_db}`.`match`
          WHERE
            `is_preliminary` = '1'
          GROUP BY
            `{a_nuc}_contig_id`
        ) AS `T`
      SET
        `C`.`hit_num` = `T`.`num`
      WHERE
        `C`.`contig_id` = `T`.`contig_id`;
      """
    )

    # update remaining properties
    cursor.execute(
      f"""
      UPDATE
        `{to_db}`.`count` AS `C`,
        (
          SELECT
            `{a_nuc}_contig_id` as `contig_id`,
            COUNT(*) AS `properties_hit_num`,
            SUM(`non_gapped`) AS `non_gapped_num`,
            SUM(`non_overlapping`) AS `non_overlapping_num`,
            SUM(`non_repeating`) AS `non_repeating_num`,
            SUM(`exceeded_clique_limit`) AS `exceeded_clique_limit_num`,
            SUM(`weakly_complete`) AS `weakly_complete_num`,
            SUM(`strongly_complete`) AS `strongly_complete_num`,
            SUM(`weakly_consecutive`) AS `weakly_consecutive_num`,
            SUM(`strongly_consecutive`) AS `strongly_consecutive_num`,
            SUM(`weakly_ordered`) AS `weakly_ordered_num`,
            SUM(`strongly_ordered`) AS `strongly_ordered_num`,
            SUM(`weakly_non_scrambled`) AS `weakly_non_scrambled_num`,
            SUM(`strongly_non_scrambled`) AS `strongly_non_scrambled_num`
          FROM
            `{to_db}`.`properties`
          GROUP BY
            `{a_nuc}_contig_id` 
        ) AS `P`
      SET
        `C`.`properties_hit_num` = `P`.`properties_hit_num`,
        `C`.`non_gapped_num` = `P`.`non_gapped_num`,
        `C`.`non_overlapping_num` = `P`.`non_overlapping_num`,
        `C`.`non_repeating_num` = `P`.`non_repeating_num`,
        `C`.`exceeded_clique_limit_num` = `P`.`exceeded_clique_limit_num`,
        `C`.`weakly_complete_num` = `P`.`weakly_complete_num`,
        `C`.`strongly_complete_num` = `P`.`strongly_complete_num`,
        `C`.`weakly_consecutive_num` = `P`.`weakly_consecutive_num`,
        `C`.`strongly_consecutive_num` = `P`.`strongly_consecutive_num`,
        `C`.`weakly_ordered_num` = `P`.`weakly_ordered_num`,
        `C`.`strongly_ordered_num` = `P`.`strongly_ordered_num`,
        `C`.`weakly_non_scrambled_num` = `P`.`weakly_non_scrambled_num`,
        `C`.`strongly_non_scrambled_num` = `P`.`strongly_non_scrambled_num`
      WHERE
        `C`.`contig_id` = `P`.`contig_id`;
      """
    )

  cursor.close()
  conn.close()

def create_alias_table(to_db: str):
  common_utils.log(f"create_alias_table {to_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`alias`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`alias` (
      `alias_id` INT NOT NULL AUTO_INCREMENT COMMENT 'primary key for the table',
      `id` INT NOT NULL COMMENT 'primary key of the table for the corresponding sequence',
      `name` VARCHAR(50) NOT NULL COMMENT 'the primary name of the sequence',
      `alias` VARCHAR(50) NOT NULL COMMENT 'alias of the sequence',
      `table` VARCHAR(50) NOT NULL COMMENT 'the table that `id` refers to',
      `type` VARCHAR(50) DEFAULT NULL COMMENT 'the type of alias',
      PRIMARY KEY (`alias_id`),
      KEY `alias` (`alias`)
    ) COMMENT='table of sequence aliases';
    """
  )

  cursor.close()
  conn.close()

def insert_alias_gene(to_db: str):
  common_utils.log(f"insert_alias_gene {to_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`alias`
    (
      `id`,
      `name`,
      `alias`,
      `table`,
      `type`
    )
    SELECT
      `G`.`gene_id` AS `id`,
      `G`.`attr_id` AS `name`,
      `G`.`attr_id` AS `alias`,
      'gene' AS `table`,
      'primary' AS `type`
    FROM `{to_db}`.`gene` AS `G`
    WHERE `G`.`type` = 'gene';
    """
  )

  cursor.close()
  conn.close()

def insert_alias_contig(to_db: str):
  common_utils.log(f"insert_alias_contig {to_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`alias`
    (
      `id`,
      `name`,
      `alias`,
      `table`,
      `type`
    )
    SELECT
      `C`.`contig_id` AS `id`,
      `C`.`name` AS `name`,
      `C`.`name` AS `name`,
      'contig' AS `table`,
      'primary' AS `type`
    FROM `{to_db}`.`contig` AS `C`;
    """
  )
  cursor.close()
  conn.close()


def insert_alias_file(to_db: str, file: str, table: str, nucleus: str):
  common_utils.log(f"insert_alias_file {to_db} {file} {table} {nucleus}")
  
  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`alias_temp`;");
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`alias_temp` (
      `name` VARCHAR(50),
      `alias` VARCHAR(50),
      `type` VARCHAR(50),
      KEY (`name`)
    );
    """
  )
  
  alias_data = common_utils.read_tsv(file)
  alias_data = alias_data.melt(
    id_vars = "primary",
    var_name = "type",
    value_name = "alias",
  )
  alias_data = alias_data.dropna(axis="index")
  alias_data = alias_data.rename({"primary": "name"}, axis="columns")
  mysql_utils.upload_in_chunks(
    alias_data,
    ["name", "alias", "type"],
    cursor,
    to_db,
    "alias_temp",
    constants.SQL_BATCH_UPLOAD_ROWS,
  )

  if table == 'contig':
    name_column = 'name'
    nucleus_column = 'nucleus'
  elif table == 'gene':
    name_column = 'attr_id'
    nucleus_column = 'contig_nucleus'
  else:
    raise Exception(f"Unexpected table: {table}")

  cursor.execute(
    f"""
    INSERT INTO
      `{to_db}`.`alias`
      (
        `id`,
        `name`,
        `alias`,
        `table`,
        `type`
      )
      SELECT
        `T`.`{table}_id`,
        `U`.`name`,
        `U`.`alias`,
        '{table}',
        `U`.`type`
      FROM `{to_db}`.`alias_temp` AS `U`
      INNER JOIN `{to_db}`.`{table}` AS `T`
      ON `T`.`{name_column}` = `U`.`name`
      AND `T`.`{nucleus_column}` = '{nucleus}';
    """
  )

  cursor.close()
  conn.close()


def create_variant_table(to_db: str):
  common_utils.log(f"create_variant_table {to_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`variant`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`variant` (
      `contig_id` INT NOT NULL AUTO_INCREMENT COMMENT 'primary key for the table',
      `name` VARCHAR(50) NOT NULL COMMENT 'the primary name of the sequence',
      `variant` VARCHAR(500) NOT NULL COMMENT 'variants/isoforms of the sequence',
      PRIMARY KEY (`contig_id`)
    ) COMMENT='Variants/isoforms of sequences';
    """
  )

def insert_variant_file(to_db: str, file: str):
  common_utils.log(f"insert_variant_file {to_db} {file}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`variant_temp`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`variant_temp` (
      `name` VARCHAR(50),
      `variant` VARCHAR(500),
      KEY (`name`)
    );
    """
  )
  
  variant_data = common_utils.read_tsv(file)
  mysql_utils.upload_in_chunks(
    variant_data,
    ["name", "variant"],
    cursor,
    to_db,
    "variant_temp",
    constants.SQL_BATCH_UPLOAD_ROWS,
  )
  
  cursor.execute(
    f"""
    INSERT INTO
      `{to_db}`.`variant`
      (
        `contig_id`,
        `name`,
        `variant`
      )
      SELECT
        `C`.`contig_id`,
        `C`.`name`,
        `V`.`variant`
      FROM `{to_db}`.`variant_temp` AS `V`
      INNER JOIN `{to_db}`.`contig` AS `C`
      ON `V`.`name` = `C`.`name`;
    """
  )

  cursor.close()
  conn.close()


def create_stats_table(to_db: str):
  common_utils.log(f"create_stats_table {to_db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`stats`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`stats` (
      `stats_id` INT NOT NULL AUTO_INCREMENT COMMENT 'Primary key for the table',
      `contig_mac_two_telomere` INT NOT NULL COMMENT 'Number of 2-telomere MAC contigs in the database',
      `contig_mac_one_telomere` INT NOT NULL COMMENT 'Number of 1-telomere MAC contigs in the database',
      `contig_mac_zero_telomere` INT NOT NULL COMMENT 'Number of 0-telomere MAC contigs in the database',
      `contig_mac_total` INT NOT NULL COMMENT 'Number of total MAC contigs in the database',
      `contig_mic_two_telomere` INT NOT NULL COMMENT 'Number of 2-telomere MIC contigs in the database',
      `contig_mic_one_telomere` INT NOT NULL COMMENT 'Number of 1-telomere MIC contigs in the database',
      `contig_mic_zero_telomere` INT NOT NULL COMMENT 'Number of 0-telomere MIC contigs in the database',
      `contig_mic_total` INT NOT NULL COMMENT 'Number of total MIC contigs in the database',
      `gene_mac_total` INT NOT NULL COMMENT 'Number of total MAC genes in the database',
      `gene_mic_total` INT NOT NULL COMMENT 'Number of total MIC genes in the database',
      `mds_mac_total` INT NOT NULL COMMENT 'Number of total MAC MDSs in the database',
      `mds_mic_total` INT NOT NULL COMMENT 'Number of total MIC MDSs in the database',
      `ies_mac_total` INT NOT NULL COMMENT 'Number of total MAC IESs (strict) in the database',
      `ies_mic_total` INT NOT NULL COMMENT 'Number of total MIC IESs (strict) in the database',
      `pointer_mac_total` INT NOT NULL COMMENT 'Number of total MAC pointers in the database',
      `pointer_mic_total` INT NOT NULL COMMENT 'Number of total MIC pointers in the database',
      `properties_non_gapped` INT NOT NULL COMMENT 'Number of non-gapped rearrangement maps in the database',
      `properties_non_overlapping` INT NOT NULL COMMENT 'Number of non-overlapping rearrangement maps in the database',
      `properties_non_repeating` INT NOT NULL COMMENT 'Number of non-repeating rearrangement maps in the database',
      `properties_exceeded_clique_limit` INT NOT NULL COMMENT 'Number of exceeded-clique-limit rearrangement maps in the database',
      `properties_weakly_complete` INT NOT NULL COMMENT 'Number of weakly-complete rearrangement maps in the database',
      `properties_strongly_complete` INT NOT NULL COMMENT 'Number of strongly-complete rearrangement maps in the database',
      `properties_weakly_consecutive` INT NOT NULL COMMENT 'Number of weakly-consecutive rearrangement maps in the database',
      `properties_strongly_consecutive` INT NOT NULL COMMENT 'Number of strongly-consecutive rearrangement maps in the database',
      `properties_weakly_ordered` INT NOT NULL COMMENT 'Number of weakly-ordered rearrangement maps in the database',
      `properties_strongly_ordered` INT NOT NULL COMMENT 'Number of strongly-ordered rearrangement maps in the database',
      `properties_weakly_non_scrambled` INT NOT NULL COMMENT 'Number of weakly-non-scrambled rearrangement maps in the database',
      `properties_strongly_non_scrambled` INT NOT NULL COMMENT 'Number of strongly-non-scrambled rearrangement maps in the database',
      `properties_total` INT NOT NULL COMMENT 'Number of total rearrangement maps in the database',
      PRIMARY KEY (`stats_id`)
    ) COMMENT='Summary statistics of the database';
    """
  )

  stats = {}

  for nucleus in ['mac', 'mic']:
    # contig query
    cursor.execute(
      f"""
      SELECT
        SUM( ( ( `telomere_five_start` > 0 ) + ( `telomere_three_start` > 0 ) ) = 2 ) AS `two_telomere`,
        SUM( ( ( `telomere_five_start` > 0 ) + ( `telomere_three_start` > 0 ) ) = 1 ) AS `one_telomere`,
        SUM( ( ( `telomere_five_start` > 0 ) + ( `telomere_three_start` > 0 ) ) = 0 ) AS `zero_telomere`,
        COUNT(*) AS `total`
      FROM `{to_db}`.`contig`
      WHERE `nucleus` = '{nucleus}';
      """
    )
    result = cursor.fetchall()[0]
    stats[f"contig_{nucleus}_two_telomere"] = result["two_telomere"]
    stats[f"contig_{nucleus}_one_telomere"] = result["one_telomere"]
    stats[f"contig_{nucleus}_zero_telomere"] = result["zero_telomere"]
    stats[f"contig_{nucleus}_total"] = result["total"]

    # gene query
    cursor.execute(
      f"""
      SELECT
        COUNT(*) AS `total`
      FROM `{to_db}`.`gene` AS `G`
      INNER JOIN `{to_db}`.`contig` AS `C`
        ON `C`.`contig_id` = `G`.`contig_id`
      WHERE `G`.`type` = 'gene'
      AND `C`.`nucleus` = '{nucleus}';
      """
    )
    result = cursor.fetchall()[0]
    stats[f"gene_{nucleus}_total"] = result["total"]

    # MDS query
    cursor.execute(
      f"""
      SELECT
        COUNT(*) AS `total`
      FROM `{to_db}`.`match`
        WHERE `is_preliminary` = 1;
      """
    )
    result = cursor.fetchall()[0]
    stats[f"mds_{nucleus}_total"] = result["total"]

    # IES strict query
    if nucleus == "mic":
      cursor.execute(
        f"""
        SELECT
          COUNT(*) AS `total`
        FROM `{to_db}`.`ies_strict`;
        """
      )
      result = cursor.fetchall()[0]
      stats[f"ies_{nucleus}_total"] = result["total"]
    elif nucleus == "mac":
      stats[f"ies_{nucleus}_total"] = 0
    else:
      raise Exception("Impossible.")

    # pointer query
    cursor.execute(
      f"""
      SELECT
        COUNT(*) AS `total`
      FROM `{to_db}`.`pointer`;
      """
    )
    result = cursor.fetchall()[0]
    stats[f"pointer_{nucleus}_total"] = {
      "mic": 2 * result["total"], # count each pointer 2x for MIC
      "mac": result["total"],
    }[nucleus]

  # properties query
  cursor.execute(
    f"""
    SELECT
      SUM(`non_gapped`) AS `non_gapped`,
      SUM(`non_overlapping`) AS `non_overlapping`,
      SUM(`non_repeating`) AS `non_repeating`,
      SUM(`exceeded_clique_limit`) AS `exceeded_clique_limit`,
      SUM(`weakly_complete`) AS `weakly_complete`,
      SUM(`strongly_complete`) AS `strongly_complete`,
      SUM(`weakly_consecutive`) AS `weakly_consecutive`,
      SUM(`strongly_consecutive`) AS `strongly_consecutive`,
      SUM(`weakly_ordered`) AS `weakly_ordered`,
      SUM(`strongly_ordered`) AS `strongly_ordered`,
      SUM(`weakly_non_scrambled`) AS `weakly_non_scrambled`,
      SUM(`strongly_non_scrambled`) AS `strongly_non_scrambled`,
      COUNT(*) AS `total`
    FROM `{to_db}`.`properties`;
    """
  )
  result = cursor.fetchall()[0]
  stats['properties_non_gapped'] = result['non_gapped']
  stats['properties_non_overlapping'] = result['non_overlapping']
  stats['properties_non_repeating'] = result['non_repeating']
  stats['properties_exceeded_clique_limit'] = result['exceeded_clique_limit']
  stats['properties_weakly_complete'] = result['weakly_complete']
  stats['properties_strongly_complete'] = result['strongly_complete']
  stats['properties_weakly_consecutive'] = result['weakly_consecutive']
  stats['properties_strongly_consecutive'] = result['strongly_consecutive']
  stats['properties_weakly_ordered'] = result['weakly_ordered']
  stats['properties_strongly_ordered'] = result['strongly_ordered']
  stats['properties_weakly_non_scrambled'] = result['weakly_non_scrambled']
  stats['properties_strongly_non_scrambled'] = result['strongly_non_scrambled']
  stats['properties_total'] = result['total']

  fields = "(" + ','.join([f"`{x}`" for x in stats.keys()]) + ")"
  values = "(" + ','.join([f"'{x}'" for x in stats.values()]) + ")"

  # insert values
  cursor.execute(
    f"""
    INSERT INTO `{to_db}`.`stats` {fields}
    VALUES {values};
    """
  )


def create_protein_table(db: str):
  common_utils.log(f"create_protein_table {db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{db}`.`protein`;")
  cursor.execute(
    f"""
    CREATE TABLE `{db}`.`protein` (
      `prot_id` int(11) NOT NULL AUTO_INCREMENT COMMENT 'Index of the `protein` table',
      `gene_id` int(11) DEFAULT NULL COMMENT 'Index of the corresponding mRNA feature in the table `gene`',
      `attr_id` varchar(50) DEFAULT NULL COMMENT 'Name of the corresponding mRNA feature',
      `length` int(11) NOT NULL COMMENT 'Length of the protein sequence in amino acids',
      `sequence` longtext COMMENT 'Sequence of amino acids for the protein',
      PRIMARY KEY (`prot_id`),
      UNIQUE KEY `gene_id` (`gene_id`),
      KEY `attr_id` (`attr_id`)
    ) COMMENT='Table of Predicted Proteins';
    """
  )

  cursor.close()
  conn.close()

def insert_protein_file(db: str, file: str):
  common_utils.log(f"insert_protein_file {db} {file}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(f"DROP TABLE IF EXISTS `{db}`.`protein_temp`;")
  cursor.execute(
    f"""
    CREATE TABLE `{db}`.`protein_temp` (
      `attr_id` VARCHAR(50),
      `sequence` longtext,
      KEY (`attr_id`)
    );
    """
  )

  data = common_utils.read_tsv(file)
  mysql_utils.upload_in_chunks(
    data,
    ["attr_id", "sequence"],
    cursor,
    db,
    "protein_temp",
    constants.SQL_BATCH_UPLOAD_ROWS,
  )

  cursor.execute(
    f"""
    INSERT INTO
    `{db}`.`protein`
    (
      `gene_id`,
      `attr_id`,
      `length`,
      `sequence`
    )
    SELECT
      `G`.`gene_id`,
      `G`.`attr_id`,
      LENGTH(`T`.`sequence`),
      `T`.`sequence`
    FROM `{db}`.`protein_temp` AS `T`
    INNER JOIN `{db}`.`gene` AS `G`
    ON `G`.`attr_id` = `T`.`attr_id`;
    """
  )

  cursor.close()
  conn.close()

def drop_temp_tables(db: str):
  common_utils.log(f"drop_temp_tables {db}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)

  temp_tables = [
    "name_temp",
    "alias_temp",
    "variant_temp",
    "protein_temp",
  ]
  for table in temp_tables:
    cursor.execute(f"DROP TABLE IF EXISTS `{db}`.`{table}`;")

def add_to_directory(
  db: str,
  name: str,
  description: str,
  organism: str,
  download_dir: str,
  assembly: str,
  url: str,
):
  common_utils.log(
    "add_to_directory" + " " +
    db + " " +
    name + " " +
    description + " " +
    organism + " " +
    download_dir + " " +
    assembly + " " +
    url
  )

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)

  cursor.execute(
    f"""
    INSERT INTO `mds_ies_db`.`db`
    (
      `name`,
      `description`,
      `sql_name`,
      `organism`,
      `download_dir`,
      `assembly`,
      `url`
    )
    VALUES
    (
      '{name}',
      '{description}',
      '{db}',
      '{organism}',
      '{download_dir}',
      '{assembly}',
      '{url}'
    );
    """
  )

  cursor.close()
  conn.close()

def dump_table(
  db: str,
  table: str,
  file: str,
):
  common_utils.log(f"dump_table {db} {table} {file}");

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  cursor.execute(
    f"""
    SELECT `COLUMN_NAME`
    FROM `INFORMATION_SCHEMA`.`COLUMNS`
    WHERE `TABLE_SCHEMA` = '{db}'
    AND `TABLE_NAME` = '{table}';
    """
  )
  columns = cursor.fetchall()
  columns = [x[0] for x in columns]

  cursor.execute(f"SELECT COUNT(*) FROM `{db}`.`{table}`;")
  num_rows = cursor.fetchall()[0][0]

  data = mysql_utils.download_in_chunks(
    db,
    table,
    columns,
    cursor,
    num_rows,
    batch_size = constants.SQL_BATCH_DOWNLOAD_ROWS,
  )

  common_utils.write_tsv(data, file)


# # #track name="prec_eliminated_sequences" description="intervals comprising complement of the precursor segments labelled comp_[left-flanking-segment-prod-id]_[left-flanking-segment-index]_[right-flanking-segment-prod-id]_[right-flanking-segment-index]" itemRgb-"On"
# # OXYTRI_MIC_33550	1	8703	none_0_Contig4386.0_-1	0	+	1	8703	255,0,153
# # sdrap_oxy_mac2012_100720_prec_eliminated_sequences.bed
# def make_ie_bed(to_db: str, from_db: str, ies_type: str):
#   handle = open(f"../output/{from_db}_prec_{ies_type}_ies_sequences.bed", "w")

#   # write header
#   if ies_type == "strict":
#     handle.write(
#       '#track name="prec_strict_ies_sequences" ' +
#       'description="Strict IESs. Intervals on the precursor segment between two consecutive matches of a product segment ' +
#       'that do not overlap matches of any other product segment. ' +
#       'Labelled [product-name]_[left-flanking-match-indexes]_[right-flanking-match-indexes]" ' +
#       'itemRgb-"On"\n'
#     )
#   elif ies_type == "weak":
#     handle.write(
#       '#track name="prec_weak_ies_sequences" ' +
#       'description="Weak IESs. Intervals on the precursor segment between two consecutive matches of a product segment.' +
#       'Labelled [product-name]_[left-flanking-match-indexes]_[right-flanking-match-indexes]" ' +
#       'itemRgb-"On"\n'
#     )
#   else:
#     raise Exception("Unknown ies type: " + str(ies_type))

#   conn = mysql_utils.get_connection()
#   cursor = conn.cursor(dictionary=True)
#   cursor.execute("SELECT COUNT(*) AS `count` FROM `{to_db}`.`ies_{$type}`;")
#   num_ies = cursor.fetchall()["count"]
  
#   batch_size_rows = 10000
#   for i in range(0, num_ies, batch_size_rows):
#     cursor.execute(
#       f"""
#       SELECT
#         `mic_name`,
#         `mic_start`,
#         `mic_end`,
#         `mac_name`,
#         `left_index`,
#         `right_index`
#       FROM `{to_db}`.`ies_{ies_type}`
#       LIMIT {i}, {batch_size_rows};
#       """
#     )

#     for result in cursor.fetchall():
#       handle.write(
#         f"{result['mic_name']}\t{result['mic_start']}\t{result['mic_end']}\t" + # chrom chromStart chromEnd
#         f"{result['mac_name']}_{result['left_index']}_{result['right_index']}\t" + # name
#         f"0\t+\t{result['mic_start']}\t{result['mic_end']}\t255,0,153\n" # score strand thickStart thickEnd itemRGB
#       )

#   cursor.close()
#   conn.close()

from_db = "sdrap_oxy_mac2012_May_30_2022"
to_db = "hello_world"
# conn = mysql_utils.get_connection()

# Retest everything
# create_name_temp_table("hello_world", from_db, "Contig.*", "OXYTRI.*")
# create_alias_table(to_db, from_db, ["main\\oxy_tri_jrb310_mac_2012_alias.csv", "main\\oxy_tri_jrb310_mic_2014_alias.csv"])
# create_contig_table(to_db, from_db)
# create_match_table(to_db, from_db)
# create_pointer_table(to_db, from_db)
# create_properties_table(to_db, from_db)
# create_parameter_table(to_db, from_db)
# create_coverage_table(to_db, from_db)
# GENES AND IES NEED TO BE DONE HERE
# create_count_table(to_db, from_db)
# create_alias_table(to_db, from_db)
# insert_alias_contig(to_db)
# insert_alias_gene(to_db)
# insert_alias_file(
#   to_db,
#   "data/oxy_tri_jrb310_mac_2012_alias.tsv",
#   "contig",
#   "mac",
# )
# insert_alias_file(
#   to_db,
#   "data/oxy_tri_jrb310_mic_2014_alias.tsv",
#   "contig",
#   "mic",
# )
# create_variant_table(to_db)
# insert_variant_file(to_db, "data/oxy_tri_jrb310_mac_2020_variant.tsv")
# create_stats_table(to_db)
create_protein_table(to_db)
insert_protein_file(to_db, "data/Oxytricha_trifallax_022112_aa.tsv")
# add_to_directory(
#   "hello_world",
#   "spacemonkey",
#   "Super duper space monkey cowboy organism.",
#   "Only the finest specimen of space monkey dna.",
#   "hello_world",
#   "my_assembly",
#   "hello_world"
# )
# dump_table("hello_world", "variant", "variant.tsv")