import mysql.connector
import pandas as pd
import numpy as np
import collections

import common_utils
import interval_utils
import mysql_utils

SQL_BATCH_UPLOAD_ROWS = 1000

def create_contig_table(
  to_db : str,
  from_db : str,
):
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
    LEFT JOIN (
      SELECT `name`, `id`
      FROM `{to_db}`.`alias`
      WHERE `table` = 'contig'
      GROUP BY `id`
    ) AS `A`
      ON `N`.`nuc_id` = `A`.`id`
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
    for i in range(0, alias_to.shape[0], SQL_BATCH_UPLOAD_ROWS):
      common_utils.log(f"{i} / {alias_to.shape[0]}")
      values = mysql_utils.make_sql_values(
        alias_to.loc[
          i : (i + SQL_BATCH_UPLOAD_ROWS - 1),
          ["id", "name", "alias", "table", "type"]
        ])
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
        VALUES {values};
        """
      )
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

from_db = "sdrap_oxy_mac2012_May_30_2022"
to_db = "hello_world"
# conn = mysql_utils.get_connection()

# Retest everything
# create_alias_table(to_db, from_db, ["main\\oxy_tri_jrb310_mac_2012_alias.csv", "main\\oxy_tri_jrb310_mic_2014_alias.csv"])
# create_contig_table(to_db, from_db)
# create_match_table(to_db, from_db)
# create_pointer_table(to_db, from_db)
# create_properties_table(to_db, from_db)
# create_parameter_table(to_db, from_db)
create_coverage_table(to_db, from_db)

# cursor = conn.cursor(dictionary=True)
# # cursor.execute(f"SELECT MAX(`mac_contig_id`) FROM `hello_world`.`properties`")
# cursor.execute(f"""
#   SELECT
#     `contig_id`,
#     `non_tel_length`
#   FROM `{to_db}`.`contig`
#   WHERE `contig_id` BETWEEN {0} AND {1000}
# """)
# cursor.execute(f"""
#   SELECT
#     `contig_id`,
#     `non_tel_length`
#   FROM `{to_db}`.`contig`
#   WHERE `contig_id` BETWEEN {0} AND {1000}
# """)
# x = cursor.fetchall()
# cursor.execute(f"""
#   INSERT INTO `{to_db}`.`tab`
#     SELECT 123 AS `x` UNION ALL SELECT 345 AS `x`;
# """)

# conn.close()