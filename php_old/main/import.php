<?php
include "../databaseConnect.php";
include "../sqlUtil.php";
include "util.php";

function checkFromDatabaseExists($link, $from) {
  printf("Starting: verify `{$from}` exists;\n");
  tryQuery($link, "SELECT * FROM `{$from}`.`alias` WHERE 1=0;");
  printf("Success: verified `{$from}` exists\n");
}

function createDatabase($link, $db) {
  tryQuery($link, "DROP DATABASE IF EXISTS `{$db}`;");
  tryQuery($link, "CREATE DATABASE `{$db}`;");
  printf("Success: created database `{$db}`\n");
}

function createTempName($link, $to) {
  printf("Starting: `{$to}`.`temp_name`\n");
  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`temp_name`;");
  tryQuery(
    $link, 
    "CREATE TABLE `{$to}`.`temp_name` (
      `contig_id` INT NOT NULL,
      `name` VARCHAR(50) NOT NULL,
      PRIMARY KEY (`contig_id`),
      KEY `name` (`name`)
    );"
  );
  printf("Success: created `{$to}`.`temp_name`\n");
}

function insertTempName($link, $to, $from, $macRegex, $micRegex) {
  printf("Starting: inserting `{$to}`.`temp_name`\n");
  tryQuery(
    $link,
    "INSERT INTO `{$to}`.`temp_name`
    (
      `contig_id`,
      `name`
    )
    SELECT
      `nuc_id`,
      `alias`
    FROM `{$from}`.`alias`
    WHERE `alias` REGEXP '{$macRegex}'
    OR `alias` REGEXP '{$micRegex}';"
  );
  printf("Success: inserted `{$to}`.`temp_name`\n");
}

function insertTempName_oxyTriJrb310_mac2020_genbankV1ToOxydb($link, $to, $from, $micRegex) {
  printf("Starting: inserting genbank.1 to oxydb `{$to}`.`temp_name`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`temp_alias`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`temp_alias` (
      `genbank.1` VARCHAR(50),
      `oxydb` VARCHAR(50),
      KEY (`genbank.1`),
      KEY (`oxydb`)
    ) COMMENT='Summary data about the hits on a contig';"
  );
  
  printf("Success: created `{$to}`.`temp_alias`\n");

  $allRows = readCSV('oxy_tri_jrb310_mac_2020_genbankV1_to_oxydb.csv');
  $step = 10000;
  for ($i = 1; $i < count($allRows); $i += $step) {
    $values = array();
    for ($j = $i; $j < min($i + $step, count($allRows)); $j++) {
      $values[] = "('{$allRows[$j][0]}','{$allRows[$j][1]}')";
    }
    $values = implode(',', $values);
    tryQuery($link, "INSERT INTO `{$to}`.`temp_alias` VALUES {$values};");
  }

  tryQuery(
    $link,
    "INSERT INTO `{$to}`.`temp_name`
    (
      `contig_id`,
      `name`
    )
    SELECT
      `F`.`nuc_id`,
      `T`.`oxydb`
    FROM `{$to}`.`temp_alias` AS `T`
    INNER JOIN `{$from}`.`alias` AS `F`
      ON `F`.`alias` = `T`.`genbank.1`;"
  );

  tryQuery(
    $link,
    "INSERT INTO `{$to}`.`temp_name`
    (
      `contig_id`,
      `name`
    )
    SELECT
      `nuc_id`,
      `alias`
    FROM `{$from}`.`alias`
    WHERE `alias` REGEXP '{$micRegex}';"
  );
  printf("Success: inserted `{$to}`.`temp_name`\n");
}

function createContig($link, $to) {
  printf("Starting: create `{$to}`.`contig`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`contig`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`contig` (
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
    );"
  );
  printf("Success: created `{$to}`.`contig`\n");
}

function insertContig($link, $to, $from) {
  printf("Starting: insert `{$to}`.`contig`\n");

  tryQuery($link,
    "INSERT INTO `{$to}`.`contig`
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
    FROM `{$from}`.`nucleotide` AS `N`
    LEFT JOIN `{$to}`.`temp_name` AS `A`
      ON `N`.`nuc_id` = `A`.`contig_id`
    LEFT JOIN `{$from}`.`feature` AS `F`
      ON `N`.`feat_id` = `F`.`feat_id`
    LEFT JOIN `{$from}`.`telomere` AS `T`
      ON `N`.`nuc_id` = `T`.`nuc_id`;"
  );

  printf("Success: inserted `{$to}`.`contig`\n");
}

function createHsp($link, $to) {
  printf("Starting: create `{$to}`.`hsp`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`hsp`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`hsp` (
      `hsp_id` INT NOT NULL,
      `mac_name` VARCHAR(50) NOT NULL COMMENT 'Name of MAC contig',
      `mac_contig_id` INT NOT NULL COMMENT 'Index of the MAC contig in the table `contig`',
      `mac_start` INT NOT NULL COMMENT 'Starting position of the MDS in the MAC contig',
      `mac_end` INT NOT NULL COMMENT 'Ending position of the MDS in the MAC contig',
      `mic_name` VARCHAR(50) NOT NULL COMMENT 'Name of MIC contig',
      `mic_contig_id` INT NOT NULL COMMENT 'Index of the MIC contig in the table `contig`',
      `mic_start` INT NOT NULL COMMENT 'Starting position of the MDS in the MAC contig',
      `mic_end` INT NOT NULL COMMENT 'Ending position of the MDS in the MAC contig',
      `mic_direct` SET('+','-','*') NOT NULL,
      `length` INT NOT NULL COMMENT 'Length of the MDS in the MAC',
      `pident` FLOAT NOT NULL COMMENT 'Percentage of identical matches',
      `mismatch` INT NOT NULL COMMENT 'Number of mismatches',
      `evalue` VARCHAR(10) NOT NULL COMMENT 'Expected value of seeing as good of a match at random',
      `bitscore` FLOAT NOT NULL COMMENT 'High scoring pair bit score',
      PRIMARY KEY (`hsp_id`),
      KEY `mac_contig_id` (`mac_contig_id`),
      KEY `mic_contig_id` (`mic_contig_id`)
    );"
  );

  printf("Success: created `{$to}`.`hsp`\n");
}

// FIXME: Should try and do this in batches of say 100,000 ?
function insertHsp($link, $to, $from) {
  printf("Starting: insert `{$to}`.`hsp`\n");

  $max = mysqli_fetch_row( tryQuery($link, "SELECT MAX(`hsp_id`) FROM `{$from}`.`hsp`;") );
  $max = $max[0];

  $step = 10000; // do 10000 at a time
  for ($start = 0; $start <= $max; $start += $step) {
    $end = $start + $step - 1;
    printf("Starting: insert into `{$to}`.`hsp` rows {$start}..{$end}\n");
    flush();
    tryQuery($link,
      "INSERT INTO `{$to}`.`hsp`
      (
        `hsp_id`,
        `mac_name`,
        `mac_contig_id`,
        `mac_start`,
        `mac_end`,
        `mic_name`,
        `mic_contig_id`,
        `mic_start`,
        `mic_end`,
        `mic_direct`,
        `length`,
        `pident`,
        `mismatch`,
        `evalue`,
        `bitscore`
      )
      SELECT
        `H`.`hsp_id`,
        `NMAC`.`name`,
        `NMAC`.`nuc_id`,
        `H`.`prod_start`,
        `H`.`prod_end`,
        `NMIC`.`name`,
        `NMIC`.`nuc_id`,
        `H`.`prec_start`,
        `H`.`prec_end`,
        `H`.`orientation`,
        `H`.`length`,
        `H`.`pident`,
        `H`.`mismatch`,
        `H`.`evalue`,
        `H`.`bitscore`
      FROM `{$from}`.`hsp` AS `H`
      LEFT JOIN `{$to}`.`contig` AS `NMIC`
        ON `H`.`prec_nuc_id` = `NMIC`.`contig_id`
      LEFT JOIN `{$to}`.`contig` AS `NMAC`
        ON `H`.`prod_nuc_id` = `NMAC`.`contig_id`
      WHERE
        {$start} <= `H`.`hsp_id`
      AND
        `H`.`hsp_id` <= {$end};"
    );
  }

  printf("Success: inserted `{$to}`.`hsp`\n");
}

function createMatch($link, $to) {
  printf("Starting: create `{$to}`.`match`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`match`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`match` (
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
      `hsp_id` TEXT CHARACTER SET latin1 NOT NULL COMMENT 'string representing a list of _ delimited primary keys of the entries in the `hsp` table that merged to become the match',
      `index` INT NOT NULL COMMENT 'index of the match or fragment in the arrangement of the corresponding product on the corresponding precursor',
      `pre_cov` DECIMAL(8,2) NOT NULL COMMENT 'percentage of product segment of preliminary match covered by product segment of additional match',
      `add_cov` DECIMAL(8,2) NOT NULL COMMENT 'percentage of product segment of additional match covered by product segment of preliminary match',
      `is_preliminary` BOOLEAN NOT NULL COMMENT 'indicates whether the match is a preliminary match (1 means yes/true, 0 means no/false)',
      `is_additional` BOOLEAN NOT NULL COMMENT 'indicates whether the match is an additional match (1 means yes/true, 0 means no/false)',
      `is_fragment` BOOLEAN NOT NULL COMMENT 'indicates whether the match is a fragment (1 means yes/true, 0 means no/false)',
      PRIMARY KEY (`match_id`),
      KEY `mac_contig_id` (`mac_contig_id`),
      KEY `mic_contig_id` (`mic_contig_id`)
    );"
  );
  printf("Success: created `{$to}`.`match`\n");
}

function insertMatch($link, $to, $from) {
  printf("Starting: insert `{$to}`.`match`\n");

  tryQuery($link,
    "INSERT INTO `{$to}`.`match`
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
      `hsp_id`,
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
      `M`.`hsp_id`,
      ABS(`M`.`index`),
      `M`.`pre_cov`,
      `M`.`add_cov`,
      `M`.`is_preliminary`,
      `M`.`is_additional`,
      `M`.`is_fragment`
    FROM `{$from}`.`match` AS `M`
    LEFT JOIN `{$to}`.`contig` AS `NMIC`
      ON `M`.`prec_nuc_id` = `NMIC`.`contig_id`
    LEFT JOIN `{$to}`.`contig` AS `NMAC`
      ON `M`.`prod_nuc_id` = `NMAC`.`contig_id`;"
  );
  printf("Success: inserted `{$to}`.`match`\n");
}

function createPointer($link, $to) {
  printf("Starting: create `{$to}`.`pointer`\n");
  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`pointer`;");
  tryQuery($link, 
    "CREATE TABLE `{$to}`.`pointer` (
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
    );"
  );
  printf("Success: created `{$to}`.`pointer`\n");
}

function insertPointer($link, $to, $from) {
  printf("Starting: insert `{$to}`.`pointer`\n");
  tryQuery($link,
    "INSERT INTO `{$to}`.`pointer`
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
    FROM `{$from}`.`pointer` AS `P`
    LEFT JOIN `{$to}`.`contig` AS `NMIC`
      ON `P`.`prec_nuc_id` = `NMIC`.`contig_id`
    LEFT JOIN `{$to}`.`contig` AS `NMAC`
      ON `P`.`prod_nuc_id` = `NMAC`.`contig_id`
    LEFT JOIN `{$to}`.`match` AS `MLEFT`
      ON `P`.`left_match_id` = `MLEFT`.`match_id`
    LEFT JOIN `{$to}`.`match` AS `MRIGHT`
      ON `P`.`right_match_id` = `MRIGHT`.`match_id`;"
  );
  printf("Success: inserted `{$to}`.`pointer`\n");
}

function createProperties($link, $to) {
  printf("Starting: create `{$to}`.`properties`\n");
  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`properties`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`properties` (
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
    ) COMMENT='table of properties for each arrangement with coverage above the threshold found in the `parameter` table in the `property_min_coverage` column';"
  );
  printf("Success: created `{$to}`.`properties`\n");
}

function insertProperties($link, $to, $from) {
  printf("Starting: insert `{$to}`.`properties`\n");
  tryQuery(
    $link, 
    "INSERT INTO `{$to}`.`properties`
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
    FROM `{$from}`.`properties` AS `P`
    LEFT JOIN `{$to}`.`contig` AS `NMIC`
      ON `P`.`prec_nuc_id` = `NMIC`.`contig_id`
    LEFT JOIN `{$to}`.`contig` AS `NMAC`
      ON `P`.`prod_nuc_id` = `NMAC`.`contig_id`;"
  );
  printf("Success: inserted `{$to}`.`properties`\n");
}

function createParameter($link, $to) {
  printf("Starting: create `{$to}`.`parameter`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`parameter`;");
  tryQuery($link, 
    "CREATE TABLE `{$to}`.`parameter` (
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
    ) COMMENT='table of user defined parameters used when running the pipeline';"
  );

  printf("Success: created `{$to}`.`parameter`\n");
}

function insertParameter($link, $to, $from) {
  tryQuery($link, "INSERT INTO `{$to}`.`parameter` SELECT * from `{$from}`.`parameter`;");
  printf("Success: inserted `{$to}`.`parameter`\n");
}

function createCoverage($link, $to) {
  printf("Starting: create `{$to}`.`coverage`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`coverage`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`coverage` (
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
    ) COMMENT='Summary data about the matches between a pair of contig.';"
  );

  printf("Success: created `{$to}`.`coverage`\n");
}

// update properties table with the pairwise base pair count and mac coverage values
// FIXME: Remove this and just stick with the properties table. If you need more properties
// run another SDRAP annotation and set the threshold higher.
function insertCoverage($link, $to) {
  printf("Starting: insert into `{$to}`.`coverage`\n");

  function overlaps($a, $b) {
    if ($a['start'] > $b['end']+1) return false;
    if ($b['start'] > $a['end']+1) return false;
    return true;
  }
  function unionIntervals($intervals) {
    usort($intervals, function($x, $y) { return $x['start'] - $y['start']; });
    $union = array();
    $last = null;
    foreach ($intervals as $x) {
      if ($last != null && overlaps($last, $x)) {
        $last['end'] = max($last['end'], $x['end']);
        $union[count($union)-1] = $last;
      } else {
        $union[] = $x;
        $last = $x;
      }
    }
    return $union;
  }
  $maxMacId = mysqli_fetch_row(tryQuery($link, "SELECT MAX(`mac_contig_id`) FROM `{$to}`.`properties`"));
  $maxMacId = $maxMacId[0];
  $step = 100;
  $covId = 1;
  for ($startMacId = 0; $startMacId <= $maxMacId; $startMacId += $step) {
    $endMacId = min($startMacId + $step - 1, $maxMacId);

    printf("Starting: insert `{$to}`.`coverage` MAC id's {$startMacId}..{$endMacId}\n");

    $telomereQuery = tryQuery($link,
      "SELECT
        `contig_id`,
        `non_tel_length`
      FROM
        `{$to}`.`contig`
      WHERE
        {$startMacId} <= `contig_id`
      AND
        `contig_id` <= {$endMacId};"
    );


    $macLength = array();
    while ($assoc = mysqli_fetch_assoc($telomereQuery)) {
      $macLength[$assoc['contig_id']] = intval($assoc['non_tel_length']);
    }

    // get the MDS's
    $matchQuery = tryQuery($link,
      "SELECT
        `mac_contig_id`,
        `mic_contig_id`,
        `mac_start`,
        `mac_end`,
        `mic_start`,
        `mic_end`
      FROM
        `{$to}`.`match`
      WHERE
        {$startMacId} <= `mac_contig_id`
      AND
        `mac_contig_id` <= {$endMacId}
      AND
        `is_preliminary` = 1;"
    );

    $allMacStart = array();
    $allMacEnd = array();
    $allMicStart = array();
    $allMicEnd = array();
    $groups = array();
    while ($assoc = mysqli_fetch_assoc($matchQuery)) {
      $assoc['mac_contig_id'] = intval($assoc['mac_contig_id']);
      $assoc['mic_contig_id'] = intval($assoc['mic_contig_id']);
      $assoc['mac_start'] = intval($assoc['mac_start']);
      $assoc['mac_end'] = intval($assoc['mac_end']);
      $assoc['mic_start'] = intval($assoc['mic_start']);
      $assoc['mic_end'] = intval($assoc['mic_end']);
      $key = "{$assoc['mac_contig_id']}_{$assoc['mic_contig_id']}";
      $allMacStart[$key] = isset($allMacStart[$key]) ? min($allMacStart[$key], $assoc['mac_start']) : $assoc['mac_start'];
      $allMacEnd[$key] = isset($allMacEnd[$key]) ? max($allMacEnd[$key], $assoc['mac_end']) : $assoc['mac_end'];
      $allMicStart[$key] = isset($allMicStart[$key]) ? min($allMicStart[$key], $assoc['mic_start']) : $assoc['mic_start'];
      $allMicEnd[$key] = isset($allMicEnd[$key]) ? max($allMicEnd[$key], $assoc['mic_end']) : $assoc['mic_end'];
      $groups[$key][] = array(
        'start' => $assoc['mac_start'],
        'end' => $assoc['mac_end']
      );
    }

    $tableRows = '';
    foreach ($groups as $key => $intervals) {
      $union = unionIntervals($intervals);
      $coverageLength = 0;
      foreach ($union as $x) {
        $coverageLength += $x['end'] - $x['start'] + 1;
      }
      
      $macMicId = explode('_', $key);
      $macId = intval($macMicId[0]);
      $micId = intval($macMicId[1]);
      $coveragePercent = 100 * $coverageLength / $macLength[$macId];
      if ($tableRows != '') {
        $tableRows .= ' UNION ALL ';
      }
      $macStart = $allMacStart[$key];
      $macEnd = $allMacEnd[$key];
      $micStart = $allMicStart[$key];
      $micEnd = $allMicEnd[$key];
      if ($tableRows == '') {
        $tableRows .=
          "SELECT " .
          "{$covId} AS `cov_id`, " .
          "{$macId} AS `mac_contig_id`, " .
          "{$macStart} AS `mac_start`, " .
          "{$macEnd} AS `mac_end`, " .
          "{$micId} AS `mic_contig_id`, " .
          "{$micStart} AS `mic_start`, " .
          "{$micEnd} AS `mic_end`," .
          "{$coverageLength} AS `length`, " .
          "{$coveragePercent} AS `coverage` ";
      } else {
        $tableRows .=
          "SELECT " .
          "{$covId}, " .
          "{$macId}, " .
          "{$macStart}, " .
          "{$macEnd}, " .
          "{$micId}, " .
          "{$micStart}, " .
          "{$micEnd}, " .
          "{$coverageLength}, " .
          "{$coveragePercent} ";
      }
      $covId += 1;
    }

    tryQuery($link,
      "INSERT INTO `{$to}`.`coverage`
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
        ({$tableRows}) AS `T`
      LEFT JOIN
        `{$to}`.`contig` AS `NMAC`
      ON
        `T`.`mac_contig_id` = `NMAC`.`contig_id`
      LEFT JOIN
        `{$to}`.`contig` AS `NMIC`
      ON
        `T`.`mic_contig_id` = `NMIC`.`contig_id`;"
    );
  }
  printf("Success: inserted `{$to}`.`coverage`\n");
}

function updateCoverage($link, $to) {
  printf("Starting: update `{$to}`.`coverage`.`mds_num`\n");
  // update the mds_num
  tryQuery($link,
    "UPDATE
      `{$to}`.`coverage` AS `C`,
      (
        SELECT
          `mac_contig_id`,
          `mic_contig_id`,
          COUNT(*) AS `mds_num`
        FROM `{$to}`.`match`
        WHERE `is_preliminary` = 1
        GROUP BY
          `mac_contig_id`,
          `mic_contig_id`
      ) AS `U`
      SET `C`.`mds_num` = `U`.`mds_num`
      WHERE `C`.`mac_contig_id` = `U`.`mac_contig_id`
      AND `C`.`mic_contig_id` = `U`.`mic_contig_id`;"
  );
  printf("Success: updated `{$to}`.`coverage`.`mds_num`\n");

  // update the pointer_num
  printf("Starting: update `{$to}`.`coverage`.`pointer_num`\n");
  tryQuery($link,
    "UPDATE
      `{$to}`.`coverage` AS `C`,
      (
        SELECT
          `mac_contig_id`,
          `mic_contig_id`,
          COUNT(*) AS `pointer_num`
        FROM `{$to}`.`pointer`
        WHERE `is_preliminary` = 1
        GROUP BY
          `mac_contig_id`,
          `mic_contig_id`
      ) AS `U`
      SET `C`.`pointer_num` = `U`.`pointer_num`
      WHERE `C`.`mac_contig_id` = `U`.`mac_contig_id`
      AND `C`.`mic_contig_id` = `U`.`mic_contig_id`;"
  );
  printf("Success: updated `{$to}`.`coverage`.`pointer_num`\n");

  // update the ies_num
  printf("Starting: update `{$to}`.`coverage`.`ies_num`\n");
  tryQuery($link,
    "UPDATE
      `{$to}`.`coverage` AS `C`,
      (
        SELECT
          `mac_contig_id`,
          `mic_contig_id`,
          COUNT(*) AS `ies_num`
        FROM `{$to}`.`ies_strict`
        GROUP BY
          `mac_contig_id`,
          `mic_contig_id`
      ) AS `U`
      SET `C`.`ies_num` = `U`.`ies_num`
      WHERE `C`.`mac_contig_id` = `U`.`mac_contig_id`
      AND `C`.`mic_contig_id` = `U`.`mic_contig_id`;"
  );
  printf("Success: updated `{$to}`.`coverage`.`ies_num`\n");
}

function createCount($link, $to) {
  printf("Starting: create `{$to}`.`count`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`count`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`count` (
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
    ) COMMENT='Summary data about the hits on a contig';"
  );

  printf("Success: created `{$to}`.`count`\n");
}

function insertCount($link, $to) {
  printf("Starting: insert `{$to}`.`count`\n");

  // initialize table
  tryQuery($link,
    "INSERT INTO
      `{$to}`.`count`
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
      `{$to}`.`contig`;"
  );
  printf("Success: intialized table `{$to}`.`count`\n");

  // get `gene`
  tryQuery(
    $link,
    "UPDATE
      `{$to}`.`count` AS `C`,
      (
        SELECT
          `contig_id` AS `contig_id`,
          COUNT(*) AS `num`
        FROM `{$to}`.`gene`
        WHERE `type` = 'gene'
        GROUP BY `contig_id`
      ) AS `G`
    SET `C`.`gene_num` = `G`.`num`
    WHERE `C`.`contig_id` = `G`.`contig_id`;"
  );

  for ($i = 0; $i < 2; $i++) {
    //update mac and mic separately
    $aNuc = $i == 0 ? 'mac' : 'mic';
    $bNuc = $i == 0 ? 'mic' : 'mac';

    // get `mds`, `ies_strict`, `pointer`
    for ($j = 0; $j < 3; $j++) {
      switch ($j) {
        case 0:
          $table = 'match';
          $field = 'mds_num';
          $where = '`is_preliminary` = 1';
          break;
        case 1:
          $table = 'ies_strict';
          $field = 'ies_num';
          $where = '1';
          break;
        case 2:
          $table = 'pointer';
          $field = 'pointer_num';
          $where = '`is_preliminary` = 1';
          break;
      }

      tryQuery(
        $link,
        "UPDATE
          `{$to}`.`count` AS `C`,
          (
            SELECT
              `{$aNuc}_contig_id` AS `contig_id`,
              COUNT(*) AS `num`
            FROM `{$to}`.`{$table}`
            WHERE {$where}
            GROUP BY `{$aNuc}_contig_id`
          ) AS `T`
        SET `C`.`{$field}` = `T`.`num`
        WHERE `C`.`contig_id` = `T`.`contig_id`;"
      );
      printf("Success: updated `{$to}`.`count`.`{$field}` ({$aNuc})\n");
    }

    // get `contig_hits`
    tryQuery($link,
      "UPDATE
        `{$to}`.`count` AS `C`,
        (
          SELECT 
            `{$aNuc}_contig_id` AS `contig_id`,
            COUNT(DISTINCT `{$bNuc}_contig_id`) AS `num`
          FROM
            `{$to}`.`match`
          WHERE
            `is_preliminary` = '1'
          GROUP BY
            `{$aNuc}_contig_id`
        ) AS `T`
      SET
        `C`.`hit_num` = `T`.`num`
      WHERE
        `C`.`contig_id` = `T`.`contig_id`;"
    );
    printf("Success: updated `{$to}`.`contig`.`contig_hits` ({$aNuc})\n");

    // update remaining properties
    tryQuery($link,
      "UPDATE
        `{$to}`.`count` AS `C`,
        (
          SELECT
            `{$aNuc}_contig_id` as `contig_id`,
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
            `{$to}`.`properties`
          GROUP BY
            `{$aNuc}_contig_id` 
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
        `C`.`contig_id` = `P`.`contig_id`;"
    );

    printf("Success: updated additional properties in `{$to}`.`count` ({$aNuc})\n");
  }
}

// #track name="prec_eliminated_sequences" description="intervals comprising complement of the precursor segments labelled comp_[left-flanking-segment-prod-id]_[left-flanking-segment-index]_[right-flanking-segment-prod-id]_[right-flanking-segment-index]" itemRgb-"On"
// OXYTRI_MIC_33550	1	8703	none_0_Contig4386.0_-1	0	+	1	8703	255,0,153
// sdrap_oxy_mac2012_100720_prec_eliminated_sequences.bed
function makeIesBed($link, $to, $from, $type) {
  $handle = fopen('../output/' . $from . "_prec_{$type}_ies_sequences.bed", 'w');

  // write header
  fwrite(
    $handle,
    $type == 'strict' ?
      "#track name=\"prec_strict_ies_sequences\" ".
      "description=\"Strict IESs. Intervals on the precursor segment between two consecutive matches of a product segment " .
      "that do not overlap matches of any other product segment. " .
      "Labelled [product-name]_[left-flanking-match-indexes]_[right-flanking-match-indexes]\" " .
      "itemRgb-\"On\"\n"
    :
      "#track name=\"prec_weak_ies_sequences\" ".
      "description=\"Weak IESs. Intervals on the precursor segment between two consecutive matches of a product segment." .
      "Labelled [product-name]_[left-flanking-match-indexes]_[right-flanking-match-indexes]\" " .
      "itemRgb-\"On\"\n"
  );

  $numIes = mysqli_fetch_row(tryQuery($link, "SELECT COUNT(*) from `{$to}`.`ies_{$type}`;"));
  $numIes = $numIes[0];
  
  $rowStep = 10000;
  for ($i = 0; $i < $numIes; $i += $rowStep) {
    $iesQuery = tryQuery($link,
      "SELECT
        `mic_name`,
        `mic_start`,
        `mic_end`,
        `mac_name`,
        `left_index`,
        `right_index`
      FROM `{$to}`.`ies_{$type}`
      LIMIT $i, {$rowStep};"
    );

    while ($assoc = mysqli_fetch_assoc($iesQuery)) {
      fwrite(
        $handle,
        "{$assoc['mic_name']}\t{$assoc['mic_start']}\t{$assoc['mic_end']}\t" . // chrom chromStart chromEnd
        "{$assoc['mac_name']}_{$assoc['left_index']}_{$assoc['right_index']}\t" . // name
        "0\t+\t{$assoc['mic_start']}\t{$assoc['mic_end']}\t255,0,153\n" //score strand thickStart thickEnd itemRGB
      );
    }
  }
}

function createAlias($link, $to) {
  printf("Starting: create `{$to}`.`alias`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`alias`;");
  tryQuery(
    $link,
    "CREATE TABLE `{$to}`.`alias` (
      `alias_id` INT NOT NULL AUTO_INCREMENT COMMENT 'primary key for the table',
      `id` INT NOT NULL COMMENT 'primary key of the table for the corresponding sequence',
      `name` VARCHAR(50) NOT NULL COMMENT 'the primary name of the sequence',
      `alias` VARCHAR(50) NOT NULL COMMENT 'alias of the sequence',
      `table` VARCHAR(50) NOT NULL COMMENT 'the table that `id` refers to',
      `type` VARCHAR(50) DEFAULT NULL COMMENT 'the type of alias',
      PRIMARY KEY (`alias_id`),
      KEY `alias` (`alias`)
    ) COMMENT='table of sequence aliases';"
  );

  printf("Success: created `alias`\n");
}

function insertAlias_gene($link, $to) {
  printf("Starting: update `{$to}`.`alias` with `{$to}`.`gene`\n");

  tryQuery($link,
    "INSERT INTO `{$to}`.`alias`
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
    FROM `{$to}`.`gene` AS `G`
    WHERE `G`.`type` = 'gene'"
  );

  printf("Success: updated `{$to}`.`alias` with `{$to}`.`gene`\n");
}

function insertAlias_contig($link, $to) {
  printf("Starting: update `{$to}`.`alias` with `{$to}`.`contig`\n");

  tryQuery($link,
    "INSERT INTO `{$to}`.`alias`
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
    FROM `{$to}`.`contig` AS `C`;"
  );

  printf("Success: updated `{$to}`.`alias` with `{$to}`.`contig`\n");
}


function updateAliasWithFile($link, $to, $file, $table, $nucleus) {
  printf("Starting: update `{$to}`.`alias` with {$file} in table '{$table}'\n");
  
  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`temp_alias`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`temp_alias` (
      `name` VARCHAR(50),
      `alias` VARCHAR(50),
      `type` VARCHAR(50),
      KEY (`name`)
    );"
  );
  
  printf("Success: created `{$to}`.`temp_alias`\n");
  
  $allRows = readCSV($file);
  $header = $allRows[0];
  $step = 10000;
  for ($i = 1; $i < count($allRows); $i += $step) {
    $values = array();
    for ($j = $i; $j < min($i + $step, count($allRows)); $j++) {
      for ($k = 1; $k < count($header); $k++) {
        if ($allRows[$j][$k] != 'NA') {
          $values[] = "('{$allRows[$j][0]}','{$allRows[$j][$k]}','{$header[$k]}')";
        }
      }
    }
    $values = implode(',', $values);
    tryQuery($link, "INSERT INTO `{$to}`.`temp_alias` VALUES {$values};");
  }

  if ($table == 'contig') {
    $name_column = 'name';
    $nucleus_column = 'nucleus';
  } else if ($table == 'gene') {
    $name_column = 'attr_id';
    $nucleus_column = 'contig_nucleus';
  } else {
    printf("Error: unknown table '{$table}'");
    exit();
  }

  tryQuery($link,
    "INSERT INTO
      `{$to}`.`alias`
      (
        `id`,
        `name`,
        `alias`,
        `table`,
        `type`
      )
      SELECT
        `T`.`{$table}_id`,
        `U`.`name`,
        `U`.`alias`,
        '{$table}',
        `U`.`type`
      FROM `{$to}`.`temp_alias` AS `U`
      INNER JOIN `{$to}`.`{$table}` AS `T`
      ON `T`.`{$name_column}` = `U`.`name`
      AND `T`.`{$nucleus_column}` = '{$nucleus}';"
  );

  printf("Success: updated `{$to}`.`alias` with {$file}\n");
}

function createVariant($link, $to) {
  printf("Starting: create `{$to}`.`variant`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`variant`;");
  tryQuery(
    $link,
    "CREATE TABLE `{$to}`.`variant` (
      `contig_id` INT NOT NULL AUTO_INCREMENT COMMENT 'primary key for the table',
      `name` VARCHAR(50) NOT NULL COMMENT 'the primary name of the sequence',
      `variant` VARCHAR(500) NOT NULL COMMENT 'variants/isoforms of the sequence',
      PRIMARY KEY (`contig_id`)
    ) COMMENT='Variants/isoforms of sequences';"
  );

  printf("Success: created `variant`\n");
}

function updateVariantWithFile($link, $to, $file) {
  printf("Starting: update `{$to}`.`variant` with {$file}\n");
  
  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`temp_variant`;");
  tryQuery($link,
  "CREATE TABLE `{$to}`.`temp_variant` (
      `name` VARCHAR(50),
      `variant` VARCHAR(500),
      KEY (`name`)
    );"
  );
  
  printf("Success: created `{$to}`.`temp_variant`\n");
  
  $allRows = readCSV($file);
  $step = 10000;
  for ($i = 1; $i < count($allRows); $i += $step) {
    $values = array();
    for ($j = $i; $j < min($i + $step, count($allRows)); $j++) {
      $values[] = "('{$allRows[$j][0]}','{$allRows[$j][1]}')";
    }
    $values = implode(',', $values);
    tryQuery($link, "INSERT INTO `{$to}`.`temp_variant` VALUES {$values};");
  }
  
  tryQuery($link,
    "INSERT INTO
      `{$to}`.`variant`
      (
        `contig_id`,
        `name`,
        `variant`
      )
      SELECT
        `C`.`contig_id`,
        `C`.`name`,
        `V`.`variant`
      FROM `{$to}`.`temp_variant` AS `V`
      INNER JOIN `{$to}`.`contig` AS `C`
      ON `V`.`name` = `C`.`name`;"
  );

  printf("Success: updated `{$to}`.`alias` with {$file}\n");
}

function createStats($link, $to) {
  printf("Starting: create `{$to}`.`stats`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`stats`;");
  tryQuery(
    $link,
    "CREATE TABLE `{$to}`.`stats` (
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
    ) COMMENT='Summary statistics of the database';"
  );

  printf("Success: created `{$to}`.`stats`\n");
}

function insertStats($link, $to) {
  printf("Starting: insert `{$to}`.`stats`\n");

  $stats = array();

  foreach (array('mac', 'mic') as $nucleus) {
    // contig query
    $assoc = mysqli_fetch_assoc(tryQuery($link,
      "SELECT
          SUM( ( ( `telomere_five_start` > 0 ) + ( `telomere_three_start` > 0 ) ) = 2 ) AS `two_telomere`,
          SUM( ( ( `telomere_five_start` > 0 ) + ( `telomere_three_start` > 0 ) ) = 1 ) AS `one_telomere`,
          SUM( ( ( `telomere_five_start` > 0 ) + ( `telomere_three_start` > 0 ) ) = 0 ) AS `zero_telomere`,
          COUNT(*) AS `total`
        FROM `{$to}`.`contig`
        WHERE `nucleus` = '{$nucleus}';"
    ));
    $stats["contig_{$nucleus}_two_telomere"] = $assoc['two_telomere'];
    $stats["contig_{$nucleus}_one_telomere"] = $assoc['one_telomere'];
    $stats["contig_{$nucleus}_zero_telomere"] = $assoc['zero_telomere'];
    $stats["contig_{$nucleus}_total"] = $assoc['total'];

    // gene query
    $assoc = mysqli_fetch_assoc(tryQuery($link,
      "SELECT
        COUNT(*) AS `total`
      FROM `{$to}`.`gene` AS `G`
      INNER JOIN `{$to}`.`contig` AS `C`
        ON `C`.`contig_id` = `G`.`contig_id`
      WHERE `G`.`type` = 'gene'
      AND `C`.`nucleus` = '{$nucleus}';"
    ));
    $stats["gene_{$nucleus}_total"] = $assoc['total'];

    // MDS query
    $assoc = mysqli_fetch_assoc(tryQuery($link,
      "SELECT
      COUNT(*) AS `total`
      FROM `{$to}`.`match`
      WHERE `is_preliminary` = 1;"
    ));
    $stats["mds_{$nucleus}_total"] = $assoc['total'];

    // IES strict query
    if ($nucleus == 'mic') {
      $assoc = mysqli_fetch_assoc(tryQuery($link,
        "SELECT
        COUNT(*) AS `total`
        FROM `{$to}`.`ies_strict`;"
      ));
      $stats["ies_{$nucleus}_total"] = $assoc['total'];
    } else {
      $stats["ies_{$nucleus}_total"] = 0;
    }

    // pointer query
    $assoc = mysqli_fetch_assoc(tryQuery($link,
      "SELECT
      COUNT(*) AS `total`
      FROM `{$to}`.`pointer`;"
    ));
    $stats["pointer_{$nucleus}_total"] = $nucleus == 'mic' ? 2*$assoc['total'] : $assoc['total'];
  }

  // properties query
  $assoc = mysqli_fetch_assoc(tryQuery($link,
    "SELECT
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
    FROM `{$to}`.`properties`;"
  ));
  $stats['properties_non_gapped'] = $assoc['non_gapped'];
  $stats['properties_non_overlapping'] = $assoc['non_overlapping'];
  $stats['properties_non_repeating'] = $assoc['non_repeating'];
  $stats['properties_exceeded_clique_limit'] = $assoc['exceeded_clique_limit'];
  $stats['properties_weakly_complete'] = $assoc['weakly_complete'];
  $stats['properties_strongly_complete'] = $assoc['strongly_complete'];
  $stats['properties_weakly_consecutive'] = $assoc['weakly_consecutive'];
  $stats['properties_strongly_consecutive'] = $assoc['strongly_consecutive'];
  $stats['properties_weakly_ordered'] = $assoc['weakly_ordered'];
  $stats['properties_strongly_ordered'] = $assoc['strongly_ordered'];
  $stats['properties_weakly_non_scrambled'] = $assoc['weakly_non_scrambled'];
  $stats['properties_strongly_non_scrambled'] = $assoc['strongly_non_scrambled'];
  $stats['properties_total'] = $assoc['total'];

  $fields = '(' . implode(',', array_map(function($x) {return "`{$x}`";}, array_keys($stats))) . ')';
  $values = '(' . implode(',', array_map(function($x) {return "'{$x}'";}, array_values($stats))) . ')';

  // insert values
  tryQuery($link,
    "INSERT INTO `{$to}`.`stats` {$fields}
    VALUES {$values};"
  );
  
  printf("Success: inserted `{$to}`.`stats`\n");
}

function createProtein($link, $to) {
  printf("Starting: create `{$to}`.`protein`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`protein`;");
  tryQuery($link, 
    "CREATE TABLE `{$to}`.`protein` (
      `prot_id` int(11) NOT NULL AUTO_INCREMENT COMMENT 'Index of the `protein` table',
      `gene_id` int(11) DEFAULT NULL COMMENT 'Index of the corresponding mRNA feature in the table `gene`',
      `attr_id` varchar(50) DEFAULT NULL COMMENT 'Name of the corresponding mRNA feature',
      `length` int(11) NOT NULL COMMENT 'Length of the protein sequence in amino acids',
      `sequence` longtext COMMENT 'Sequence of amino acids for the protein',
      PRIMARY KEY (`prot_id`),
      UNIQUE KEY `gene_id` (`gene_id`),
      KEY `attr_id` (`attr_id`)
    ) COMMENT='Table of Predicted Proteins'"
  );

  printf("Success: created `{$to}`.`protein`\n");
}

function insertProtein($link, $to, $file) {
  printf("Starting: insert `{$to}`.`protein` with '{$file}'\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`temp_protein`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`temp_protein` (
      `attr_id` VARCHAR(50),
      `sequence` longtext,
      KEY (`attr_id`)
    );"
  );

  $allRows = readCSV($file);
  $step = 1000;
  for ($i = 1; $i < count($allRows); $i += $step) {
    $values = array();
    for ($j = $i; $j < min($i + $step, count($allRows)); $j++) {
      $values[] = "('{$allRows[$j][0]}','{$allRows[$j][1]}')";
    }
    $values = implode(',', $values);
    tryQuery($link, "INSERT INTO `{$to}`.`temp_protein` VALUES {$values};");
  }

  tryQuery($link,
    "INSERT INTO
    `{$to}`.`protein`
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
    FROM `{$to}`.`temp_protein` AS `T`
    INNER JOIN `{$to}`.`gene` AS `G`
    ON `G`.`attr_id` = `T`.`attr_id`;"
  );

  printf("Success: inserted `{$to}`.`protein` with '{$file}'\n");
}

function createGoMap($link, $to) {
  printf("Starting: create `{$to}`.`go_map`\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`go_map`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`go_map` (
      `go_map_id` INT NOT NULL AUTO_INCREMENT COMMENT 'Index of the GO terms map',
      `gene_id` INT DEFAULT NULL COMMENT 'Primary key of the `gene` table for the corresponding gene',
      `gene_attr_id` VARCHAR(50) NOT NULL COMMENT 'Name of the gene',
      `accession` VARCHAR(10) NOT NULL,
      PRIMARY KEY (`go_map_id`),
      KEY `gene_id` (`gene_id`),
      KEY `accession` (`accession`)
    ) COMMENT='GO annotations for genes';"
  );

  printf("Success: created `{$to}`.`go_map`\n");
}

function insertGoMap($link, $to, $file) {
  printf("Starting: insert `{$to}`.`go_map` from file '{$file}'\n");

  tryQuery($link, "DROP TABLE IF EXISTS `{$to}`.`temp_go_map`;");
  tryQuery($link,
    "CREATE TABLE `{$to}`.`temp_go_map` (
      `gene_attr_id` VARCHAR(50) NOT NULL,
      `accession` VARCHAR(10) NOT NULL,
      KEY `gene_attr_id` (`gene_attr_id`)
    )"
  );

  $allRows = readCSV($file);
  $step = 100000;
  for ($i = 1; $i < count($allRows); $i += $step) {
    $rows = makeSqlRows(array_slice($allRows, $i, $step));
    tryQuery($link, "INSERT INTO `{$to}`.`temp_go_map` VALUES {$rows};");
  }

  tryQuery($link,
    "INSERT INTO `{$to}`.`go_map`
      (`gene_id`, `gene_attr_id`, `accession`)
    SELECT
      `G`.`gene_id`, `TGO`.`gene_attr_id`, `TGO`.`accession`
    FROM `{$to}`.`temp_go_map` AS `TGO`
    INNER JOIN `{$to}`.`gene` AS `G` ON
      `G`.`attr_id` = `TGO`.`gene_attr_id`;"
  );

  // tryQuery($link, "DROP TABLE `{$to}`.`temp_go_map`;");

  printf("Success: inserted `{$to}`.`go_map` from file '{$file}'\n");
}

function createGoTerms($link, $to) {
  printf("Starting: create `{$to}`.`go_map`\n");

  tryQuery($link,
    "CREATE TABLE `go_terms` (
      `accession` VARCHAR(10) NOT NULL COMMENT 'Gene ontology accession',
      `name` VARCHAR(200) DEFAULT NULL COMMENT 'Gene ontology term name',
      `ontology` VARCHAR(100) DEFAULT NULL COMMENT 'GO term ontology',
      PRIMARY KEY (`accession`)
    ) COMMENT='Gene Ontology classifications;"
  );

  printf("Success: created `{$to}`.`go_map`\n");
}

function insertGoTerms($link, $to, $file) {
  printf("Starting: create `{$to}`.`go_terms`\n");

  $allRows = readCSV($file);
  $step = 10000;
  for ($i = 1; $i < count($allRows); $i += $step) {
    $values = array();
    for ($j = $i; $j < min($i + $step, count($allRows)); $j++) {
      $values[] = '(' . implode(',', $allRows[$j]) . ')';
    }
    $values = implode(',', $values);
    tryQuery($link, "INSERT INTO `{$to}`.`go_terms` VALUES {$values};");
  }

  printf("Success: created `{$to}`.`go_terms`\n");
}

function dropTempTables($link, $to) {
  $tempTables = array(
    "`{$to}`.`temp_name`",
    "`{$to}`.`temp_alias`",
    "`{$to}`.`temp_variant`",
    "`{$to}`.`temp_protein`",
  );
  foreach ($tempTables as $table) {
    printf("Starting: drop {$table}\n");
    tryQuery($link, "DROP TABLE IF EXISTS {$table};");
    printf("Success: dropped {$table}\n");
  }
}

function addToDatabaseDirectory($link, $to, $name, $description, $organism, $downloadDir, $assembly) {
  printf("Starting: add `{$to}` to `mds_ies_db`.`db` directory\n");

  tryQuery(
    $link,
    "INSERT INTO `db` (`name`, `description`, `sql_name`, `organism`, `download_dir`, `assembly`)
    VALUES ('{$name}', '{$description}', '{$to}', '{$organism}', '{$downloadDir}', '{$assembly}');"
  );

  printf("Success: add `{$to}` to `mds_ies_db` directory\n");
}

function dumpTable($link, $to, $table, $fileName) {
  printf("Starting: dump `{$to}`.`$table` to '{$fileName}'\n");

  $rowsPerQuery = 10000;
  $numRows = tryQuery($link, "SELECT COUNT(*) FROM `{$to}`.`{$table}`");
  $numRows = mysqli_fetch_row($numRows);
  $numRows = $numRows[0];
  
  $fileHandle = fopen($fileName, 'w');
  $writeHeader = true;
  $currRow = 0;
  while ($currRow < $numRows) {
    $query = tryQuery($link, "SELECT * FROM `{$to}`.`{$table}` LIMIT {$currRow}, {$rowsPerQuery};");
    while ($assoc = mysqli_fetch_assoc($query)) {
      if ($writeHeader) {
        fputcsv($fileHandle, array_keys($assoc));
        $writeHeader = false;
      }
      fputcsv($fileHandle, array_values($assoc));
      $currRow++;
    }
  }
  fclose($fileHandle);

  printf("Success: dumped `{$to}`.`$table` to '{$fileName}'\n");
}

function main() {
  header("Content-type: text/plain");

  $link = tryConnectDatabase('mds_ies_db');
  $to = $_GET['to'];
  $from = $_GET['from'];
  $macRegex = $_GET['macRegex'];
  $micRegex = $_GET['micRegex'];
  $stage = $_GET['stage'];

  switch ($stage) {
    case 'create':
      checkFromDatabaseExists($link, $from);
      createDatabase($link, $to);
      break;

    case 'tempName':
      createTempName($link, $to);
      insertTempName($link, $to, $from, $macRegex, $micRegex);
      break;

    case 'insertTempName_mac2020_genbankV1ToOxydb':
      createTempName($link, $to);
      insertTempName_oxyTriJrb310_mac2020_genbankV1ToOxydb($link, $to, $from, $micRegex);
      break;

    case 'contig':
      createContig($link, $to);
      insertContig($link, $to, $from);
      break;

    case 'hsp':
      createHsp($link, $to);
      insertHsp($link, $to, $from);
      break;

    case 'match':
      createMatch($link, $to);
      insertMatch($link, $to, $from);
      break;

    case 'coverage':
      createCoverage($link, $to);
      insertCoverage($link, $to);
      updateCoverage($link, $to);
      break;

    case 'pointer':
      createPointer($link, $to);
      insertPointer($link, $to, $from);
      break;

    case 'properties':
      createProperties($link, $to);
      insertProperties($link, $to, $from);
      break;

    case 'parameter':
      createParameter($link, $to);
      insertParameter($link, $to, $from);
      break;

    case 'count':
      createCount($link, $to);
      insertCount($link, $to);
      break;

    case 'protein|oxyTriJrb310_mac2012':
      createProtein($link, $to);
      insertProtein($link, $to, 'Oxytricha_trifallax_022112_aa.csv');
      break;
    
    case 'protein|oxyTriJrb310_mac2020':
      createProtein($link, $to);
      insertProtein($link, $to, 'oxytrijrb310pacbio.faa.csv');
      break;
    
    case 'makeIesBed_strict':
      makeIesBed($link, $to, $from, 'strict');
      break;
    
    case 'makeIesBed_weak':
      makeIesBed($link, $to, $from, 'weak');
      break;

    case 'alias':
      createAlias($link, $to);
      insertAlias_contig($link, $to);
      insertAlias_gene($link, $to);
      break;

    case 'updateAlias_oxyTriJrb310_mac2012':
      updateAliasWithFile($link, $to, 'oxy_tri_jrb310_mac_2012_alias.csv', 'contig', 'mac');
      updateAliasWithFile($link, $to, 'oxy_tri_jrb310_mac_gene_2012_alias.csv', 'gene', 'mac');
      break;

    case 'updateAlias_oxyTriJrb310_mac2020':
      updateAliasWithFile($link, $to, 'oxy_tri_jrb310_mac_2020_alias.csv', 'contig', 'mac');
      updateAliasWithFile($link, $to, 'oxy_tri_jrb310_mac_gene_2020_alias.csv', 'gene', 'mac');
      break;

    case 'updateAlias_oxyTriJrb310_mic2014':
      updateAliasWithFile($link, $to, 'oxy_tri_jrb310_mic_2014_alias.csv', 'contig', 'mic');
      break;

    case 'variant':
      createVariant($link, $to);
      break;

    case 'updateVariant_oxyTriJrb310_mac2020':
      updateVariantWithFile($link, $to, 'oxy_tri_jrb310_mac_2020_variant.csv');
      break;

    case 'stats':
      createStats($link, $to);
      insertStats($link, $to);
      break;

    case 'goMap_2012':
      createGoMap($link, $to);
      insertGoMap($link, $to, 'oxy_2012_go_map.csv');
      break;

    case 'goTerms_2012':
      createGoMap($link, $to);
      insertGoMap($link, $to, 'oxy_2012_go_terms.csv');
      break;
    
    case 'addToDatabaseDirectory_oxyTriJrb310_2012':
    case 'addToDatabaseDirectory_oxyTriJrb310_2020':
      $year = explode('_', $stage);
      $year = $year[1];
      $name = "Oxytricha Trifallax JRB310 (MAC {$year}/MIC 2014)";
      $description = "Oxytricha Trifallax JRB310 - SDRAP Annotation - MAC {$year} - MIC 2014";
      $organism = 'oxy_tri_jrb_310';
      $downloadDir = $to;
      $assembly = "mac_{$year},mic_2014";

      addToDatabaseDirectory($link, $to, $name, $description, $organism, $downloadDir, $assembly);
      break;

    case 'dropTempTables':
      dropTempTables($link, $to);
      break;
    
    case 'dumpTables':
      if (!is_dir("../output/{$to}")) {
        mkdir("../output/{$to}", 0777, true);
      }
      dumpTable($link, $to, 'alias', "../output/{$to}/alias.csv");
      dumpTable($link, $to, 'contig', "../output/{$to}/contig.csv");
      dumpTable($link, $to, 'count', "../output/{$to}/count.csv");
      dumpTable($link, $to, 'coverage', "../output/{$to}/coverage.csv");
      dumpTable($link, $to, 'gene', "../output/{$to}/gene.csv");
      dumpTable($link, $to, 'ies_strict', "../output/{$to}/ies_strict.csv");
      dumpTable($link, $to, 'ies_weak', "../output/{$to}/ies_weak.csv");
      dumpTable($link, $to, 'match', "../output/{$to}/match.csv");
      dumpTable($link, $to, 'parameter', "../output/{$to}/parameter.csv");
      dumpTable($link, $to, 'pointer', "../output/{$to}/pointer.csv");
      dumpTable($link, $to, 'properties', "../output/{$to}/properties.csv");
      dumpTable($link, $to, 'protein', "../output/{$to}/protein.csv");
      dumpTable($link, $to, 'stats', "../output/{$to}/stats.csv");
      dumpTable($link, $to, 'variant', "../output/{$to}/variant.csv");
      break;

    default:
      printf("Error: unknown stage '{$stage}'\n");
      break;
  }
}

main();

?>