<?php

include "../databaseConnect.php";
include "../sqlUtil.php";
include "gffUtil.php";

function makeGeneTable($link) {
  tryQuery($link, "DROP TABLE IF EXISTS `gene`;");
  tryQuery(
    $link,
    "CREATE TABLE `gene` (
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
     ) COMMENT='Table of genes and gene features';"
  );
}

function getContigs($link) {
  $contigQuery =
    "SELECT
      `contig_id`,
      `nucleus` AS `contig_nucleus`,
      `name` AS `contig_name`
    FROM `contig`;";

  $result = tryQuery($link, $contigQuery);
  $contigs = array();
  while ($assoc = mysqli_fetch_assoc($result)) {
    $contigs[$assoc['contig_name']] = $assoc;
  }
  return $contigs;
}

function insertGeneTable($link, $values) {
  $queryStr =
    "INSERT INTO `gene`
    (
      `contig_id`,
      `contig_name`,
      `contig_nucleus`,
      `source`,
      `type`,
      `start`,
      `end`,
      `length`,
      `score`,
      `strand`,
      `phase`,
      `attr_id`,
      `attr_parent`,
      `attr_name`,
      `attr_note`,
      `attr_parent_gene`
    )
    VALUES ";
  $first = true;
  foreach ($values as $row) {
    if (!$first) {
      $queryStr .= ',';
    } else {
      $first = false;
    }
    $queryStr .=
      '(' .
        makeSqlValue($row['contig_id']) . ',' .
        makeSqlValue($row['contig_name']) . ',' .
        makeSqlValue($row['contig_nucleus']) . ',' .
        makeSqlValue($row['source']) . ',' .
        makeSqlValue($row['type']) . ',' .
        makeSqlValue($row['start']) . ',' .
        makeSqlValue($row['end']) . ',' .
        makeSqlValue($row['length']) . ',' .
        makeSqlValue($row['score']) . ',' .
        makeSqlValue($row['strand']) . ',' .
        makeSqlValue($row['phase']) . ',' .
        makeSqlValue($row['attr_id']) . ',' .
        makeSqlValue($row['attr_parent']) . ',' .
        makeSqlValue($row['attr_name']) . ',' .
        makeSqlValue($row['attr_note']) . ',' .
        makeSqlValue($row['attr_parent_gene']) .
      ')';
  }
  $queryStr .= ';';
  tryQuery($link, $queryStr);
}

function parseGffFile($link, $gffFile) {
  $contigs = getContigs($link); // store this to lookup the contig info

  $handle = fopen($gffFile, 'r');

  if ($handle === false) {
    printf("Error: unable to open file '$gffFile'\n");
    exit();
  }

  $values = array();

  $prevGene_id = 0; // for logging
  $gene_id = 1;
  while ($line = fgets($handle)) {
    if ($gffLine = parseGffLine($line)) {

      $contig_name = $gffLine['contig_name'];
      $source = $gffLine['source'];
      $type = $gffLine['type'];
      $start = intval($gffLine['start']) - 1;
      $end = intval($gffLine['end']) - 1;
      $score = $gffLine['score'] == '.' ? null : $gffLine['score'];
      $strand = $gffLine['strand'];
      $phase = $gffLine['phase'];
      $attrs= $gffLine['attrs'];
      
      if (!in_array($gffLine['type'], array('gene', 'mRNA', 'exon', 'intron', 'CDS')))
        continue;

      $attr_id = $attrs['ID'];
      $attr_parent = isset($attrs['Parent']) ? $attrs['Parent'] : null;     
      $attr_name = isset($attrs['Name']) ? $attrs['Name'] : null;     
      $attr_note = isset($attrs['Note']) ? $attrs['Note'] : null;
      $attr_parent_gene = isset($attrs['Parent_gene']) ? $attrs['Parent_gene'] : null;

      if (!isset($attr_id)) {
        printf("Error: invalid ID in attributes '$attrs'\n");
        var_export($gffLine);
        printf("\n");
        exit();
      }

      $contigInfo = isset($contigs[$gffLine['contig_name']]) ?
        $contigs[$gffLine['contig_name']] : null;

      if (!isset($contigInfo)) {
        //printf("Warning: unable to find contig '$contig_name'\n");
        continue;
      }

      $contig_id = $contigInfo['contig_id'];
      $contig_nucleus = $contigInfo['contig_nucleus'];

      $values[] =
        array(
          'contig_id' => $contig_id,
          'contig_name' => $contig_name,
          'contig_nucleus' => $contig_nucleus,
          'source' => $source,
          'type' => $type,
          'start' => $start,
          'end' => $end,
          'length' => $end - $start + 1,
          'score' => $score,
          'strand' => $strand,
          'phase' => $phase,
          'attr_id' => $attr_id,
          'attr_parent' => $attr_parent,
          'attr_name' => $attr_name,
          'attr_note' => $attr_note,
          'attr_parent_gene' => $attr_parent_gene,
        );

      // upload to SQL every 10000 genes
      if (($gene_id % 10000) === 0) {
        printf("Starting: upload gene features %d..%d to database\n", $prevGene_id, $gene_id);
        flush();
        $prevGene_id = $gene_id + 1;
        insertGeneTable($link, $values);
        $values = array();
      }
      
      $gene_id++;
    }
  }

  if ($gene_id > $prevGene_id) {
    printf("Starting: upload gene features %d..%d to database\n", $prevGene_id, $gene_id-1);
    insertGeneTable($link, $values);
  }

  fclose($handle);
}

function main() {
  header("Content-type: text/plain");

  $db = $_GET['db'];

  printf("Starting: connect to database `{$db}`\n");
  $link = tryConnectDatabase($db);
  printf("Success: connected to database `{$db}\n");

  if ($_GET['stage'] == 'createGene') {
    printf("Starting: create table `gene`\n");
    makeGeneTable($link);
    printf("Success: created table `gene`\n");
    exit();
  }

  switch ($_GET['stage']) {
    case 'oxyTriJrb310_mac2012':
      $file = 'Oxytricha_trifallax_022112.fix.gff3';
      break;
    case 'oxyTriJrb310_mac2020':
      $file = 'O_trifallax_2020-upd3.fix.gff3';
      break;
    case 'oxyTriJrb310_mic2014':
      $file = 'oxy_tri_jrb310_mic_2014.fix.gff';
      break;
    default:
      print("Error: unknown stage: '{$_GET['stage']}'\n");
      exit();
  }

  printf("Starting: parse and upload GFF file '{$file}'\n");
  parseGffFile($link, $file);
  printf("Success: parsed and uploaded GFF file '{$file}'\n");

  mysqli_close($link);
}

main();

?>