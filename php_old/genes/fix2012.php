<?php

header("Content-type: text/plain");

include 'gffUtil.php';

$inHandle = fopen('Oxytricha_trifallax_022112.gff3', 'r');
$outHandle = fopen('Oxytricha_trifallax_022112.fix.gff3', 'w');

$featureCounts = array();

$lines = 0;

if (!$inHandle) {
  printf("Error: unable to open input\n");
  exit();
}

if (!$outHandle) {
  printf("Error: unable to open output\n");
  exit();
}

$mRNA = array();

while ($line = fgets($inHandle)) {
  $gffLine = parseGffLine($line);

  if (is_null($gffLine)) continue;

  if (!in_array($gffLine['type'], array('gene', 'mRNA', 'exon', 'intron', 'CDS'))) continue;

  if ($gffLine['type'] == 'exon' || $gffLine['type'] == 'intron') {
    $id = $gffLine['attrs']['ID'];
    $parent = $gffLine['attrs']['Parent'];

    $idPrefix = substr($id, 0, strpos($id, '_'));
    $idSuffix = substr($id, strpos($id, '_') + 1);
    
    if (!isset($featureCounts[$parent]))
      $featureCounts[$parent] = array('exon' => 0, 'intron' => 0);

    $idNum = ++$featureCounts[$parent][$gffLine['type']];

    $id = "{$idPrefix}{$idNum}_{$idSuffix}";

    $gffLine['attrs']['ID'] = $id;
  }

  fputs($outHandle, makeGffLine($gffLine));
  fputs($outHandle, "\n");
  $lines++;
  if (($lines % 1000) == 0) {
    printf("Lines: %d\n", $lines);
    flush();
  }
}

fclose($inHandle);
fclose($outHandle);

printf("Success\n");

?>