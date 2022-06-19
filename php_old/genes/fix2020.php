<?php

header("Content-type: text/plain");

include 'gffUtil.php';

$inHandle = fopen('O_trifallax_2020-upd3.gff3', 'r');
$outHandle = fopen('O_trifallax_2020-upd3.fix.gff3', 'w');

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

while ($line = fgets($inHandle)) {
  $gffLine = parseGffLine($line);

  if (is_null($gffLine)) continue;

  if (!in_array($gffLine['type'], array('gene', 'mRNA', 'exon', 'intron', 'CDS'))) continue;

  if ($gffLine['type'] == 'mRNA')  {
    $gffLine['attrs']['Parent'] = explode('.', $gffLine['attrs']['ID'])[0];
  } else if ($gffLine['type'] == 'exon' || $gffLine['type'] == 'intron') {
    $id = $gffLine['attrs']['ID'];

    $idPrefix = substr($id, 0, strpos($id, '.'));
    $parent = substr($id, strpos($id, '.') + 1);
    
    if (!isset($featureCounts[$parent]))
      $featureCounts[$parent] = array('exon' => 0, 'intron' => 0);

    $idNum = ++$featureCounts[$parent][$gffLine['type']];
    

    $id = "{$idPrefix}{$idNum}.{$parent}";

    $gffLine['attrs']['ID'] = $id;
    $gffLine['attrs']['Parent'] = $parent;
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