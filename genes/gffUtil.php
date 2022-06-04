<?php

function parseGffLine($line) {
  $line = trim($line);

  if (strlen($line) == 0) return null;

  if ($line[0] == '#') return null; // comment line

  $fields = explode("\t", rtrim($line));

  $contig_name = $fields[0];
  $source = $fields[1];
  $type = $fields[2];
  $start = $fields[3];
  $end = $fields[4];
  $score = $fields[5];
  $strand = $fields[6];
  $phase = $fields[7];
  $attrs = $fields[8];

  $attrAssoc = array();
  foreach (explode(';', $attrs) as $attrRecord) {
    if (strpos($attrRecord, '=') === false) continue;
    $attrKeyValue = explode('=', $attrRecord);
    if (count($attrKeyValue) !== 2) continue;
    $attrAssoc[$attrKeyValue[0]] = $attrKeyValue[1];
  }

  return array(
    'contig_name' => $contig_name,
    'source' => $source,
    'type' => $type,
    'start' => $start,
    'end' => $end,
    'score' => $score,
    'strand' => $strand,
    'phase' => $phase,
    'attrs' => $attrAssoc,
  );
}

function makeGffLine($assoc) {
  $attrStrs = array();
  foreach ($assoc['attrs'] as $key => $value) {
    $attrStrs[] = "{$key}={$value}";
  }
  $attrWholeStr = implode(';', $attrStrs);

  $contig_name = $assoc['contig_name'];
  $source = $assoc['source'];
  $type = $assoc['type'];
  $start = $assoc['start'];
  $end = $assoc['end'];
  $score = $assoc['score'];
  $strand = $assoc['strand'];
  $phase = $assoc['phase'];

  return
    "{$contig_name}\t{$source}\t{$type}\t{$start}\t{$end}\t" .
    "{$score}\t{$strand}\t{$phase}\t{$attrWholeStr}";
}

?>