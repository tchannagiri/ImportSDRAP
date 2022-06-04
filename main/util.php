<?php
function readCSV($file) {
  $handle = fopen($file, 'r');
  $allRows = array();
  while ($row = fgetcsv($handle)) {
    $allRows[] = $row;
  }
  return $allRows;
}
?>