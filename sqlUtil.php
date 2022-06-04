<?php

function tryQuery($link, $queryStr) {
  $result = mysqli_query($link, $queryStr);
  if ($result === false || mysqli_errno($link)) {
    printf("Error: %s\n", mysqli_error($link));
    printf("Query: $queryStr\n");
    exit();
  }
  return $result;
}

function makeSqlValue($x) {
  if (is_null($x)) return 'NULL';
  $x = str_replace("'", "''", $x); // escape single quotes
  return "'$x'";
}

function makeSqlRow($values) {
  return '(' . implode(',', array_map('makeSqlValue', $values)) . ')';
}

function makeSqlRows($rows) {
  return implode(',', array_map('makeSqlRow', $rows));
}

?>