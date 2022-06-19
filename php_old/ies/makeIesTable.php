<?php

include '../databaseConnect.php';
include '../sqlUtil.php';

function getMicContigs($link) {
  $result = tryQuery(
    $link,
    "SELECT
      `contig_id`,
      `name`
    FROM `contig`
    WHERE `nucleus` = 'mic';"
  );
  $micContigs = array();
  while ($assoc = mysqli_fetch_assoc($result)) {
    $micContigs[] = $assoc;
  }
  return $micContigs;
}

function makeMatch($start, $end, $mac) {
  return
    array(
      'start' => $start,
      'end' => $end,
      'mac' => $mac
    );
}

function makeMac($contig_id, $name, $index, $orientation) {
  return
    array(
      'contig_id' => $contig_id,
      'name' => $name,
      'index' => $index,
      'orientation' => $orientation
    );
}

function getMatches($link, $micContig) {
  $result = tryQuery(
    $link, 
    "SELECT
      `mic_start`,
      `mic_end`,
      `mac_contig_id`,
      `mac_name`,
      ABS(`index`) AS `index`,
      `orientation`
    FROM
      `match`
    WHERE
      `mic_contig_id` = '{$micContig['contig_id']}'
    AND
      `is_preliminary` = 1;"
  );
  $matches = array();
  while ($assoc = mysqli_fetch_assoc($result)) {
    $matches[] =
      makeMatch(
        $assoc['mic_start'],
        $assoc['mic_end'],
        array(
          $assoc['mac_name'] => makeMac(
            $assoc['mac_contig_id'],
            $assoc['mac_name'],
            array($assoc['index']),
            array($assoc['orientation'])
          )
        )
      );
  }
  return $matches;
}

function makeMergedMatch($start, $end, $start_macs, $end_macs) {
  return
    array(
      'start' => $start,
      'end' => $end,
      'start_macs' => $start_macs,
      'end_macs' => $end_macs
    );
}

function toMergedMatch($match) {
  return
    makeMergedMatch(
      $match['start'],
      $match['end'],
      $match['mac'],
      $match['mac']
    );
}

function mergeable($start1, $end1, $start2, $end2) {
  if ($end1 + 1 < $start2) return false;
  if ($end2 + 1 < $start1) return false;
  return true;
}

function mergeMacs($mac1, $mac2) {
  if (($mac1['name'] != $mac2['name']) || ($mac1['contig_id'] != $mac2['contig_id'])) {
    printf(strval(__FILE__) . ': ' . strval(__LINE__) . ': MACs are not mergable');
    exit();
  }

  $indexMap1 = array_combine($mac1['index'], $mac1['orientation']);
  $indexMap2 = array_combine($mac2['index'], $mac2['orientation']);
  $indexMapBoth = $indexMap1 + $indexMap2;
  ksort($indexMapBoth, SORT_NUMERIC);

  return
    makeMac(
      $mac1['contig_id'], // should be the same for both
      $mac1['name'], // should be the same for both
      array_keys($indexMapBoth),
      array_values($indexMapBoth)
      // array_unique(array_merge($mac1['index'], $mac2['index']))
    );
}

function mergeAllMacs($macs1, $macs2) {
  $allNames = array_unique(array_merge(array_keys($macs1), array_keys($macs2)));
  $merged = array();
  foreach ($allNames as $name) {
    if (isset($macs1[$name]) && isset($macs2[$name])) {
      $merged[$name] = mergeMacs($macs1[$name], $macs2[$name]);
    } else if (isset($macs1[$name])) {
      $merged[$name] = $macs1[$name];
    } else if (isset($macs2[$name])) {
      $merged[$name] = $macs2[$name];
    } else {
      printf("Error: unreachable case %s:%d\n", __FUNCTION__, __LINE__);
      exit(); 
    }
  }
  return $merged;
}

// For the non-strict matching separate out by the MAC contig
function groupMatchesByMac($matches) {
  $groupedMatches = array();
  foreach ($matches as $match) {
    $names = array_keys($match['mac']);
    $groupedMatches[$names[0]][] = $match;
  }
  return array_values($groupedMatches);
}

function mergeMatches($merged1, $merged2) {
  if ($merged1['start'] < $merged2['start']) {
    $start = $merged1['start'];
    $start_macs = $merged1['start_macs'];
  } else if ($merged1['start'] > $merged2['start']) {
    $start = $merged2['start'];
    $start_macs = $merged2['start_macs'];
  } else if ($merged1['start'] == $merged2['start']) {
    $start = $merged1['start'];
    $start_macs =
      mergeAllMacs(
        $merged1['start_macs'],
        $merged2['start_macs']
      );
  }

  if ($merged1['end'] < $merged2['end']) {
    $end = $merged2['end'];
    $end_macs = $merged2['end_macs'];
  } else if ($merged1['end'] > $merged2['end']) {
    $end = $merged1['end'];
    $end_macs = $merged1['end_macs'];
  } else if ($merged1['end'] == $merged2['end']) {
    $end = $merged1['end'];
    $end_macs =
      mergeAllMacs(
        $merged1['end_macs'],
        $merged2['end_macs']
      );
  }
  
  return
    array(
      'start' => $start,
      'end' => $end,
      'start_macs' => $start_macs,
      'end_macs' => $end_macs
    );
}

function cmpByStart($x, $y) {
  if ($x['start'] < $y['start']) return -1;
  if ($x['start'] > $y['start']) return 1;
  return 0;
}

function mergeAllMatches($matches) {
  usort($matches, 'cmpByStart');
  $mergedMatches = array();
  foreach ($matches as $match) {
    $merged = toMergedMatch($match);
    if (!empty($mergedMatches)) {
      $lastMerged = $mergedMatches[count($mergedMatches) - 1];
      if (mergeable($lastMerged['start'], $lastMerged['end'], $merged['start'], $merged['end'])) {
        $mergedMatches[count($mergedMatches) - 1] = mergeMatches($lastMerged, $merged);
      } else {
        $mergedMatches[] = $merged;
      }
    } else {
      $mergedMatches[] = $merged;
    }
  }
  return $mergedMatches;
}

function makeGroupedIes($start, $end, $macs) {
  return array('start' => $start, 'end' => $end, 'macs' => $macs);
}

function getGroupedIes($mergedMatches) {
  $ies = array();
  for ($i = 0; $i < count($mergedMatches)-1; $i++) {
    $intersectNames =
      array_intersect(
        array_keys($mergedMatches[$i]['end_macs']),
        array_keys($mergedMatches[$i+1]['start_macs'])
      );
    if (!empty($intersectNames)) {
      $intersect = array();
      foreach ($intersectNames as $name) {
        $intersect[] = array(
          'name' => $name,
          'contig_id' => $mergedMatches[$i]['end_macs'][$name]['contig_id'],
          'left_index' => $mergedMatches[$i]['end_macs'][$name]['index'],
          'right_index' => $mergedMatches[$i+1]['start_macs'][$name]['index'],
          'left_orientation' => $mergedMatches[$i]['end_macs'][$name]['orientation'],
          'right_orientation' => $mergedMatches[$i+1]['start_macs'][$name]['orientation'],
        );
      }
      $ies[] = makeGroupedIes($mergedMatches[$i]['end'] + 1, $mergedMatches[$i+1]['start'] - 1, $intersect);
    }
  }
  return $ies;
}

function quoteVal($val) {
  return "'{$val}'";
}

// function getSqlIesValues(&$values, &$currIesId, $micContig, $groupedIes) {
//   foreach ($groupedIes as $ies) {
//     foreach ($ies['macs'] as $mac) {
//       $leftIndexStr = implode(',', $mac['left_index']);
//       $rightIndexStr = implode(',', $mac['right_index']);
//       $values[] =
//         '(' .
//           quoteVal($currIesId) . ',' .
//           quoteVal($micContig['contig_id']) . ',' .
//           quoteVal($micContig['name']) . ',' .
//           quoteVal($ies['start']) . ',' .
//           quoteVal($ies['end']) . ',' .
//           quoteVal($ies['end'] - $ies['start'] + 1) . ',' .
//           quoteVal($mac['contig_id']) . ',' .
//           quoteVal($mac['name']) . ',' .
//           quoteVal($leftIndexStr) . ',' .
//           quoteVal($rightIndexStr) .
//         ')';
//       $currIesId++;
//     }
//   }
// }

function getSqlIesValues(&$values, &$currIesId, $micContig, $groupedIes) {
  foreach ($groupedIes as $ies) {
    foreach ($ies['macs'] as $mac) {
      $leftIndexStr = implode(',', $mac['left_index']);
      $rightIndexStr = implode(',', $mac['right_index']);
      $leftOrientationStr = implode(',', $mac['left_orientation']);
      $rightOrientationStr = implode(',', $mac['right_orientation']);
      $values[] =
        '(' .
          quoteVal($currIesId) . ',' .
          quoteVal($micContig['contig_id']) . ',' .
          quoteVal($micContig['name']) . ',' .
          quoteVal($ies['start']) . ',' .
          quoteVal($ies['end']) . ',' .
          quoteVal($ies['end'] - $ies['start'] + 1) . ',' .
          quoteVal($mac['contig_id']) . ',' .
          quoteVal($mac['name']) . ',' .
          quoteVal($leftIndexStr) . ',' .
          quoteVal($rightIndexStr) . ',' .
          quoteVal($leftOrientationStr) . ',' .
          quoteVal($rightOrientationStr) .
        ')';
      $currIesId++;
    }
  }
}

function test() {
  $matches = array(
    makeMatch(0, 100, 'A'), 
    makeMatch(0, 100, 'B'), 
    makeMatch(0, 1000, 'C'), 
    makeMatch(0, 1000, 'E'),
    makeMatch(1002, 1004, 'E'),
  );

  printf("matches:\n");
  var_export($matches);
  printf("\n");

  $mergedMatches = mergeAllMatches($matches);

  printf("mergedMatches:\n");
  var_export($mergedMatches);
  printf("\n");

  $groupedIes = getGroupedIes($mergedMatches);

  printf("groupedIes:\n");
  var_export($groupedIes);
  printf("\n");

  $sqlIesValues = array();
  $currIesId = 0;
  $micContig = array(
    'contig_id' => 1,
    'name' => 'Mic_Contig_1'
  );
  getSqlIesValues($sqlIesValues, $currIesId, $micContig, $groupedIes);

  printf("sqlIesValues:\n");
  var_export($sqlIesValues);
  printf("\n");
}

function main() {
  header("Content-type: text/plain");

  // prevent script from dying after 30 sec
  set_time_limit(0); 
  //ignore_user_abort(true);
  ini_set('max_execution_time', 0);


  $db = $_GET['db'];

  $type = $_GET['type'];

  if (!in_array($type, array('strict', 'weak'))) {
    printf("Error: unknown type; {$type}\n");
    exit();
  }

  $link = tryConnectDatabase($db);

  $micContigs = getMicContigs($link);

  $comment = array(
    'strict' => 'segments on the precursor contigs that lie between 2 product matches that does not overlap an other matches from other contigs',
    'weak' => 'segments on the precursor contigs that lie between 2 product matches',
  );
  // create SQL table
  printf("Starting: create database `ies_{$type}`\n");
  
  tryQuery($link, "DROP TABLE IF EXISTS `ies_{$type}`;");
  tryQuery(
    $link,
    "CREATE TABLE `ies_{$type}` (
      `ies_id` int NOT NULL COMMENT 'primary key for the table',
      `mic_contig_id` int NOT NULL COMMENT 'primary key of the `contig` table for the corresponding precursor sequence',
      `mic_name` varchar(50) NOT NULL COMMENT 'MIC contig name',
      `mic_start` int NOT NULL COMMENT 'position of the first base pair of the precursor segment of the match on the corresponding precursor sequence',
      `mic_end` int NOT NULL COMMENT 'position of the last base pair of the precursor segment of the match on the corresponding precursor sequence',
      `length` int NOT NULL COMMENT 'length of the ies',
      `mac_contig_id` int NOT NULL COMMENT 'primary key of the `contig` table for the corresponding product sequence',
      `mac_name` varchar(50) NOT NULL COMMENT 'MAC contig name',
      `left_index` varchar(200) NOT NULL COMMENT 'Indices of the MDS on the MAC contig flanking the IES on the left',
      `right_index` varchar(200) NOT NULL COMMENT 'Indices of the MDS on the MAC contig flanking the IES on the right',
      `left_orientation` varchar(200) NOT NULL COMMENT 'Orientations of the MDS on the MAC contig flanking the IES on the left',
      `right_orientation` varchar(200) NOT NULL COMMENT 'Orientations of the MDS on the MAC contig flanking the IES on the right',
      PRIMARY KEY (`ies_id`),
      KEY `mic_contig_id` (`mic_contig_id`),
      KEY `mac_contig_id` (`mac_contig_id`)
    ) COMMENT '{$comment[$type]}';"
  );
  printf("Success: created database `ies_{$type}`\n");

  printf("Starting: computing IES's and uploading to database\n");
  $values = array();
  $currIesId = 1;
  $firstUpload_i = 0; // for logging
  $iesCount = 0; // for logging
  for ($i = 0; $i < count($micContigs);  $i++) {
    $micContig = $micContigs[$i];
    $matches = getMatches($link, $micContig);

    if ($type == 'weak') {
      // separate out by MAC and then calculate IES
      $matchesByMac = groupMatchesByMac($matches);
      $groupedIes = array();
      foreach ($matchesByMac as $singleMacMatches) {
        $singleMacMergedMatches = mergeAllMatches($singleMacMatches);
        $singleMacGroupedIes = getGroupedIes($singleMacMergedMatches);
        foreach ($singleMacGroupedIes as $ies) {
          $groupedIes[] = $ies;
        }
      }
    } else { // $type == 'strict'
      // pool all MAC matches together so that there are no overlapping IES/MDS between different MAC contigs
      $mergedMatches = mergeAllMatches($matches);
      $groupedIes = getGroupedIes($mergedMatches);
    }

    getSqlIesValues($values, $currIesId, $micContig, $groupedIes);
    
    // every 100 contigs write to SQL
    if ((($i+1) % 100 === 0) || ($i == count($micContigs) - 1)) {
      printf("Starting: upload ies in contigs %d..%d to database ({$type})\n", $firstUpload_i, $i);
      flush();
      $firstUpload_i = $i + 1;
      $iesCount += count($values);
      $insertQuery = "INSERT INTO `ies_{$type}` VALUES " . implode(', ', $values) . ";";
      $values = array();
      tryQuery($link, $insertQuery);
    }
  }

  printf("Success: uploaded ies in %d contigs, %d IES's total ({$type})\n", $i + 1, $iesCount);
}

main();

?>