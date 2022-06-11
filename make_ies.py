
import collections

import constants
import common_utils
import mysql_utils


def getMicContigs(cursor, db):
  cursor.execute(
    f"""
    SELECT
      `contig_id`,
      `name`
    FROM `{db}`.`contig`
    WHERE `nucleus` = 'mic';
    """
  )
  return cursor.fetchall()

# function makeMatch(start, end, mac) {
#   return
#     array(
#       'start' => start,
#       'end' => end,
#       'mac' => mac
#     );
# }

# function makeMac(contig_id, name, index, orientation) {
#   return
#     array(
#       'contig_id' => contig_id,
#       'name' => name,
#       'index' => index,
#       'orientation' => orientation
#     );
# }

def getMacMatches(cursor, db: str, micId: str):
  cursor.execute(
    f"""
    SELECT
      `mic_start`,
      `mic_end`,
      `mac_contig_id`,
      `mac_name`,
      ABS(`index`) AS `index`,
      `orientation`
    FROM
      `{db}`.`match`
    WHERE
      `mic_contig_id` = '{micId}'
    AND
      `is_preliminary` = 1;
    """
  )

  matches = []
  for assoc in cursor.fetchall():
    matches.append({
      'start': assoc['mic_start'],
      'end': assoc['mic_end'],
      'mac': {
          assoc['mac_name']: {
            'contig_id': assoc['mac_contig_id'],
            'name': assoc['mac_name'],
            'index': [assoc['index']],
            'orientation': [assoc['orientation']],
          }
      }
    })
  return matches

# function makeMergedMatch(start, end, start_macs, end_macs) {
#   return
#     array(
#       'start' => start,
#       'end' => end,
#       'start_macs' => start_macs,
#       'end_macs' => end_macs
#     );
# }

# function toMergedMatch(match) {
#   return
#     makeMergedMatch(
#       match['start'],
#       match['end'],
#       match['mac'],
#       match['mac']
#     );
# }

def mergeable(start_1, end_1, start_2, end_2):
  if (end_1 <= (start_2 - 2)):
     return False
  if (end_2 <= (start_1 - 2)):
     return False
  return True

def mergeMacs(mac1, mac2):
  if ((mac1['name'] != mac2['name']) or (mac1['contig_id'] != mac2['contig_id'])):
    raise Exception('MACs are not mergable')

  index_orient = (
    list(zip(mac1['index'], mac1['orientation'])) +
    list(zip(mac2['index'], mac2['orientation']))
  )
  index_orient = list(sorted(index_orient, key = lambda x: x[0]))

  return {
    'contig_id': mac1['contig_id'],
    'name': mac1['name'],
    'index': [x[0] for x in index_orient],
    'orientation': [x[1] for x in index_orient],
  }

def mergeAllMacs(macs1: dict, macs2: dict):
  allNames = set(macs1.keys()) | set(macs2.keys())
  merged = {}
  for name in allNames:
    if (name in macs1) and (name in macs2):
      merged[name] = mergeMacs(macs1[name], macs2[name])
    elif name in macs1:
      merged[name] = macs1[name]
    elif name in macs2:
      merged[name] = macs2[name]
    else:
      raise Exception("Impossible")
  return merged

# For the non-strict ies, separate out by the MAC contig
def groupMatchesByMac(matches: dict):
  groupedMatches = collections.defaultdict(list)
  for match in matches:
    name = next(iter(match['mac'].keys()))
    groupedMatches[name].append(match)
  return list(groupedMatches.values())


def mergeMatches(merged1, merged2):
  if merged1['start'] < merged2['start']:
    start = merged1['start']
    start_macs = merged1['start_macs']
  elif merged1['start'] > merged2['start']:
    start = merged2['start']
    start_macs = merged2['start_macs']
  elif merged1['start'] == merged2['start']:
    start = merged1['start']
    start_macs = mergeAllMacs(merged1['start_macs'], merged2['start_macs'])

  if merged1['end'] < merged2['end']:
    end = merged2['end']
    end_macs = merged2['end_macs']
  elif merged1['end'] > merged2['end']:
    end = merged1['end']
    end_macs = merged1['end_macs']
  elif merged1['end'] == merged2['end']:
    end = merged1['end']
    end_macs = mergeAllMacs(merged1['end_macs'], merged2['end_macs'])
  
  return {
    'start': start,
    'end': end,
    'start_macs': start_macs,
    'end_macs': end_macs,
  }

# function cmpByStart(x, y) {
#   if (x['start'] < y['start']) return -1;
#   if (x['start'] > y['start']) return 1;
#   return 0;
# }

def mergeAllMatches(matches):
  matches = list(sorted(matches, key=lambda x: x['start']))
  mergedMatches = []
  for match in matches:
    merged = {
      "start": match["start"],
      "end": match["end"],
      "start_macs": match["mac"],
      "end_macs": match["mac"],
    }

    if len(mergedMatches) > 0:
      lastMerged = mergedMatches[-1]
      if mergeable(lastMerged['start'], lastMerged['end'], merged['start'], merged['end']):
        lastMerged = mergeMatches(mergedMatches[-1], merged)
      else:
        mergedMatches.append(merged)
    else:
      mergedMatches.append(merged)
  return mergedMatches

# def makeGroupedIes(start, end, macs):
#   return array('start' => start, 'end' => end, 'macs' => macs)

def getGroupedIes(mergedMatches):
  ies = []
  for i in range(len(mergedMatches) - 1):
    intersectNames = (
      set(mergedMatches[i]['end_macs'].keys()) &
      set(mergedMatches[i + 1]['start_macs'].keys())
    )
    
    if len(intersectNames) > 0:
      intersect = []
      for name in intersectNames:
        intersect.append({
          'name': name,
          'contig_id': mergedMatches[i]['end_macs'][name]['contig_id'],
          'left_index': mergedMatches[i]['end_macs'][name]['index'],
          'right_index': mergedMatches[i + 1]['start_macs'][name]['index'],
          'left_orientation': mergedMatches[i]['end_macs'][name]['orientation'],
          'right_orientation': mergedMatches[i + 1]['start_macs'][name]['orientation'],
        })
      ies.append({
        'start': mergedMatches[i]['end'] + 1,
        'end': mergedMatches[i + 1]['start'] - 1,
        'macs': intersect
      })
  return ies

# function quoteVal(val) {
#   return "'{val}'";
# }

# // function getSqlIesValues(&values, &currIesId, micContig, groupedIes) {
# //   foreach (groupedIes as ies) {
# //     foreach (ies['macs'] as mac) {
# //       leftIndexStr = implode(',', mac['left_index']);
# //       rightIndexStr = implode(',', mac['right_index']);
# //       values[] =
# //         '(' .
# //           quoteVal(currIesId) . ',' .
# //           quoteVal(micContig['contig_id']) . ',' .
# //           quoteVal(micContig['name']) . ',' .
# //           quoteVal(ies['start']) . ',' .
# //           quoteVal(ies['end']) . ',' .
# //           quoteVal(ies['end'] - ies['start'] + 1) . ',' .
# //           quoteVal(mac['contig_id']) . ',' .
# //           quoteVal(mac['name']) . ',' .
# //           quoteVal(leftIndexStr) . ',' .
# //           quoteVal(rightIndexStr) .
# //         ')';
# //       currIesId++;
# //     }
# //   }
# // }

def getSqlIesValues(values, currIesId, micContig, groupedIes):
  for ies in groupedIes:
    for mac in ies['macs']:
      leftIndexStr = ','.join([str(x) for x in mac['left_index']])
      rightIndexStr = ','.join([str(x) for x in mac['right_index']])
      leftOrientationStr = ','.join([str(x) for x in mac['left_orientation']])
      rightOrientationStr = ','.join([str(x) for x in mac['right_orientation']])
      values.append(
        "(" +
          f"'{currIesId}'" + ',' +
          f"'{micContig['contig_id']}'" + ',' +
          f"'{micContig['name']}'" + ',' +
          f"'{ies['start']}'" + ',' +
          f"'{ies['end']}'" + ',' +
          f"'{ies['end'] - ies['start'] + 1}'" + ',' +
          f"'{mac['contig_id']}'" + ',' +
          f"'{mac['name']}'" + ',' +
          f"'{leftIndexStr}'" + ',' +
          f"'{rightIndexStr}'" + ',' +
          f"'{leftOrientationStr}'" + ',' +
          f"'{rightOrientationStr}'" +
        ")"
      )
      currIesId += 1
  return currIesId

# function test() {
#   matches = array(
#     makeMatch(0, 100, 'A'), 
#     makeMatch(0, 100, 'B'), 
#     makeMatch(0, 1000, 'C'), 
#     makeMatch(0, 1000, 'E'),
#     makeMatch(1002, 1004, 'E'),
#   );

#   printf("matches:\n");
#   var_export(matches);
#   printf("\n");

#   mergedMatches = mergeAllMatches(matches);

#   printf("mergedMatches:\n");
#   var_export(mergedMatches);
#   printf("\n");

#   groupedIes = getGroupedIes(mergedMatches);

#   printf("groupedIes:\n");
#   var_export(groupedIes);
#   printf("\n");

#   sqlIesValues = array();
#   currIesId = 0;
#   micContig = array(
#     'contig_id' => 1,
#     'name' => 'Mic_Contig_1'
#   );
#   getSqlIesValues(sqlIesValues, currIesId, micContig, groupedIes);

#   printf("sqlIesValues:\n");
#   var_export(sqlIesValues);
#   printf("\n");
# }

def make_ies_table(db: str, ies_type: str):
  common_utils.log(f"{db} {ies_type}")

  if ies_type not in ["strict", "weak"]:
    raise Exception("Impossible.") 

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)

  micContigs = getMicContigs(cursor, db)
  comment = {
    "strict": (
      "segments on the precursor contigs that lie between 2 product matches" +
      " that does not overlap an other matches from other contigs"
    ),
    "weak": "segments on the precursor contigs that lie between 2 product matches",
  }

  cursor.execute(f"DROP TABLE IF EXISTS `{db}`.`ies_{ies_type}`;")
  cursor.execute(
    f"""
    CREATE TABLE `{db}`.`ies_{ies_type}` (
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
    ) COMMENT '{comment[ies_type]}';
    """
  )
  values = []
  currIesId = 1
  for i in range(len(micContigs)):
    micContig = micContigs[i]
    matches = getMacMatches(cursor, db, micContig["contig_id"])

    if ies_type == "weak":
      # separate out by MAC and then calculate IES
      matchesByMac = groupMatchesByMac(matches)
      groupedIes = []
      for singleMacMatches in matchesByMac:
        singleMacMergedMatches = mergeAllMatches(singleMacMatches)
        singleMacGroupedIes = getGroupedIes(singleMacMergedMatches)
        for ies in singleMacGroupedIes:
          groupedIes.append(ies)
    elif ies_type == "strict":
      # pool all MAC matches together so that there are no overlapping IES/MDS between different MAC contigs
      mergedMatches = mergeAllMatches(matches)
      groupedIes = getGroupedIes(mergedMatches)
    else:
      raise Exception("Impossible.")

    currIesId = getSqlIesValues(values, currIesId, micContig, groupedIes);
    
    # every 100 contigs write to SQL
    if ((i + 1) % 100 == 0) or (i == len(micContigs)):
      common_utils.log(f"{i} / {len(micContigs)}")
      insertQuery = f"INSERT INTO `{db}`.`ies_{ies_type}` VALUES " + ', '.join(values) + ";"
      values = []
      cursor.execute(insertQuery)

  cursor.close()
  conn.close()

if __name__ == "__main__":
  make_ies_table("hello_world", "strict")

# TEST THE WEAK AND CLEEEEEAAAAAAN THIS STUFF UP!!!