
import mysql.connector
import collections

import pandas as pd

import constants
import common_utils
import mysql_utils
import interval_utils


IES_COLUMNS = [
  "mic_contig_id",
  "mic_name",
  "mic_start",
  "mic_end",
  "length",
  "mac_contig_id",
  "mac_name",
  "left_index",
  "right_index",
  "left_orientation",
  "right_orientation",
]

BATCH_SIZE_CONTIGS = 100

def get_mic_contigs(
  cursor,
  db: str,
):
  cursor.execute(
    f"""
    SELECT
      `contig_id`,
      `name`
    FROM `{db}`.`contig`
    WHERE `nucleus` = 'mic';
    """
  )
  return pd.DataFrame.from_records(
    cursor.fetchall(),
    columns = ["contig_id", "name"]
  )

def get_mac_matches(
  cursor,
  db: str,
  mic_contig_id: str,
):
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
      `mic_contig_id` = '{mic_contig_id}'
    AND
      `is_preliminary` = 1;
    """
  )

  return pd.DataFrame.from_records(
    cursor.fetchall(),
    columns = [
      "mic_start",
      "mic_end",
      "mac_contig_id",
      "mac_name",
      "index",
      "orientation"
    ],
  )

def get_mic_ies(
  cursor,
  db: str,
  ies_type: str,
  mic_name: str,
  mic_contig_id: str,
):
  match_data = get_mac_matches(cursor, db, mic_contig_id)

  if ies_type == "strict":
    # matches from all MAC contigs are pooled together
    match_data_list = [match_data]
  elif ies_type == "weak":
    # matches from different MAC contigs are separated
    match_data_list = [data for _, data in match_data.groupby('mac_contig_id')]
  else:
    raise Exception(f"Unknown IES type: {ies_type}")

  ies_data_list = []
  for match_data in match_data_list:
    # Get the match intervals
    intervals = match_data[["mic_start", "mic_end"]]
    intervals = intervals.rename({"mic_start": "start", "mic_end": "end"}, axis="columns")
    intervals = intervals.to_dict("records")
    intervals = interval_utils.get_union(intervals)

    # Get the IESs as the intervals between the unions of the matches
    intervals_ies = []
    for i in range(len(intervals) - 1):
      intervals_ies.append({"start": intervals[i]["end"] + 1, "end": intervals[i + 1]["start"] - 1})
    intervals_ies = pd.DataFrame.from_records(intervals_ies, columns=["start", "end"])

    # In the rare case that a MAC contigs has multiple matches at a single
    # MIC locus, they should be merged into a single row
    def collapse_matches(data: pd.DataFrame):
      data = data.sort_values("index")
      return pd.Series({
        "index": ",".join([str(i) for i in data["index"]]),
        "orientation": ",".join(data["orientation"])
      })

    # Get the matches which flank an IES on the left
    matches_left = match_data.drop("mic_start", axis="columns")
    matches_left = matches_left.rename({"mic_end": "start"}, axis="columns")
    matches_left["start"] += 1
    matches_left = pd.merge(intervals_ies, matches_left, on="start")
    matches_left = matches_left.groupby(["start", "end", "mac_contig_id", "mac_name"])
    matches_left = matches_left.apply(collapse_matches)

    # Get the matches which flank an IES on the right
    matches_right = match_data.drop("mic_end", axis="columns")
    matches_right = matches_right.rename({"mic_start": "end"}, axis="columns")
    matches_right["end"] -= 1
    matches_right = pd.merge(intervals_ies, matches_right, on="end")
    matches_right = matches_right.groupby(["start", "end", "mac_contig_id", "mac_name"])
    matches_right = matches_right.apply(collapse_matches)

    # Get the IESs which are flanked on both the right and left
    # by a MAC match from the same MAC contig.
    ies_data = pd.merge(
      matches_left,
      matches_right,
      how = "inner",
      left_index = True,
      right_index = True,
      suffixes = ("_left", "_right")
    )
    ies_data = ies_data.reset_index()
    ies_data_list.append(ies_data)
  
  if len(ies_data_list) > 0:
    ies_data = pd.concat(ies_data_list, axis="index")
    ies_data["length"] = ies_data["end"] - ies_data["start"] + 1
    ies_data["mic_name"] = mic_name
    ies_data["mic_contig_id"] = mic_contig_id

    ies_data = ies_data.rename(
      {
        "start": "mic_start",
        "end": "mic_end",
        "index_left": "left_index",
        "index_right": "right_index",
        "orientation_left": "left_orientation",
        "orientation_right": "right_orientation",
      },
      axis = "columns",
    )
    ies_data = ies_data[IES_COLUMNS]
  else:
    ies_data = pd.DataFrame([], columns=IES_COLUMNS)
  return ies_data


def create_ies_table(db: str, ies_type: str):
  common_utils.log(f"create_ies_table {db} {ies_type}")

  if ies_type not in ["strict", "weak"]:
    raise Exception("Impossible.") 

  conn = mysql_utils.get_connection()
  cursor = conn.cursor()

  mic_contig_list = get_mic_contigs(cursor, db)

  comment = {
    "strict": (
      "segments on the precursor contigs that lie between 2 product matches" +
      " that does not overlap an other matches from other contigs"
    ),
    "weak": "segments on the precursor contigs that lie between 2 product matches",
  }

  cursor = conn.cursor()
  cursor.execute(f"DROP TABLE IF EXISTS `{db}`.`ies_{ies_type}`;")
  cursor.execute(
    f"""
    CREATE TABLE `{db}`.`ies_{ies_type}` (
      `ies_id` int NOT NULL AUTO_INCREMENT COMMENT 'primary key for the table',
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
  sql_columns = ", ".join(f"`{x}`" for x in IES_COLUMNS)
  ies_data_list = []
  for mic_num, mic_contig in enumerate(mic_contig_list.to_records(), 1):
    ies_data_list.append(
      get_mic_ies(
        cursor,
        db,
        ies_type,
        mic_contig["name"],
        mic_contig["contig_id"],
      )
    )
    
    if ((mic_num % BATCH_SIZE_CONTIGS) == 0) or (mic_num == len(mic_contig_list)):
      common_utils.log(f"{mic_num} / {len(mic_contig_list)}")
      ies_data = pd.concat(ies_data_list, axis="index")
      ies_data_list = []
      values = mysql_utils.make_sql_values(ies_data[IES_COLUMNS])
      cursor.execute(
        f"""
        INSERT INTO `{db}`.`ies_{ies_type}`
        ({sql_columns})
        VALUES {values};
        """
      )
      cursor.close()
      conn.close()
      conn = mysql_utils.get_connection()
      cursor = conn.cursor()

  cursor.close()
  conn.close()

# make_ies_table("hello_world", "strict")
# make_ies_table("hello_world", "weak")

# conn = mysql_utils.get_connection()
# ies = get_mic_ies_strict(conn, "hello_world", "OXYTRI_MIC_67642", 26284)
# print("hi")
  # mic_name = "OXYTRI_MIC_67642"
  # mic_contig_id = 26284