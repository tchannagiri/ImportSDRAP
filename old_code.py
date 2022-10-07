# # #track name="prec_eliminated_sequences" description="intervals comprising complement of the precursor segments labelled comp_[left-flanking-segment-prod-id]_[left-flanking-segment-index]_[right-flanking-segment-prod-id]_[right-flanking-segment-index]" itemRgb-"On"
# # OXYTRI_MIC_33550	1	8703	none_0_Contig4386.0_-1	0	+	1	8703	255,0,153
# # sdrap_oxy_mac2012_100720_prec_eliminated_sequences.bed
# def make_ie_bed(to_db: str, from_db: str, ies_type: str):
#   handle = open(f"../output/{from_db}_prec_{ies_type}_ies_sequences.bed", "w")

#   # write header
#   if ies_type == "strict":
#     handle.write(
#       '#track name="prec_strict_ies_sequences" ' +
#       'description="Strict IESs. Intervals on the precursor segment between two consecutive matches of a product segment ' +
#       'that do not overlap matches of any other product segment. ' +
#       'Labelled [product-name]_[left-flanking-match-indexes]_[right-flanking-match-indexes]" ' +
#       'itemRgb-"On"\n'
#     )
#   elif ies_type == "weak":
#     handle.write(
#       '#track name="prec_weak_ies_sequences" ' +
#       'description="Weak IESs. Intervals on the precursor segment between two consecutive matches of a product segment.' +
#       'Labelled [product-name]_[left-flanking-match-indexes]_[right-flanking-match-indexes]" ' +
#       'itemRgb-"On"\n'
#     )
#   else:
#     raise Exception("Unknown ies type: " + str(ies_type))

#   conn = mysql_utils.get_connection()
#   cursor = conn.cursor(dictionary=True)
#   cursor.execute("SELECT COUNT(*) AS `count` FROM `{to_db}`.`ies_{$type}`;")
#   num_ies = cursor.fetchall()["count"]
  
#   batch_size_rows = 10000
#   for i in range(0, num_ies, batch_size_rows):
#     cursor.execute(
#       f"""
#       SELECT
#         `mic_name`,
#         `mic_start`,
#         `mic_end`,
#         `mac_name`,
#         `left_index`,
#         `right_index`
#       FROM `{to_db}`.`ies_{ies_type}`
#       LIMIT {i}, {batch_size_rows};
#       """
#     )

#     for result in cursor.fetchall():
#       handle.write(
#         f"{result['mic_name']}\t{result['mic_start']}\t{result['mic_end']}\t" + # chrom chromStart chromEnd
#         f"{result['mac_name']}_{result['left_index']}_{result['right_index']}\t" + # name
#         f"0\t+\t{result['mic_start']}\t{result['mic_end']}\t255,0,153\n" # score strand thickStart thickEnd itemRGB
#       )

#   cursor.close()
#   conn.close()



def create_alias_table(
  to_db: str,
  from_db: str,
  alias_file_list: list[str],
):
  common_utils.log(f"create_alias {to_db} {from_db} {alias_file_list}")

  conn = mysql_utils.get_connection()
  cursor = conn.cursor(dictionary=True)
  cursor.execute(f"SELECT * FROM `{from_db}`.`alias`;")
  alias_from = cursor.fetchall()
  alias_from = pd.DataFrame.from_records(alias_from)

  cursor.execute(f"DROP TABLE IF EXISTS `{to_db}`.`alias`;")
  cursor.execute(
    f"""
    CREATE TABLE `{to_db}`.`alias` (
      `alias_id` INT NOT NULL AUTO_INCREMENT COMMENT 'primary key for the table',
      `id` INT NOT NULL COMMENT 'primary key of the table for the corresponding sequence',
      `name` VARCHAR(50) NOT NULL COMMENT 'the primary name of the sequence',
      `alias` VARCHAR(50) NOT NULL COMMENT 'alias of the sequence',
      `table` VARCHAR(50) NOT NULL COMMENT 'the table that `id` refers to',
      `type` VARCHAR(50) DEFAULT NULL COMMENT 'the type of alias',
      PRIMARY KEY (`alias_id`),
      KEY `name` (`name`),
      KEY `alias` (`alias`)
    ) COMMENT='table of sequence aliases';
    """
  )

  for alias_file in alias_file_list:
    common_utils.log(alias_file)
    alias_to = pd.read_csv(alias_file)
    alias_to = alias_to.merge(
      alias_from[["nuc_id", "alias"]],
      left_on = "primary",
      right_on = "alias",
      how = "inner",
    )
    alias_to = alias_to.drop("alias", axis="columns")
    alias_to = alias_to.melt(
      id_vars = ["nuc_id", "primary"],
      var_name = "type",
      value_name = "alias",
    )
    alias_to = alias_to.dropna(axis='index')
    alias_to = alias_to.rename({"nuc_id": "id", "primary": "name"}, axis="columns")
    alias_to["table"] = "contig"
    mysql_utils.upload_in_chunks(
      alias_to,
      ["id", "name", "alias", "table", "type"],
      cursor,
      to_db,
      "alias",
      constants.SQL_BATCH_UPLOAD_ROWS,
    )
    # for i in range(0, alias_to.shape[0], SQL_BATCH_UPLOAD_ROWS):
    #   common_utils.log(f"{i} / {alias_to.shape[0]}")
    #   values = mysql_utils.make_sql_values(
    #     alias_to.loc[
    #       i : (i + SQL_BATCH_UPLOAD_ROWS - 1),
    #       ["id", "name", "alias", "table", "type"]
    #     ])
    #   cursor.execute(
    #     f"""
    #     INSERT INTO `{to_db}`.`alias`
    #     (
    #       `id`,
    #       `name`,
    #       `alias`,
    #       `table`,
    #       `type`
    #     )
    #     VALUES {values};
    #     """
    #   )
  # cursor.close()
  # conn.close()


# Problems to fix:
# 1. Remove "CDS" features which seem to be duplicate of "exon".
#    Only keep "gene", "mRNA", "exon" features.
# 2. The names used are genbank "LASUXXX" format, so convert them to primary
#    format "ContigXXX".
# def fix_tetsp_mac2015(data):
#   data = data.loc[data["type"].isin(["gene", "mRNA", "exon"])].reset_index(drop=True).copy()
#   alias = pd.read_csv(constants.TETSP_MAC_2015_ALIAS, sep="\t").set_index("genbank")

#   data["contig_name_new"] = alias.loc[data["contig_name"], "primary"].reset_index(drop=True)
#   data = data.to_dict("records")
#   for row in data:
#     contig_name_old = row["contig_name"]
#     contig_name_new = row["contig_name_new"]
#     row["attr_id"] = row["attr_id"].replace(row["contig_name"], row["contig_name_new"])
#     for key in ["contig_name", "attr_id", "attr_name", "attr_parent"]:
#       if row[key] is not None:
#         row[key] = row[key].replace(contig_name_old, contig_name_new)

#   data = pd.DataFrame.from_records(data)
#   data = data.drop("contig_name_new", axis="columns")
#   return data
