import pandas as pd
import datetime
import csv

def log(x):
  print(datetime.datetime.now().strftime("%H:%M:%S: ") + str(x))

def read_tsv(file: str, **args) -> pd.DataFrame:
  return pd.read_csv(file, sep="\t", **args)

def write_tsv(data: pd.DataFrame, file: str, **args):
  return data.to_csv(file, sep="\t", index=False, quoting=csv.QUOTE_NONNUMERIC, **args)

def get_contig_info(cursor, db: str):
  cursor.execute(
    f"""
    SELECT
      `{to_db}`.`contig_id`,
      `nucleus` AS `contig_nucleus`,
      `name` AS `contig_name`
    FROM `contig`;
    """
  )
  return cursor.fetchall()
