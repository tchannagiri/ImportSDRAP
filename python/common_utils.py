import os
import pandas as pd
import datetime
import csv

def log(x):
  print(datetime.datetime.now().strftime("%H:%M:%S: ") + str(x))

def make_parent_dir(file_name):
  os.makedirs(os.path.dirname(file_name), exist_ok=True)

def read_tsv(file: str, **args) -> pd.DataFrame:
  return pd.read_csv(file, sep="\t", **args)

def write_tsv(data: pd.DataFrame, file: str, **args):
  make_parent_dir(file)
  return data.to_csv(file, sep="\t", index=False, quoting=csv.QUOTE_ALL, **args)

def get_contig_info(cursor, db: str):
  cursor.execute(
    f"""
    SELECT
      `{db}`.`contig_id`,
      `nucleus` AS `contig_nucleus`,
      `name` AS `contig_name`
    FROM `contig`;
    """
  )
  return cursor.fetchall()
