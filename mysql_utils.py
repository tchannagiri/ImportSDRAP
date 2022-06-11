import pandas as pd

import common_utils
import mysql.connector

def get_connection() -> mysql.connector.MySQLConnection:
  conn = mysql.connector.connect(
    host = "localhost",
    port = "8888",
    user = "web",
    password = "RZhRwsau6HZrMUXf",
    autocommit = True,
  )
  return conn


def make_sql_values(data):
  data = data.applymap(lambda x: "NULL" if pd.isna(x) else f"'{x}'")
  data = data.apply(
    lambda row: "(" + ",".join(row) + ")",
    axis = "columns"
  )
  return ",".join(data)

def upload_in_chunks(
  data: pd.DataFrame,
  columns: list[str],
  cursor,
  table_name: str,
  batch_size: int,
):
  columns_sql = ", ".join([f"`{x}`" for x in columns])
  for i in range(0, data.shape[0], batch_size):
    common_utils.log(f"{i} / {data.shape[0]}")
    values = make_sql_values(data.loc[i : (i + batch_size - 1), columns])
    cursor.execute(
      f"""
      INSERT INTO {table_name}
      ({columns_sql})
      VALUES {values}
      """
    )