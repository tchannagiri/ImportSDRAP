import pandas as pd

import common_utils
import mysql.connector

# This is the default MYSQL port for running on the server.
# If running on a local machine with ports forwarded to the main
# server, the appropriate port must be passed on the command line.
# If the host, user, or password must be changed they should be changed
# in the get_connection() function below.
PORT = "3306"

def set_port(new_port: str):
  global PORT
  PORT = new_port

def get_connection() -> mysql.connector.MySQLConnection:
  conn = mysql.connector.connect(
    host = "localhost",
    port = PORT,
    user = "web",
    password = "RZhRwsau6HZrMUXf",
    autocommit = True,
  )
  return conn

def make_sql_values(data) -> str:
  data = data.applymap(lambda x: "NULL" if pd.isna(x) else "'" + str(x).replace("'", "''") + "'")
  data = data.apply(
    lambda row: "(" + ",".join(row) + ")",
    axis = "columns"
  )
  return ",".join(data)

def upload_in_chunks(
  data: pd.DataFrame,
  columns: list[str],
  cursor,
  db: str,
  table_name: str,
  batch_size: int,
):
  columns_sql = ", ".join([f"`{x}`" for x in columns])
  for i in range(0, data.shape[0], batch_size):
    common_utils.log(f"{i} / {data.shape[0]}")
    values = make_sql_values(data.loc[i : (i + batch_size - 1), columns])
    cursor.execute(
      f"""
      INSERT INTO `{db}`.`{table_name}`
      ({columns_sql})
      VALUES {values};
      """
    )

def download_in_chunks(
  db: str,
  table_name: str,
  columns: list[str],
  cursor, # expects a non-dictionary cursor
  num_rows: int,
  batch_size: int,
):
  sql_columns = ",".join(f"`{x}`" for x in columns)
  row_list = []
  for start_row in range(0, num_rows, batch_size):
    common_utils.log(f"{start_row} / {num_rows}")
    cursor.execute(
      f"""
      SELECT {sql_columns}
      FROM `{db}`.`{table_name}`
      LIMIT {start_row}, {batch_size};
      """
    )
    row_list += cursor.fetchall()
  
  return pd.DataFrame.from_records(row_list, columns=columns)
