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


def create_database(
  conn : mysql.connector.MySQLConnection,
  to_db : str,
):
  cursor = conn.cursor()
  cursor.execute(f"DROP DATABASE IF EXISTS `{to_db}`;")
  cursor.execute(f"CREATE DATABASE `{to_db}`;")
  log(f"Create `{to_db}`")

def make_sql_values(data):
  data = data.applymap(lambda x: "NULL" if x is None else f"'{x}'")
  data = data.apply(
    lambda row: "(" + ",".join(row) + ")",
    axis = "columns"
  )
  return ",".join(data)
