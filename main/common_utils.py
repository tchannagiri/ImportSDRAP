import pandas as pd
import datetime

def log(x):
  print(datetime.datetime.now().strftime("%H:%M:%S: ") + str(x))

def read_tsv(file: str, **args):
  return pd.read_csv(file, sep="\t", **args)
