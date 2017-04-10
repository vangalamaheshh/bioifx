#!/usr/bin/env python
#vim: syntax=python tabstop=2 expandtab

__author__ = "Mahesh Vangala"
__email__ = "<vangalamaheshh@gmail.com>"
__date__ = "Apr 10, 2017"
__doc__ = """
  Print Variant Annotation into CSV

  Just a formatting script
"""

import pandas as pd
import numpy as np

def formatDB(dbTemp):
  db = pd.DataFrame()
  db[["id", "nuc_change"]] = dbTemp[["query", "_id"]][dbTemp["notfound"] != True]
  db["aa_change"] = [dbTemp["dbnsfp.aa.ref"][index] + '>' + dbTemp["dbnsfp.aa.alt"][index] + dbTemp["dbnsfp.aa.pos"][index] \
                      if pd.notnull(dbTemp["dbnsfp.aa.ref"][index]) else getAminoChange(dbTemp["dbnsfp.aa"][index]) \
                      for index in dbTemp[dbTemp["notfound"] != True].index]
  db["gene"] = [dbTemp["clinvar.gene.symbol"][index] if pd.notnull(dbTemp["clinvar.gene.symbol"][index]) \
                else '-' for index in dbTemp[dbTemp["notfound"] != True].index]
  db["desc"] = [dbTemp["dbnsfp.clinvar.trait"][index] if pd.notnull(dbTemp["dbnsfp.clinvar.trait"][index]) \
                else '-' for index in dbTemp[dbTemp["notfound"] != True].index]
  return db

def getAminoChange(val):
  try:
    val = eval(val)
  except NameError:
    return "-"
  except TypeError:
    return "-"
  return_val = []
  for elem in val:
    if all(k in elem for k in ["ref", "alt", "pos"]):
      if isinstance(elem["pos"], int):
        vals = list()
        vals.append(str(elem["pos"]))
      else:
        vals = elem["pos"]
      
      return_val.append(elem["ref"] + ">" + elem["alt"] + repr(vals))
  
  return ";".join(return_val)

if __name__ == "__main__":
  db = pd.read_csv("variants.csv", sep = ",", header = 0)
  db = formatDB(db)
  db.to_csv("variants.format.csv", header = True, index = False)
