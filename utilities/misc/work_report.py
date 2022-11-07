#!/usr/bin/env python

#------------------------
# @author: Mahesh Vangala
# @Date: Nov, 06, 2022
#------------------------

import sys
from collections import defaultdict
import re
import pandas as pd

def get_info(file):
  info = defaultdict(lambda: defaultdict(int))
  days = 0
  with open(file, "r") as fh:
    for line in fh:
      line = line.strip()
      if not line:
        continue
      day = re.split(r",\s*", line)[0]
      if len(day) == 3:
        days += 1
        for line in fh:
          if line[0] == '_':
            break
          elif line[0].isdigit():
            _, project, effort = re.search(r"(\d+).\s+(.+)\:\s*(\d+)", line).groups()
            info[project]['effort'] += effort
            info[project]['days'] += 1
  
  for proj in info.keys():
    info[proj]['effort%'] = f"{info[proj]['effort']/days:.2f}" 
  
  return info

if __name__ == "__main__":
  infile = sys.argv[1]
  info = get_info(infile)
  df = pd.DataFrame.from_dict(info, orient = 'index')
  df.to_csv("work-summary.csv", sep = ",", header = True, index = False)
