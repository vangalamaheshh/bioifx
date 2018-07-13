#!/usr/bin/env python

import sys
import re
from collections import OrderedDict

def get_icd9_info(stroke_file):
  data = open(stroke_file, 'r').read()
  data = data.split('\n')
  d = OrderedDict()
  trigger_on = False
  for line in data:
    if line[:8] == 'ICD-9-DM':
      trigger_on = True
      continue 
    if line[:9] == 'ICD-10-CM':
      trigger_on = False
      break
    if trigger_on:
      match_obj = re.search(r'(\d+.{0,1}\d*)\s+(.+)$', line)
      if match_obj and len(match_obj.groups()) == 2:
        d[match_obj.group(1)] = match_obj.group(2)
  return d

def get_icd10_info(stroke_file):
  data = open(stroke_file, 'r').read()
  data = data.split('\n')
  d = OrderedDict()
  trigger_on = False
  for line in data:
    if line[:9] == 'ICD-10-CM':
      trigger_on = True 
      continue
    if line[:25] == 'Note for Coverdell users:': 
      trigger_on = False
      break
    if trigger_on:
      match_obj = re.search(r'([IGO]\d+.{0,1}\d{0,}\s{0,}\-{0,1}\s{0,}[IGO]{0,1}\d{0,}.{0,1}\d{0,})\s+(.+)\s+(I\d+)\s+(.+)', line)
      if match_obj:
        if len(match_obj.groups()) == 4:
          (icd10_code1, icd10_text1, icd10_code2, icd10_text2) = match_obj.groups()
          d[icd10_code1] = icd10_text1
          d[icd10_code2] = icd10_text2
        else:
          raise "Error in regexp logic for ICD10 code fetch."
      else:
        match_obj = re.search(r'([IGO]\d+.{0,1}\d{0,}\s{0,}\-{0,1}\s{0,}[IGO]{0,1}\d{0,}.{0,1}\d{0,})\s+(.+)', line)
        if match_obj:
          if len(match_obj.groups()) == 2:
            d[match_obj.group(1)] = match_obj.group(2)
          else:
            raise "Error in regexp logic for ICD10 code fetch."
  return d

if __name__ == '__main__':
  stroke_file = sys.argv[1]
  icd9_info = get_icd9_info(stroke_file)
  icd10_info = get_icd10_info(stroke_file)
  for key, val in icd10_info.items():
    print(key, "\t", val)
