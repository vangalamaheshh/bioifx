#!/usr/bin/env python

import sys
import re
from collections import OrderedDict
from ruamel.yaml import YAML

def get_icd9_info(stroke_file):
  data = open(stroke_file, 'r', encoding = 'utf-8').read()
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
  data = open(stroke_file, 'r', encoding = 'utf-8').read()
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
      if line and line[0] == '#':
        continue
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

def print_yaml_with_comments(info):
	yaml = YAML()
	input = """
inclusion:
  icd9s:
    -
  icd10s:
    -
"""
	data = yaml.load(input)
	for icd_type in ['icd9s', 'icd10s']:
		start = True
		index = 0
		for key, val in info[icd_type].items():
			key = key.replace(' ', '')
			key = key.replace('-', ' - ')
			if start:
				start = False
				data['inclusion'][icd_type][0] = key
			else:
				data['inclusion'][icd_type].append(key)
			data['inclusion'][icd_type].yaml_add_eol_comment(val, index)
			index += 1
	yaml.dump(data, sys.stdout)
			
		
  
if __name__ == '__main__':
  stroke_file = sys.argv[1]
  icd9_info = get_icd9_info(stroke_file)
  icd10_info = get_icd10_info(stroke_file)
  print_yaml_with_comments(OrderedDict({"icd9s": icd9_info, "icd10s": icd10_info}))
  
