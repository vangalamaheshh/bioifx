#!/usr/bin/env python

#vim: syntax=python tabstop=2 expandtab

from random import randint
import sys
import math

def generateNums(num, numRange):
  numList = list()
  count = 0
  while count < num:
    randNum = randint(1, numRange)
    if randNum not in numList:
      count = count + 1
      numList.append(randNum)
  
  return numList

def generateMM():
  numList = generateNums(5, 70)
  numList.extend(generateNums(1, 25))
  return numList

def generatePB():
  numList = generateNums(5, 69)
  numList.extend(generateNums(1, 26))
  return numList

if __name__ == "__main__":
  numList = generatePB() if len(sys.argv) > 1 and sys.argv[1] == "PB" else generateMM()
  print(", ".join([str(num) for num in numList[:-1]]), " - ", str(numList[-1]))
