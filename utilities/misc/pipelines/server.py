#!/usr/bin/env python
#vim: syntax=python tabstop=2 expandtab

#---------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 27, 2018
#---------------------------

from flask import Flask, request

import json
import os
from datetime import datetime

#------  GLOBALS  --------#

#------------------------#
#    Flask methods       #
#------------------------#
app = Flask(__name__)

@app.route("/", methods = ['GET', 'POST'])
def hello():
  return "Hello World!"

