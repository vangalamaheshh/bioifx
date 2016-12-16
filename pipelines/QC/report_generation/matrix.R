# vim: syntax=r tabstop=4 expandtab
# coding: utf-8

#!/bin/env R

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

make_option( c("-h", "--help"), action="store_true", default=FALSE, help="print this help message and exit" )

option_list = list(
    make_option(c("-f", "--file"), action="store", default=NA, type='character', help="read distribution output file. Multiple files can be given as comma separated names.")
)



opt = parse_args(OptionParser(option_list=option_list))

opt$f
