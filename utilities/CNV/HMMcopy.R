#!/usr/bin/env Rscript

# vim: syntax=r tabstop=2 expandtab

#------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Oct, 4, 2016
#------------------------------------

library(HMMcopy)

run_HMMcopy <- function(r_file, g_file, m_file, qc_out_png, cnv_out_pdf) {
  normal_reads <- wigsToRangedData(r_file, g_file, m_file)
  normal_copy <- correctReadcount(normal_reads)
  # GC bias plot
  png(qc_out_png, height = 8, width = 8, unit = "in", res = 300)
  par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))
  plotBias(normal_copy, pch = 20, cex = 0.5)
  dev.off()
  # pdf with corrected read counts 
  # one chromosome per page
  pdf(cnv_out_pdf)
  for (chr in c(seq(1, 22), "X")) {
    par(mar = c(4, 4, 2, 0)) 
    plotCorrection(normal_copy, pch = ".", chr = paste("chr", chr, sep = "")) 
  }
  dev.off()
}

#------- COMMANDLINE ARGS --------------#
args <- commandArgs(trailingOnly = TRUE)
read_file <- args[1]
gc_file <- args[2]
map_file <- args[3]
corr_map_gc_png <- args[4]
cnv_out_pdf <- args[5]
#---------------------------------------#

run_HMMcopy(read_file, gc_file, map_file, corr_map_gc_png, cnv_out_pdf)
