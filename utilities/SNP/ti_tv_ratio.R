#!/usr/bin/env Rscript
# vim: syntax=r tabstop=2 expandtab

library(ggplot2)
library(scales)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

data <- read.csv(args[1], header = T)
chrom_names <- unique(data[, 1])

ggimages <- list()

for (chrom_name in chrom_names) {
  
  sub_data <- data[data[, 1] == chrom_name, ]
  cur_plot <- ggplot(sub_data, aes(x=window_start, y=titv)) +
    geom_point() +
    stat_smooth() +
    scale_x_continuous(labels=comma) +
    xlab("Genomic Position") +
    ylab("Ti/Tv") +
    ggtitle(paste("Ti/Tv by 100,000 base pair windows - ", chrom_name, sep = ""))
  ggimages <- c(ggimages, list(cur_plot))
}

png(args[2], width = 8, height = 8, unit = "in", res = 300)
do.call("grid.arrange", c(ggimages, ncol = 1))
dev.off()
