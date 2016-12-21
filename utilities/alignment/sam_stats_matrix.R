#!/usr/bin/env Rscript
# vim: syntax=r tabstop=4 expandtab

#----------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Aug, 1, 2015
#-----------------------------------

library(ggplot2)
library(reshape2)

args <- commandArgs( trailingOnly = TRUE )

data <- read.csv( args[1], sep=",", header=TRUE, check.names=F, row.names = 1 )

x <- data
x["id"] <- rownames(data)
x_melt <- melt(x)

png( args[2], width = 8, height = 8, unit="in",res=300 )

upper_limit <- max(x$TotalReadCount)
limits <- seq( 0, upper_limit, length.out=10)

cust_labels <- vector("character",length=length(limits))

if( nchar(upper_limit) < 7 ) {
  cust_labels <- paste(round(limits/1000),"K",sep="") 
  limits <- round(limits/1000) * 1000
} else {
  cust_labels <- paste(round(limits/1000000),"M",sep="") 
  limits <- round(limits/1000000) * 1000000
}


colors <- c(TotalReadCount="Grey", UniqueReadCount="steelblue")

ggplot(x_melt, aes(x=id, y=value, fill=variable)) + 
    geom_bar( stat = "identity", position="identity") + 
    scale_y_continuous("",limits=c(0,upper_limit), labels=cust_labels, breaks=limits) + 
    scale_fill_manual(values=colors) + 
    labs( title="Read Alignment Report\n\n", x = "Sample Names", y="") + 
    guides(fill=guide_legend(title=NULL)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5, size=10))

junk <- dev.off()

