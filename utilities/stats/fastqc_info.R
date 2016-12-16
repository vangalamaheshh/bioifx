#!/usr/bin/env Rscript

# vim: syntax=r tabstop=4 expandtab

library("ggplot2")
library("pheatmap")

options(error = function() traceback(2))
args <- commandArgs( trailingOnly = TRUE )
fastqc_csv <- args[1]
fastqc_png <- args[2]

qc <- read.csv (fastqc_csv,header = TRUE, row.names=1,sep=",")

x <- data.frame(Sample = colnames(qc), Reads = as.numeric(as.matrix((qc["Total_Sequences",]))))



upper_limit <- max(x$Reads)

limits <- seq( 0, upper_limit, length.out=10)

cust_labels <- vector("character",length=length(limits))

if( nchar(upper_limit) < 7 ) {

  cust_labels <- paste(round(limits/1000),"K",sep="")

  limits <- round(limits/1000) * 1000

} else {

  cust_labels <- paste(round(limits/1000000),"M",sep="")

  limits <- round(limits/1000000) * 1000000

}



png(args[2],width = 12, height = 8, unit='in',res=300)

ggplot(x, aes(Sample, Reads,)) + 

  geom_bar(aes(Sample),stat = "identity", position="identity",fill="midnightblue") + 

  scale_y_continuous("",limits=c(0,upper_limit), labels=cust_labels, breaks=limits) + 

  labs( title="\nDemultiplexed Reads\n\n", x = "", y=" Reads") + 

  scale_fill_gradient()+

  theme_bw() + 

  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12),

        axis.text.y = element_text(face="bold",size=14),

        panel.border = element_rect(colour = 'black', fill=NA, size=1),

        plot.title = element_text(lineheight=.8, size = 20, face=c("italic","bold")),

        axis.text.y = element_text(""),

        axis.title = element_text(color="black",size = 18, face = "bold")

  )

dev.off()


png(args[3], width = 12, height = 8, unit='in',res=300)

vqc <- qc[!rownames(qc) %in% c("Total_Sequences","Total_Deduplicated_Percentage"),]
vqc <- as.matrix(apply(vqc, 1, function(x) ifelse (x == 'PASS',0, ifelse (x=='WARN',1,2))))

pheatmap(vqc, legend = FALSE, cluster_rows = F, cluster_cols = F, main = "\nFASTQC Summary",
        display_numbers = F, fontsize = 9, border_color = 'black', 
        color = colorRampPalette(c("lightgreen", "lightyellow", "lightpink"))(100))

dev.off()
