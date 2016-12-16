#!/usr/bin/env Rscript

#--------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: Sep, 22, 2016
#--------------------------------

stat_list <- list.files(pattern = "*.stats.txt", recursive = T)
df <- read.table(stat_list[1], sep = ":", header = F, colClasses = c("character","numeric"), row.names = 1)
for(stat_file in stat_list[2:length(stat_list)]){
  df_temp <- read.table(stat_file, sep = ":", header = F, colClasses = c("character","numeric"), row.names = 1)
  df[, stat_file] <- df_temp[,1]
}
names(df) <- stat_list
write.csv(file = "matrix_file.csv", df, quote = F)