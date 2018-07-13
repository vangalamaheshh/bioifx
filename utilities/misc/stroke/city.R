#!/usr/bin/env Rscript

library("ggplot2")
library("ggrepel")
library("plyr")

get_processed_df <- function(df) {
  names(df) <- c("CITY", "COUNT")
  df$CITY <- toupper(df$CITY)
  df$PERCENT <- (df$COUNT/sum(df$COUNT))*100
  df <- df[df$PERCENT > 10, ]
  df_other <- data.frame(
    CITY = c("OTHER"), 
    COUNT = c(sum(df[df$PERCENT <= 10, "SUM"])),
    PERCENT = c(sum(df[df$PERCENT <= 10, "PERCENT"]))
  )
  df <- rbind(df, df_other)
  rownames(df) <- NULL
  print(df)
  df <- df[order(df$PERCENT), ]
  df$ymax <- cumsum(df$PERCENT)
  df$ymin <- c(0, head(df$ymax, n = -1))
  return (df)
}

plot_donut <- function(df) {
  ggplot(df, aes(fill = CITY, x = (ymax + ymin)/2, y = ymax)) +
  geom_rect(
    aes(xmin = 30, xmax = 40, ymin = ymin, ymax = ymax, fill = CITY),
    color = "White"
  ) +
  geom_label_repel(
    aes(label = paste(CITY, paste("[", round(PERCENT, digits = 3), "%]", sep = ""))),
    arrow = NULL,
    segment.color = NULL,
    segment.alpha = NULL,
    segment.size = 0, 
    size = 3
  ) +
  coord_polar(theta = "y") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Patient Cohort Demographics Distribution")
  ggsave("city.png")
}

args = commandArgs(trailingOnly = T)
city_file = args[1]

df <- read.csv(city_file, header = TRUE, sep = ",")
df <- get_processed_df(df)
plot_donut(df)
