df <- read.csv("stdin")

df_filtered <- df[df$id>=65 & df$seq_len <= 25000 & df$clip <= 1.3 & df$seq_len>=300,]

write.csv(df_filtered, "filtered_by_stats.csv", row.names=F)
