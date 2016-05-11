df <- read.csv("stdin")

attach(df)
match <- M - (edit_dist - I - D)
df_filtered <- df[match/(edit_dist+S_CLIP+H_CLIP+match) > 0.63 & read_len <= 30000,]

write.csv(df_filtered, "filtered_by_stats.csv", row.names=F)
