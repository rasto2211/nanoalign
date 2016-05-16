args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Args: [data.csv] [output_box_plot.pdf]", call.=FALSE)
} 

identities <- read.csv(args[1], header=T)
  
samples <- identities[1:(nrow(identities)-2),]
viterbi <- identities[nrow(identities)-1,]
metrichor <- identities[nrow(identities),]

pdf(args[2])

cols <- ncol(identities)
min <- min(identities, na.rm=T)
max <- max(identities, na.rm=T)

boxplot(samples, notch=T, names=1:cols, ylim=c(min,max))
points(y=viterbi, x=1:cols, col="red", pch=4)
points(y=metrichor, x=1:cols , col="blue", pch=4)

legend("bottomright", c("Viterbi", "Metrichor"), 
    col=c("red", "blue"), pch=c(4,4))

barplot(colSums(is.na(samples)), names.arg=1:cols)
 
dev.off()
