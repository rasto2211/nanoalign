library(ggplot2)
library(grid)
# library(xtable)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Args: [intersection.csv] [baseline.csv] [read_name]", call.=FALSE)
} 

data_samples <- read.csv(args[1], header=T)
colnames(data_samples) <- c("k","samples","inter")
data_baselines <- read.csv(args[2], header=F)
colnames(data_baselines) <- c("k", "inter") 

viterbi <- data_baselines[1:22,]
metrichor <- data_baselines[23:44,]
step_size <- 15

read_name <- args[3]

pdf(paste(read_name,"_separate.pdf", sep=""))
for (i in 9:30) {
  data_k <- data_samples[data_samples$k==i,]
  viterbi_k <- viterbi[viterbi$k==i,]
  metrichor_k <- metrichor[metrichor$k==i,]

  breaks_x <- seq(min(data_k$samples),max(data_k$samples),step_size)

  # Set ticks
  # min_y <- min(data_k$inter, viterbi_k$inter, metrichor_k$inter)*100
  # max_y <- max(data_k$inter, viterbi_k$inter, metrichor_k$inter)*100
  # lim_y_low <- round(min_y)
  # lim_y_upper <- round(max_y)
  # y_step <- (lim_y_upper - lim_y_low) / 22
  # breaks_y <- seq(max(lim_y_low-y_step,0), min(lim_y_upper+y_step,100), y_step)

  gg <- ggplot(data_k, aes(x=samples,y=inter*100)) + 
    geom_point(aes(col="samples"), col="blue") + 
    geom_line(aes(col="samples"), show.legend=TRUE, col="blue") +
    scale_x_continuous(breaks=breaks_x) +
    # scale_y_continuous(breaks=breaks_y) +
    geom_hline(data=viterbi_k, aes(yintercept=inter*100, fill="viterbi_k", group=k), 
		show.legend=TRUE, linetype="dashed", col="red") +
    geom_hline(data=metrichor_k,aes(yintercept=inter*100, fill="metrichor_k"),
		show.legend=TRUE, linetype="dashed", col="black") +
    xlab("Number of samples") +
    ylab("Intersection of kmers with reference (%)") +
    guides(colour = guide_legend(title="Length of kmer"), 
	    fill = guide_legend(title="Baselines", 
				  override.aes = list(colour=c("black", "red")), 
				  labels=c("metrichor_k", "viterbi_k"), order=2)) +
    ggtitle(paste(basename(read_name),"\n K=",i,sep=""))

  print(gg)
}

dev.off()
