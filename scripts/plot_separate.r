library(ggplot2)
library(grid)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Args: [intersection.csv] [baseline.csv] [read_name]", call.=FALSE)
} 

data_samples <- read.csv(args[1], header=T)
data_baselines <- read.csv(args[2], header=T)

# Every odd row contains Viterbi data and even contains Metrichor.
viterbi <- data_baselines[seq(1,nrow(data_baselines),2),]
metrichor <- data_baselines[seq(2,nrow(data_baselines),2),]
step_size <- 15

read_name <- args[3]

pdf(paste(read_name,"_separate.pdf", sep=""))

# Kmer range 9...30.
for (i in 9:30) {
  # Take data for the given kmers of size i.
  # Take only number of 
  data_k <- data_samples[data_samples$k==i,]
  viterbi_k <- viterbi[viterbi$k==i,]
  metrichor_k <- metrichor[metrichor$k==i,]

  breaks_x <- seq(min(data_k$num_samples),max(data_k$num_samples),step_size)

  gg <- ggplot(data_k,
	       aes(x=num_samples, 
		   y=(true_positive/(true_positive+false_negative))*100)) + 
    #geom_point(aes(col="num_samples"), col="blue") + 
    geom_line(aes(col="num_samples"), show.legend=TRUE, col="blue") +
    scale_x_continuous(breaks=breaks_x) +
    # scale_y_continuous(breaks=breaks_y) +
    geom_hline(data=viterbi_k, 
	       aes(yintercept=(true_positive/(true_positive+false_negative))*100, 
		   fill="viterbi_k", group=k), 
	       show.legend=TRUE, linetype="dashed", col="red") +
    geom_hline(data=metrichor_k,
	       aes(yintercept=(true_positive/(true_positive+false_negative))*100,
		   fill="metrichor_k"),
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
