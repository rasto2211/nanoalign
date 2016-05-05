library(ggplot2)
library(grid)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Args: [compound_samples_intersection.csv] [compound_baselines_intersection.csv] ", 
       call.=FALSE)
} 

data_samples <- read.csv(args[1], header=T)
data_baselines <- read.csv(args[2], header=T)

# Every odd row contains Viterbi data and even contains Metrichor.
viterbi <- data_baselines[seq(1,nrow(data_baselines),2),]
metrichor <- data_baselines[seq(2,nrow(data_baselines),2),]
step_size <- 15

pdf("compound_intersection.pdf")

# Kmer range 9...30.
for (i in 9:30) {
  # Take data for the given kmers of size i.
  data_k <- data_samples[data_samples$k==i,]
  viterbi_k <- viterbi[viterbi$k==i,]
  metrichor_k <- metrichor[metrichor$k==i,]

  agg_by_num_samples <- function(df) {
    aggregate(cbind(true_positive,true_negative, false_positive,
		  false_negative) ~ num_samples, data=df, FUN=sum)
  }
  data_agg <- agg_by_num_samples(data_k)

  # Sum over columns and keep names of columns.
  sum_cols_df <- function(df) {
    res <- data.frame(c(sum(df$true_positive)),
		      c(sum(df$true_negative)),
		      c(sum(df$false_positive)),
		      c(sum(df$false_negative)))
    colnames(res) <- c("true_positive", "true_negative", "false_positive", 
		    "false_negative")
    res
  }
  viterbi_agg <- sum_cols_df(viterbi_k)
  metrichor_agg <- sum_cols_df(metrichor_k)

  breaks_x <- seq(min(data_agg$num_samples),max(data_agg$num_samples),step_size)

  # Sensitivity(recall - hit rate) plot. - TP/(TP+FN)
  gg <- ggplot(data_agg,
	       aes(x=num_samples, 
		   y=(true_positive/(true_positive+false_negative))*100)) + 
    geom_line(aes(col="num_samples"), show.legend=TRUE, col="blue") +
    scale_x_continuous(breaks=breaks_x) +
    geom_hline(data=viterbi_agg, 
	       aes(yintercept=(true_positive/(true_positive+false_negative))*100, 
		   fill="viterbi_agg"), 
	       show.legend=TRUE, linetype="dashed", col="red") +
    geom_hline(data=metrichor_agg,
	       aes(yintercept=(true_positive/(true_positive+false_negative))*100,
		   fill="metrichor_agg"),
	       show.legend=TRUE, linetype="dashed", col="black") +
    xlab("Number of samples") +
    ylab("Compound sensitivity - TP/(TP+FN)") +
    guides(colour = guide_legend(title="Length of kmer"), 
	   fill = guide_legend(title="Baselines", 
			       override.aes = list(colour=c("black", "red")), 
			       labels=c("metrichor_agg", "viterbi_agg"), 
			       order=2)) + ggtitle(paste("K=",i,sep=""))
  print(gg)

  # Second plot -- (TP+FP)/(TP+FN).
  gg2 <- ggplot(data_agg, aes(x=num_samples,
	    y=(true_positive+false_positive)/(true_positive+false_negative))) +
    geom_line(aes(col="num_samples"), show.legend=TRUE, col="blue") +
    scale_x_continuous(breaks=breaks_x) +
    geom_hline(data=viterbi_agg, aes(
      yintercept=(true_positive+false_positive)/(true_positive+false_negative), 
      fill="viterbi_agg"), show.legend=TRUE, linetype="dashed", col="red") +
    geom_hline(data=metrichor_agg,
	       aes(yintercept=(true_positive+false_positive)/
		      (true_positive+false_negative), fill="metrichor_agg"),
	       show.legend=TRUE, linetype="dashed", col="black") +
    xlab("Number of samples") +
    ylab("found kmers/ref. kmers") +
    guides(colour = guide_legend(title="Length of kmer"), 
	   fill = guide_legend(title="Baselines", 
			       override.aes = list(colour=c("black", "red")), 
			       labels=c("metrichor_agg", "viterbi_agg"), order=2)) + ggtitle(paste("K=",i,sep=""))

  print(gg2)

  # Precision TP/(TP+FP)
  gg3 <- ggplot(data_agg, aes(x=num_samples, 
			  y=(true_positive/(true_positive+false_negative))*100)) +
    geom_line(aes(col="num_samples"), show.legend=TRUE, col="blue") +
    scale_x_continuous(breaks=breaks_x) +
    geom_hline(data=viterbi_agg, aes(
      yintercept=(true_positive/(true_positive+false_negative))*100, 
      fill="viterbi_agg"), show.legend=TRUE, linetype="dashed", col="red") +
    geom_hline(data=metrichor_agg,
	    aes(yintercept=(true_positive/(true_positive+false_negative))*100, 
		   fill="metrichor_agg"),
		   show.legend=TRUE, linetype="dashed", col="black") +
    xlab("Number of samples") +
    ylab("Compound precision - TP/(TP+FP)") +
    guides(colour = guide_legend(title="Length of kmer"), 
	   fill = guide_legend(title="Baselines", 
			       override.aes = list(colour=c("black", "red")), 
			       labels=c("metrichor_agg", "viterbi_agg"), order=2)) + ggtitle(paste("K=",i,sep=""))

  print(gg3)
}

dev.off()
