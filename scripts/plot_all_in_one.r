require(RColorBrewer)
library(ggplot2)
# library(xtable)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Args: [intersection.csv] [baseline.csv] [read_name]", call.=FALSE)
} 

data_samples <- read.csv(args[1], header=T)
data_baselines <- read.csv(args[2], header=T)

n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
			    rownames(qual_col_pals)))

# Every odd row contains Viterbi data and even contains Metrichor.
viterbi <- data_baselines[seq(1,nrow(data_baselines),2),]
metrichor <- data_baselines[seq(2,nrow(data_baselines),2),]
step_size <- 15

read_name <- args[3]

identities <-(data_samples$true_positive/
	      (data_samples$true_positive+data_samples$false_negative))*100
identities <- cbind(identities, (data_baselines$true_positive/
	      (data_baselines$true_positive+data_baselines$false_negative))*100)
x_breaks <- seq(min(data_samples$num_samples), 
		max(data_samples$num_samples),step_size)
y_breaks <- seq(round(min(identities)), round(max(identities)),5)

gg <- ggplot(data_samples, 
	     aes(x=num_samples,
		 y=(true_positive/(true_positive+false_negative))*100, 
		 group=factor(k))) + 
  #geom_point(aes(col=factor(k))) + 
  geom_line(aes(group=factor(k),col=factor(k)), show.legend=TRUE) +
  scale_x_continuous(breaks=x_breaks) +
  scale_y_continuous(breaks=y_breaks) +
  geom_hline(data=viterbi, 
	     aes(yintercept=(true_positive/(true_positive+false_negative))*100, 
		 fill="Viterbi", group=k), 
	     show.legend=TRUE, linetype="dashed", col="red") +
  geom_hline(data=metrichor,
	     aes(yintercept=(true_positive/(true_positive+false_negative))*100, 
		 fill="Metrichor"),
	      show.legend=TRUE, linetype="dashed", col="black") +
  xlab("Number of samples") +
  ylab("Sensitivity (%)") +
  ggtitle(basename(read_name)) +
  scale_colour_manual(values = col_vector) +
  guides(colour = guide_legend(title="Length of kmer"), 
	 fill = guide_legend(title="Baselines", 
			     override.aes = list(colour=c("black", "red")), 
			     labels=c("Metrichor", "Viterbi"), order=2)) +
  theme(plot.title = element_text(size=12))

ggsave(paste(read_name,"_all_in_one.pdf", sep=""), gg)
