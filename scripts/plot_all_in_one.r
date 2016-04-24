require(RColorBrewer)
library(ggplot2)
# library(xtable)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Args: [intersection.csv] [baseline.csv] [read_name]", call.=FALSE)
} 

data_samples <- read.csv(args[1], header=F)
colnames(data_samples) <- c("k","samples","inter")
data_baselines <- read.csv(args[2], header=F)
colnames(data_baselines) <- c("k", "inter") 

n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
			    rownames(qual_col_pals)))

viterbi <- data_baselines[1:22,]
metrichor <- data_baselines[23:44,]
step_size <- 15

read_name <- args[3]

gg <- ggplot(data_samples, aes(x=samples,y=inter*100, group=factor(k))) + 
  geom_point(aes(col=factor(k))) + 
  geom_line(aes(group=factor(k),col=factor(k)), show.legend=TRUE) +
  scale_x_continuous(breaks=seq(min(data_samples$samples),max(data_samples$samples),step_size)) +
  scale_y_continuous(breaks=seq(round(min(data_samples$inter*100)),
		      round(max(data_samples$inter*100)),5)) +
  geom_hline(data=viterbi, aes(yintercept=inter*100, fill="Viterbi", group=k), 
	      show.legend=TRUE, linetype="dashed", col="red") +
  geom_hline(data=metrichor,aes(yintercept=inter*100, fill="Metrichor"),
	      show.legend=TRUE, linetype="dashed", col="black") +
  xlab("Number of samples") +
  ylab("Intersection of kmers with reference (%)") +
  ggtitle(basename(read_name)) +
  scale_colour_manual(values = col_vector) +
  guides(colour = guide_legend(title="Length of kmer"), 
	  fill = guide_legend(title="Baselines", 
				override.aes = list(colour=c("black", "red")), 
				labels=c("Metrichor", "Viterbi"), order=2)) +
 theme(plot.title = element_text(size=12))

ggsave(paste(read_name,"_all_in_one.pdf", sep=""), gg)
