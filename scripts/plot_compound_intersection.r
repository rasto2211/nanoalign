# Plots multiple plots for every $k \in \{9 \dots 30\}$ and various compound 
# statistics.
# Args: [compound_samples_intersection.csv] [compound_baselines_intersection.csv]  
# - `compound_samples_intersection.csv` contains rows produced by 
# `create_compound_csv.py` from SAM files of samples. 
# - `compound_baselines_intersection.csv` contains also rows produced 
# by the script mentioned above
# but from Viterbi and Metrichor sequences for various values of $k$ mentioned
# above.

library(ggplot2)
library(grid)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Args: [compound_samples_intersection.csv] [compound_baselines_intersection.csv] ", 
       call.=FALSE)
} 

# Source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##########################################################################
# end of multiplot

data_samples <- read.csv(args[1], header=T)
data_baselines <- read.csv(args[2], header=T)

# Every odd row contains Viterbi data and even contains Metrichor.
viterbi <- data_baselines[seq(1,nrow(data_baselines),2),]
metrichor <- data_baselines[seq(2,nrow(data_baselines),2),]
step_size <- 50

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
    ylab("Compound sensitivity (%)") +
    guides(colour = guide_legend(title="Length of kmer"), 
	   fill = guide_legend(title="Baselines", 
			       override.aes = list(colour=c("black", "red")), 
			       labels=c("metrichor_agg", "viterbi_agg"), 
			       order=2)) + ggtitle(paste("K=",i,sep=""))

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
    ylab("found kmers / reference kmers") +
    guides(colour = guide_legend(title="Length of kmer"), 
	   fill = guide_legend(title="Baselines", 
			       override.aes = list(colour=c("black", "red")), 
			       labels=c("metrichor_agg", "viterbi_agg"), order=2)) + ggtitle(paste("K=",i,sep=""))


  # Precision TP/(TP+FP)
  gg3 <- ggplot(data_agg, aes(x=num_samples, 
			  y=(true_positive/(true_positive+false_positive))*100)) +
    geom_line(aes(col="num_samples"), show.legend=TRUE, col="blue") +
    scale_x_continuous(breaks=breaks_x) +
    geom_hline(data=viterbi_agg, aes(
      yintercept=(true_positive/(true_positive+false_positive))*100, 
      fill="viterbi_agg"), show.legend=TRUE, linetype="dashed", col="red") +
    geom_hline(data=metrichor_agg,
	    aes(yintercept=(true_positive/(true_positive+false_positive))*100, 
		   fill="metrichor_agg"),
		   show.legend=TRUE, linetype="dashed", col="black") +
    xlab("Number of samples") +
    ylab("Compound precision (%)") +
    guides(colour = guide_legend(title="Length of kmer"), 
	   fill = guide_legend(title="Baselines", 
			       override.aes = list(colour=c("black", "red")), 
			       labels=c("metrichor_agg", "viterbi_agg"), order=2)) + ggtitle(paste("K=",i,sep=""))

  # Specificity TN/(FP+TN)
  gg4 <- ggplot(data_agg, aes(x=num_samples, 
			y=100*(true_negative/(false_positive+true_negative)))) +
    geom_line(aes(col="num_samples"), show.legend=TRUE, col="blue") +
    scale_x_continuous(breaks=breaks_x) +
    geom_hline(data=viterbi_agg, aes(
      yintercept=100*(true_negative/(false_positive+true_negative)), 
      fill="viterbi_agg"), show.legend=TRUE, linetype="dashed", col="red") +
    geom_hline(data=metrichor_agg,
	    aes(yintercept=100*(true_negative/(false_positive+true_negative)), 
		   fill="metrichor_agg"),
		   show.legend=TRUE, linetype="dashed", col="black") +
    xlab("Number of samples") +
    ylab("Compound specificity (%)") +
    guides(colour = guide_legend(title="Length of kmer"), 
	   fill = guide_legend(title="Baselines", 
			       override.aes = list(colour=c("black", "red")), 
			       labels=c("metrichor_agg", "viterbi_agg"), 
			       order=2)) + ggtitle(paste("K=",i,sep=""))

    multiplot(gg, gg2, gg3, gg4, cols=2)
}

dev.off()
