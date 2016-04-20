#!/usr/bin/Rscript

x <- read.csv('stdin', header = T); 
no_na <- na.omit(x)[,1]
cat(sprintf("%f\t%f\t%f\t%d\t%f\t%f\n", 
      sd(no_na), 
      mean(no_na), 
      median(no_na), 
      sum(is.na(x)), 
      quantile(no_na, 0.25), 
      quantile(no_na, 0.75)));
