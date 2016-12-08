#!/usr/bin/env Rscript

#Rscript -e 'install.packages("ggplot2", dep = TRUE, repos="http://cran.rstudio.com/")'
library(ggplot2)
library(tools)

args<-commandArgs(TRUE)

fn=args[1]
name <- file_path_sans_ext(args[1])
n <- basename(name)
cov <- read.delim(fn, header=FALSE)
names(cov) <- c("CHR", "POSITION", "COVERAGE")

options(bitmapType='cairo')
png(filename=args[2])
ggplot(cov, aes(POSITION, COVERAGE)) + geom_point() + geom_hline(yintercept = 30, color = "green") +
ggtitle(paste("Depth of Coverage for",n))
junk <- dev.off()
