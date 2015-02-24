#!/local/gensoft2/exe/R/3.1.2/scripts/Rscript

#plots a coverage graphe from a .mpileup coverage file

library(ggplot2)
library(tools)

args <- commandArgs(trailingOnly = TRUE)

d=read.table(pipe(paste("cut -f 1,2,4", args[1])),sep="\t")
names(d)=c("chr","pos","cov")

out=paste(file_path_sans_ext(args[1]),"_coverage_graph.jpeg",sep="")

jpeg(out)
ggplot(d, aes(x=pos,y=cov)) + geom_histogram(stat="identity",binwidth=100, colour="black", fill="white")+ facet_wrap(~chr, scales = "free")
dev.off()
