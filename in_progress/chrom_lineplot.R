#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3) stop("Usage: Rscript chrom_lineplot.R file1.txt file2.txt N")

library(ggplot2)

f1 <- args[1]; f2 <- args[2]; n <- as.integer(args[3])

p1 <- sub("\\.txt$","",basename(f1))
p2 <- sub("\\.txt$","",basename(f2))
out <- paste0(p1,"_vs_",p2)

d1 <- read.table(f1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
d2 <- read.table(f2, header=TRUE, sep="\t", stringsAsFactors=FALSE)

d1 <- d1[!grepl("_region$", d1$chrom), ]
d2 <- d2[!grepl("_region$", d2$chrom), ]

chrs <- unique(d1$chrom)[1:min(n, length(unique(d1$chrom)))]

d1 <- d1[d1$chrom %in% chrs, ]; d2 <- d2[d2$chrom %in% chrs, ]
d1$set <- p1; d2$set <- p2

d <- rbind(d1, d2)
d$chrom <- factor(d$chrom, levels=chrs)
d$x <- as.numeric(d$chrom)

p <- ggplot(d, aes(x=x, y=mean, color=set, group=set)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks=seq_along(chrs), labels=chrs) +
  labs(x="Chromosome", y="Mean") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1))

ggsave(paste0(out,".png"), p, width=12, height=5, dpi=150)
ggsave(paste0(out,".pdf"), p, width=12, height=5)
