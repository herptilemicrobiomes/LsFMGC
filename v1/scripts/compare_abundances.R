#!/usr/bin/env Rscript

library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)

library(DESeq2)
# Define the command line options
#abundance <- commandArgs(trailingOnly = TRUE)[1]
#output <- commandArgs(trailingOnly = TRUE)[2]
abundance = 'results/quantify/LsFMGC_coverage_95.tsv.gz'

metadata = 'wood_frog_2022_sample_metdata.tsv'
outfile = "basidiobolus_MA_comparison.pdf"
pdf(outfile)
# sample metadata, fix treatment to not have spaces

sample_metadata = read_tsv(metadata,col_names=T) %>% 
  mutate_at(c("treatment"),~ gsub(' ','_',.x,fixed=TRUE)) %>%  
  column_to_rownames('Biosample') %>% arrange(treatment)

sample_metadata$treatment <- factor(sample_metadata$treatment)

# Read the abundance table using 'read_tsv' function
abundance.counts <- read_tsv(abundance,skip=1,col_names=TRUE) %>% 
  select(-c(Chr,Start,End,Strand)) %>% filter(Length >= 100)

abundance.counts = rename_with(abundance.counts,~ gsub("results/mapping/","", .x, fixed=TRUE)) %>% 
    rename_with(~ gsub("\\.\\d+_rep\\.bam","", .x, perl=TRUE)) %>% 
  select(-c(Length)) %>% column_to_rownames('Geneid')

# check that sample names are overlapping and same but don't have to be same order
all(rownames(sample_metadata) %in% colnames(abundance.counts))
all(colnames(abundance.counts) %in% rownames(sample_metadata))
# reorder columns based on the sample_metadata
abundance.counts <- abundance.counts[, rownames(sample_metadata)]
all(rownames(sample_metadata) == colnames(abundance.counts))

# force to integer as if we run with --fraction option it will generate fractional counts
abundance.counts <- abundance.counts %>% mutate(across(everything(), as.integer))

dds <- DESeqDataSetFromMatrix(countData = abundance.counts,
                              colData = sample_metadata,
                              design = ~ treatment)
dds

smallestGroupSize <- 20
keep <- rowSums(counts(dds) >= 100) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
res
resNo_vs_Basid1717 <- results(dds, contrast=c("treatment","control","STP1717.1_pilot"))
resNo_vs_Basid1710 <- results(dds, contrast=c("treatment","control","STP1710.7_pilot"))
resBasid1717_vs_Basid1710 <- results(dds, contrast=c("treatment","STP1717.1_pilot","STP1710.7_pilot"))

overrep_No_vs_Basid1710 <- subset(resNo_vs_Basid1710,resNo_vs_Basid1710$log2FoldChange > 2)
underrep_No_vs_Basid1710 <- subset(resNo_vs_Basid1710,resNo_vs_Basid1710$log2FoldChange < -2)

overrep_No_vs_Basid1717 <- subset(resNo_vs_Basid1717,resNo_vs_Basid1717$log2FoldChange > 2)
underrep_No_vs_Basid1717 <- subset(resNo_vs_Basid1717,resNo_vs_Basid1717$log2FoldChange < -2)

overrep_Basid1717_vs_Basid1710 <- subset(resBasid1717_vs_Basid1710,resBasid1717_vs_Basid1710$log2FoldChange > 2)
underrep_Basid1717_vs_Basid1710 <- subset(resBasid1717_vs_Basid1710,resBasid1717_vs_Basid1710$log2FoldChange < -2)

pdf("results/quantify/LsFMGC_DESeq_plots.pdf")
plotMA(resNo_vs_Basid1717, main="Control vs STP1717.1")
plotMA(resNo_vs_Basid1710, main="Control vs STP1710.7")
plotMA(resBasid1717_vs_Basid1710, main="STP1717.1 vs STP1710.7")

resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

resAsh <- lfcShrink(dds, coef=2, type="ashr")
plotMA(resAsh, ylim=c(-3,3), main="ashr")

dev.off()

resfile = 'results/quantify/LsFMGC_DESeq_Ctl_vs_Basid1717_95.csv'
write.csv(as.data.frame(resNo_vs_Basid1717),file=resfile)

resfile = 'results/quantify/LsFMGC_DESeq_Ctl_vs_Basid1717_95_over.csv'
write.csv(as.data.frame(overrep_No_vs_Basid1717),file=resfile)
resfile = 'results/quantify/LsFMGC_DESeq_Ctl_vs_Basid1717_95_under.csv'
write.csv(as.data.frame(underrep_No_vs_Basid1717),file=resfile)

##

resfile = 'results/quantify/LsFMGC_DESeq_Ctl_vs_Basid1710_95.csv'
write.csv(as.data.frame(resNo_vs_Basid1710),file=resfile)

resfile = 'results/quantify/LsFMGC_DESeq_Ctl_vs_Basid1710_95_over.csv'
write.csv(as.data.frame(overrep_No_vs_Basid1710),file=resfile)
resfile = 'results/quantify/LsFMGC_DESeq_Ctl_vs_Basid1710_95_under.csv'
write.csv(as.data.frame(underrep_No_vs_Basid1710),file=resfile)

## 
resfile = 'results/quantify/LsFMGC_DESeq_Basid1717_vs_Basid1710_95.csv'
write.csv(as.data.frame(resBasid1717_vs_Basid1710),file=resfile)

resfile = 'results/quantify/LsFMGC_DESeq_Basid1717_vs_Basid1710_95_over.csv'
write.csv(as.data.frame(overrep_Basid1717_vs_Basid1710),file=resfile)
resfile = 'results/quantify/LsFMGC_DESeq_Basid1717_vs_Basid1710_95_under.csv'
write.csv(as.data.frame(underrep_Basid1717_vs_Basid1710),file=resfile)

