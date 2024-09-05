#!/usr/bin/env Rscript

library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)

library(DESeq2)
# Define the command line options
#abundance <- commandArgs(trailingOnly = TRUE)[1]
#output <- commandArgs(trailingOnly = TRUE)[2]
abundance = 'results/quantify/LsFMGC_coverage_95.tsv'
metadata = 'wood_frog_2022_sample_metdata.tsv'
# sample metadata, fix treatment to not have spaces

sample_metadata = read_tsv(metadata,col_names=T) %>% 
  mutate_at(c("treatment"),~ gsub(' ','_',.x,fixed=TRUE))  %>% column_to_rownames('Biosample')

sample_metadata$treatment <- factor(sample_metadata$treatment)

# Read the abundance table using 'read_tsv' function
abundance_table <- read_tsv(abundance,skip=1,col_names=TRUE) %>% 
  select(-c(Chr,Start,End,Strand)) %>% filter(Length >= 100)

abundance.counts = rename_with(abundance_table,~ gsub("results/mapping/","", .x, fixed=TRUE)) %>% 
    rename_with(~ gsub("\\.\\d+_rep\\.bam","", .x, perl=TRUE)) %>% 
  select(-c(Length)) %>% column_to_rownames('Geneid')

all(rownames(sample_metadata) == colnames(abundance.counts))

abundance.counts <- abundance.counts[, rownames(sample_metadata)]

dds <- DESeqDataSetFromMatrix(countData = abundance.counts,
                              colData = sample_metadata,
                              design = ~ treatment)
dds

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

outfile = paste0(output,"/basidiobolus_comparison.pdf")
ggsave(outputfile,p,width=12,height=8)