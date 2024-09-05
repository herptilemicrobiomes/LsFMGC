#!/usr/bin/env Rscript

library(tidyverse)
library(readr)
library(dplyr)

frog_samples = read_csv("wood_frog_samples.csv",col_names=T) %>% select(c(Biosample)) %>%
  separate(Biosample, c('sample_name','biosample'))
all_samples = read_csv("lib/wood_frog_2022_all_samples.csv",col_names=T) %>% 
  select(c(sample_name, treatment)) %>% distinct()

combined_md = frog_samples %>% inner_join(all_samples,by="sample_name") %>% 
  unite("Biosample",sample_name:biosample, sep=".")
write_tsv(combined_md,"wood_frog_2022_sample_metdata.tsv")
