

Pre-processing Sample list
===
Problem: we have the large 90+ UHM fecal datasets, but only want woodfrog experiments.

How to take list of sample names to keep and join to the samples.csv data
Datasets need to be sorted:
I have made these file copies in the `lib` folder
```
sort wood_frog_samples.txt > s
mv s wood_frog_samples.txt
```

Then run join, it will default to using first column in each dataset as join, columns are delimited by ','

```
cd lib
join -t,  wood_frog_samples.txt  samples.csv > ../wood_frog_samples.csv
```

Adding metadata
===
We need to add metadata to the wood_frog_samples.csv file eventually (or another metadata.csv file if you prefer) to provide the data we will use to stratify plots later on.

Algorithm
===
1. Assemble metagenomes (metashot pipeline)
2. Predict genes in these metagenomes see [10_predict_by_assembl.sh](https://github.com/herptilemicrobiomes/MAG_Fecal/blob/main/pipeline/10_predict_by_assembl.sh). Note this set needs to have the genes renamed so that the UHM biosample prefixes these. -- [`pipeline/00_setup.sh`](pipeline/00_setup.sh)
3. filter / trim reads (FastP - this was done in step 1 but we threw away these files so re-run) -- [`pipeline/01_filter_trim.sh`](pipeline/01_filter_trim.sh)
4. Build non-redundant gene catalog across all samples with CD-HIT. This is LsFMGC set -- `pipeline/02_build_gene_catalog.sh`
5. Map clean/filter reads from step 2 to the non-redundant CDS LsFMGC -- `pipeline/03_map_reads_to_catalog.sh`
6. Quantify the per-gene depth with featureCounts -- `pipeline/04_quantify_gene_abundance.sh` to produce a catalog with depths matrix (eg rows are the non-redudnant genes from step 5/CD-HIT and columns are each UHM Biosample, the values are the RPM (reads-per-million) normalized read-depth)).
7. Plot these in R -- need a script
