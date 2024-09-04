
This works by taking prodigal results from annotation of metagenome contigs (not MAGs) stored in `by_assembly` in the form of BIOSAMPLE.cds.fa.gz and BIOSAMPLE.aa.fa.gz. These sequences in each file are renamed so the BIOSAMPLE is prefixed changing the name from `>k141_103518_1` to `>UHM20.35743__k141_103518_1` for example. This is done by the script [00_setup.sh](pipeline/00_setup.sh) and uses the [wood_frog_samples.csv](wood_frog_samples.csv).

The folder `input` has the raw reads the naming of BioSample to fastq files is based on the `sample.csv` or in this case [wood_frog_samples.csv](wood_frog_samples.csv) file.  These raw reads are trimmed and QCed with fastp and these are stored in a folder `work/BIOSAMPLE/BIOSAMPLE_R[12].fastq.gz`. This step is done by the script [01_filter_trim.sh](pipeline/01_filter_trim.sh).

To build the non-redundant gene catalog from the ~20-30 M predicted proteins/genes we use [mmseqs2](https://github.com/soedinglab/MMseqs2) for clustering the proteins. One can also use cd-hit but I think memory requirements are smaller with mmseqs2, it would be worth comparing and profiling results.

First the proteins are collapsed into clusters based on 100%, 95%, 90%, and 50% identity and a representative chosen from each cluster. We can use these different clustering cutoffs to see how they impact our interpretation of the gene catalog.
I am trying to cluster the CDS directly but it seems to take too long. But it will be helpful to contrast that. This clustering with mmseqs2 is run with code [02_build_gene_catalog.sh](pipeline/02_build_gene_catalog.sh).

Once the clusters are created we can get the corresponding CDS sequence for the named representative sequence back out to make a DNA database of the representative genes. The abundance of each gene in each metagenome sample is computed by using [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) (wich should be faster) or just [BWA](https://github.com/lh3/bwa). We can also try bowtie. This produces BAM files which can then be processed with [featureCounts](https://subread.sourceforge.net/featureCounts.html) part of (subread)[https://subread.sourceforge.net/].


Some things I had to do to make this work
===
**Metadata**
(note we need to build the metadata complete CSV that has the info about each sample eg Basidiobolus status, date, Basidiobolus strain,etc).

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
