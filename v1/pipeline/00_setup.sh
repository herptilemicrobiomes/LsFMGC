#!/usr/bin/bash -l
#SBATCH -p short -c 16 --mem 8gb --out logs/00_setup_rename.log

#INDIR=/bigdata/stajichlab/shared/projects/Herptile/Metagenome/Fecal/Proteins/by_assembly
INDIR=by_assembly # this folder contains files named XX.cds.fa.gz XX.gff.gz which are prodigal output from a per-metagenome run
CDSDIR=input_cds
PEPDIR=input_pep
mkdir -p $CDSDIR $PEPDIR
IFS=,
SAMPLES=wood_frog_samples.csv

tail -n +2 $SAMPLES | while read BIOSAMPLE READFILE
do
    if [[ ! -f $INDIR/${BIOSAMPLE}.cds.fa.gz ]]; then
        echo "skipping $BIOSAMPLE cannot find $INDIR/${BIOSAMPLE}.cds.fa.gz"
        continue
    fi
    if [[ ! -f $INDIR/${BIOSAMPLE}.aa.fa.gz ]]; then
        echo "skipping $BIOSAMPLE cannot find $INDIR/${BIOSAMPLE}.aa.fa.gz"
        continue
    fi
    export BIOSAMPLE
    pigz -dc $INDIR/${BIOSAMPLE}.cds.fa.gz | perl -p -e 's/^>(\S+).+/>$ENV{BIOSAMPLE}__$1/' > $CDSDIR/${BIOSAMPLE}.fasta
    pigz -dc $INDIR/${BIOSAMPLE}.aa.fa.gz | perl -p -e 's/^>(\S+).+/>$ENV{BIOSAMPLE}__$1/' > $PEPDIR/${BIOSAMPLE}.fasta

done
