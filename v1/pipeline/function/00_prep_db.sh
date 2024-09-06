#!/usr/bin/bash -l
#SBATCH -p short -c 2 -N 1  --mem 2gb --out=logs/prep_db_95.log 

module load bioperl
INFILE=db/LsFMGC_AA_95_rep.fasta
TEMP=db/$(basename $INFILE .fasta)__split
mkdir -p $TEMP

bp_dbsplit --size 10000 --prefix $TEMP/LsFMGC -if fasta -of fasta $INFILE

ls $TEMP | wc -l
