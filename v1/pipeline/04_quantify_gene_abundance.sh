#!/usr/bin/bash -l
#SBATCH -p short --mem 64gb -c 16 -N 1 -n 1 --out logs/quantify_gene_abundance.log

module load subread

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

# we need to construct a GTF file from the consensus mmseqs/cd-hit fasta CDS
GTF=$(echo -n $DB | perl -p -e 's/\.fasta/.gtf/')
OUTDIR=results/quantify
MAPPING=results/mapping
mkdir -p $OUTDIR
GENOME=db/LsFMGC_cluster_rep.fasta
if [[ ! -f $DB.gtf ]]; then
    python scripts/fasta_to_gtf.py $DB > $GTF
fi
if [ ! -s $OUTDIR/LsFMGC_coverage.tsv ]; then
    featureCounts -p  -o $OUTDIR/LsFMGC_coverage.tsv -G $GENOME -J -g gene_id -F GTF --countReadPairs $(find $MAPPING -size +0 -name "*.bam")
fi