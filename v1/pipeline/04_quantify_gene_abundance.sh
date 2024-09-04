#!/usr/bin/bash -l
#SBATCH -p short --mem 384gb -c 64 -N 1 -n 1 --out logs/quantify_gene_abundance_95.log

module load subread

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

# we need to construct a GTF file from the consensus mmseqs/cd-hit fasta CDS

OUTDIR=results/quantify
MAPPING=results/mapping
IDENTITY=95
GENOME=db/LsFMGC_AA_${IDENTITY}_rep.to_CDS.fasta
# this file was make in step 02 
GTF=$(echo -n $GENOME | perl -p -e 's/\.fasta/.gtf/')

mkdir -p $OUTDIR

if [ ! -s $OUTDIR/LsFMGC_coverage.tsv ]; then
    # not really paired end perhaps because the CDS are too short anyways
    featureCounts -p -o $OUTDIR/LsFMGC_coverage_${IDENTITY}.tsv -a $GTF -O -M \
        -G $GENOME -g gene_id -F GTF -t ORF --tmpDir $SCRATCH  --fraction \
        --countReadPairs $(find $MAPPING -size +0 -name "*.${IDENTITY}_rep.bam")
fi
