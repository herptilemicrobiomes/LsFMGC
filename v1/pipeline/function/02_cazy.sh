#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 64 --mem 64gb --out logs/cazy.log

CPU=2
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

#N=${SLURM_ARRAY_TASK_ID}
#if [ -z $N ]; then
#    N=$1
#    if [ -z $N ]; then
#        echo "need to provide a number by --array or cmdline"
#        exit
#    fi
#fi

module load dbcanlight
module load workspace/scratch

INFILE=db/LsFMGC_AA_95_rep.fasta
TEMP=db/$(basename $INFILE .fasta)__split
PREFIX=LsFMGC
OUTDIR=results/function/cazy
mkdir -p $OUTDIR

IN=$TEMP/${PREFIX}.$N
time dbcanlight search -i $INFILE -m cazyme -o $OUTDIR -t $CPU -b 100000
time dbcanlight search -i $INFILE -m sub -o $OUTDIR -t $CPU -b 100000
time dbcanlight search -i $INFILE -m diamond -o $OUTDIR -t $CPU -b 100000