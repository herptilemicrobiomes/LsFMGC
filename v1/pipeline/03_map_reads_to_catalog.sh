#!/usr/bin/bash -l
#SBATCH --mem 96gb -c 48 -N 1 -n 1 --out logs/map_reads_to_catalog_95.%a.log

module load workspace/scratch
module load bwa-mem2
module load samtools

INPUT=input_cds
IDENTITY=95
DB=db/LsFMGC_AA_${IDENTITY}_rep.to_CDS.fasta
SAMPFILE=wood_frog_samples.csv
WORK=working
OUT=results/mapping
TEMP=$SCRATCH
mkdir -p $OUT

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
    SORTTHREAD=8
    CPU=$(expr $CPU - $SORTTHREAD)
    if $CPU -lt 1; then
        CPU=$SLURM_CPUS_ON_NODE
        SORTTHREAD=1
    fi
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
fi
if [ -z $N ]; then
    echo "cannot run without a number provided either cmdline or --array in sbatch"
    exit
fi
IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN SHOTGUN 
do  
    LEFT=$WORK/$STRAIN/${STRAIN}_R1.fq.gz
    RIGHT=$WORK/$STRAIN/${STRAIN}_R2.fq.gz
    SRTED=$OUT/$STRAIN.${IDENTITY}_rep.bam
    READGROUP="@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$STRAIN\tPL:illumina\tCN:UHM"
    bwa-mem2 mem -Y -t $CPU $DB $LEFT $RIGHT | samtools sort --threads $SORTTHREAD -O bam -o $SRTED -T $TEMP -
    samtools index $SRTED
done
