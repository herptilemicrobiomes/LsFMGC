#!/usr/bin/bash -l
#SBATCH --mem 16gb -c 48 -N 1 -n 1 --out logs/map_reads_to_catalog.%a.log

module load workspace/scratch
module load bwa-mem2
module load samtools

INPUT=input_cds
IDENTITY=95
DB=db/LsFMGC_${IDENTITY}_cluster_rep.fasta # or whatever name you chose, and indicating perhaps what clustering method you used
SAMPFILE=wood_frog_samples.csv
WORK=working
OUT=results/mapping
TEMP=$SCRATCH
mkdir -p $OUT
# index the db for bwa mapping
if [[ ! -f $DB.pac ]]; then
    bwa-mem2 index $DB
fi

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
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
    SRTED=$OUT/$STRAIN.cram
    READGROUP="@RG\tID:$STRAIN\tSM:$STRAIN\tLB:$STRAIN\tPL:illumina\tCN:UHM"
    bwa-mem2 mem -t $CPU -o $DB $LEFT $RIGHT | samtools sort --threads 8 -O bam -o $SRTED -T $TEMP -
done
