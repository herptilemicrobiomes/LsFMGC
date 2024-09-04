#!/usr/bin/bash -l
#SBATCH -p epyc --mem 384gb -c 96 -N 1 -n 1 --out logs/cluster_cds_pep.%A.log
MEM=350G
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
#module load cd-hit
#cd-hit

DOPEP=1
DOCDS=1

module load mmseqs2
module load samtools
module load bcftools

# or use mmseqs2
mkdir -p db
hostname
date
cat input_cds/*.fasta > db/LsFMGC_CDS_all.fasta
bgzip --threads $CPU db/LsFMGC_CDS_all.fasta
cat input_pep/*.fasta > db/LsFMGC_AA_all.fasta
bgzip --threads $CPU db/LsFMGC_AA_all.fasta
samtools faidx db/LsFMGC_CDS_all.fasta.gz
samtools faidx db/LsFMGC_AA_all.fasta.gz

date

if [ ! -z $DOPEP ]; then
    IN=input_pep
    DB=db/LsFMGC_AA
    if [ ! -f $DB ]; then
        mmseqs createdb db/LsFMGC_AA_all.fasta.gz $DB --compressed 1
    fi

    if [ ! -f $DB.idx ]; then
        mmseqs createindex $DB $SCRATCH
    fi
    # may need to mess with this clustering step using cd-hit parameters perhaps
    for identity in 1 0.95 0.9 0.5
    do
        NAME=$(perl -e "print $identity * 100")
        echo "Clustering AA at $identity"
        if [ ! -f ${DB}_${NAME}_cluster ]; then
            mmseqs cluster $DB ${DB}_${NAME}_cluster $SCRATCH --min-seq-id $identity -c 0.8 --cov-mode 1 --threads $CPU --split-memory-limit $MEM --kmer-per-seq 80
        fi
        mmseqs createsubdb ${DB}_${NAME}_cluster $DB ${DB}_${NAME}_cluster_rep
        mmseqs convert2fasta ${DB}_${NAME}_cluster_rep ${DB}_${NAME}_rep.fasta
        # perhaps cleanup aferwards   
    done
    date
fi

if [ ! -z $DOCDS ]; then
    IN=input_cds
    DB=db/LsFMGC_CDS

    if [ ! -f $DB ]; then
        mmseqs createdb db/LsFMGC_CDS_all.fasta.gz $DB --compressed 1
    fi
    if [ ! -f $DB.idx ]; then
        mmseqs createindex $DB $SCRATCH
    fi
    # may need to mess with this clustering step using cd-hit parameters perhaps
    for identity in 1 0.95 0.9 0.5
    do
        NAME=$(perl -e "print $identity * 100")
        echo "Clustering CDS at $identity"
#        if [ ! -f ${DB}_${NAME}_lincluster ]; then
#            mmseqs linclust $DB ${DB}_${NAME}_lincluster $SCRATCH --min-seq-id $identity -c 0.8 --cov-mode 1 --threads $CPU --split-memory-limit $MEM
#        fi
#        mmseqs createsubdb ${DB}_${NAME}_lincluster $DB ${DB}_${NAME}_lincluster_rep
#        mmseqs convert2fasta ${DB}_${NAME}_lincluster_rep ${DB}_${NAME}_rep.fasta   
    done
fi

