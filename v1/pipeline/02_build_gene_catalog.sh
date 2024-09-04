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
module load bwa-mem2

# or use mmseqs2
mkdir -p db
hostname
date

CDSDB=db/LsFMGC_CDS_all.fasta
AADB=db/LsFMGC_AA_all.fasta
if [ ! -f $CDSDB.gz ]; then
    cat input_cds/*.fasta > $CDSDB
    bgzip --threads $CPU $CDSDB
    samtools faidx $CDSDB.gz
fi
if [ ! -f $AADB.gz ]; then
    cat input_pep/*.fasta > $AADB
    bgzip --threads $CPU $AADB
    samtools faidx $AADB.gz
fi


date

if [ ! -z $DOPEP ]; then
    IN=input_pep
    DB=db/LsFMGC_AA
    if [ ! -f $DB ]; then
        mmseqs createdb $AADB.gz $DB --compressed 1
    fi

    if [ ! -f $DB.idx ]; then
        mmseqs createindex $DB $SCRATCH
    fi
    # may need to mess with this clustering step using cd-hit parameters perhaps
    for identity in 1 0.95 0.9 0.5
    do
        NAME=$(perl -e "print $identity * 100")
        AA2CDS=${DB}_${NAME}_rep.to_CDS.fasta
        echo "Clustering AA at $identity"
        # won't overwrite if it already exists
        mmseqs cluster $DB ${DB}_${NAME}_cluster $SCRATCH --min-seq-id $identity -c 0.8 --cov-mode 1 --threads $CPU --split-memory-limit $MEM --kmer-per-seq 80

        if [ ! -f ${DB}_${NAME}_cluster_rep ]; then
            mmseqs createsubdb ${DB}_${NAME}_cluster $DB ${DB}_${NAME}_cluster_rep
        fi
        if [ ! -f ${DB}_${NAME}_cluster_rep.fasta ]; then
            mmseqs convert2fasta ${DB}_${NAME}_cluster_rep ${DB}_${NAME}_rep.fasta
        fi
        if [ ! -s $AA2CDS ]; then
            grep '^>' ${DB}_${NAME}_rep.fasta | sed 's/>//' > $SCRATCH/names.txt
            samtools faidx $CDSDB.gz -r $SCRATCH/names.txt > $AA2CDS
        fi
        # index the db for bwa mapping
        if [[ ! -f $AA2CDS.pac ]]; then
            bwa-mem2 index $AA2CDS
        fi
        GTF=$(echo -n $AA2CDS | perl -p -e 's/\.fasta/.gtf/')
        if [[ ! -f $GTF ]]; then
            python scripts/fasta_to_gtf.py -i $AA2CDS -o $GTF
            bgzip -k $GTF
            tabix $GTF.gz
        fi

        # perhaps cleanup aferwards   
    done
    date
fi

if [ ! -z $DOCDS ]; then
    IN=input_cds
    DB=db/LsFMGC_CDS

    if [ ! -f $DB ]; then
        mmseqs createdb $CDSDB.gz $DB --compressed 1
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


