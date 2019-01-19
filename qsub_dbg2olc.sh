#!/bin/bash
#$ -N dbg13o30k2a.01km25
#$ -pe openmp 2
#$ -R Y
#$ -q bio

#module load blasr/20140724
module load jje/jjeutils/0.1a
module load mchakrab/dbg2olc
module load smrtanalysis/2.3.0p5 
module load canu/2016-01-13

##Variables for editing based on assembly
#platanus contigs
CONTIGFILE="iso1.pool_contig.fa"
#raw reads long reads
FULLPBREADS="chardonnay.fq"
#kmer value from canu
MYKVAL=25
#expected genome size
GENOMESIZE=130000000
#desired coverage 30 is a good value
COVERAGE=30
###Do not edit below
BACKBONERAWFA="backbone_raw.fasta"
DBG2OLCCONS="DBG2OLC_Consensus_info.txt"
PACBIOREADS="${PREFIX}_${COVERAGE}x.u.fastq"
#test for fq or fastq file
if [ $(echo ${FULLPBREADS} | awk -F . '{print $NF}') = "fastq" ]; then
   PREFIX=$(basename ${FULLPBREADS} .fastq)
   ln -sf ${FULLPBREADS} ${PREFIX}.u.fastq
elif [ $(echo ${FULLPBREADS} | awk -F . '{print $NF}') = "fq" ]; then
   PREFIX=$(basename ${FULLPBREADS} .fq)
   ln -sf ${FULLPBREADS} ${PREFIX}.u.fastq
else
   echo "please give a fastq input file. run gunzip ${FULLPBREADS} if necessary"
   exit 1
fi

FULLPBREADS=${PREFIX}.u.fastq
#test for coverage
echo "calculate top ${COVERAGE}x for genomesize: ${GENOMESIZE} for ${PREFIX}"
MYCOVERAGE=$(($(bioawk -cfastx '{sum+=length($seq)} END {print sum}' $FULLPBREADS)/$GENOMESIZE))
if [ $MYCOVERAGE -gt $COVERAGE ]; then
   echo "coverage less that ${COVERAGE}. Coverage = ${MYCOVERAGE}"
   ln -sf ${PREFIX}.u.fastq ${PACBIOREADS}
else
   echo "Begin calculating longest reads for ${COVERAGE}x coverage"
   fastqSample -I ${PREFIX} -U -O ${PREFIX}_${COVERAGE}x -max -g ${GENOMESIZE} -c ${COVERAGE}
   echo "end calculate"
fi
bioawk -cfastx '{print ">"$name"\n"$seq}' ${PACBIOREADS} > $(basename ${PACBIOREADS} .fastq).fasta
PACBIOREADS="${PREFIX}_${COVERAGE}x.u.fasta"
#to here

###overlap
###using contig file as input example
#DBG2OLC k KmerSize AdaptiveTh THRESH_VALUE1 KmerCovTh THRESH_VALUE2 MinOverlap THRESH_VALUE3 \
#   Contigs NGS_CONTIG_FILE \
#   f LONG_READS.FASTA RemoveChimera 1

###parameter tuning
##For 10x/20x PacBio data:
#KmerCovTh 2-5, MinOverlap 10-30, AdaptiveTh 0.001~0.01
##For 50x-100x PacBio data:
#KmerCovTh 2-10, MinOverlap 50-150, AdaptiveTh 0.01-0.02

MYLD=0
MYKMERCOV=2
MYADAPT=0.01
MYMINOVL=35
MYRMCH=1

echo "begin dbg2olc step 1 overlap"
#DBG2OLC_Linux k 17 AdaptiveTh 0.001 KmerCovTh 2 MinOverlap 10 RemoveChimera 1 \
#   Contigs ${CONTIGFILE} \
#   f ${PACBIOREADS}
#
DBG2OLC_Linux LD ${MYLD} k ${MYKVAL} AdaptiveTh ${MYADAPT} KmerCovTh ${MYKMERCOV} \
   MinOverlap ${MYMINOVL} RemoveChimera ${MYRMCH} \
   Contigs ${CONTIGFILE} \
   f ${PACBIOREADS}
echo "end dbg2olc step 1 overlap"

###consensus
#check n50 of backbone file first; proceed if ok
#concat pbreads and contigs

cat ${CONTIGFILE} ${PACBIOREADS} > ctg_pb.fasta

#ulimit -n unlimited

#consensus script SPARC
#split_and_run_sparc.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt \
#   ctg_pb.fasta \
#   ./consensus_dir 2 > cns_log.txt

#consensus script PBDAGCON
echo "begin dbg2olc step 2 consensus"
split_and_run_pbdagcon.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt \
    ctg_pb.fasta \
    ./consensus_dir 2>cns_log.err 1>cns_log.out
echo "end dbg2olc step 2 consensus"

#test if consensus fasta file exists
if [ -f consensus_dir/final_assembly.fasta ]
then
   cp consensus_dir/final_assembly.fasta ./${PREFIX}_${COVERAGE}x_LD${MYLD}_K${MYKVAL}_KCOV${MYKMERCOV}_ADAPT${MYADAPT}_MINOVL${MYMINOVL}_RMCHIM${MYRMCH}_asm.fasta
else
   echo "Consensus failed. consensus_dir/final_assembly.fasta does not exist. Please check for errors and rerun script"
   exit 1
fi
