#!/bin/bash
# TIMEOR_load.sh
# Ashley Mae Conard
# Last Mod. 9/27/19
# Purpose: Run 3 first steps to gather data, run FastQC and align using HISAT2
# NOTE: This script is useful for the command line version of TIMEOR

if [ $# -ne 4 ]; then
	echo $0: "Usage: ./get_fastq_files.sh NUM_THREADS /FASTQ/DIR/ /PATH/TO/SraAccList.txt /OUT/DIR/results/"
	exit 1
fi

# Naming parameters
NUM_THREADS=$1
FASTQ=$2
SRALIST=$3
OUTDIR=$4

# Defining timestamp function
timestamp() {
	date +"%Y-%m-%d_%H-%M-%S"
}

# Getting data
time /data/compbio/timeor/scripts/get_fastq_files.sh ${NUM_THREADS} ${FASTQ} ${SRALIST} &> ${FASTQ}/fastq_file_stderr.txt
wait

# Running fastq quality control
time /data/compbio/timeor/scripts/run_fastQC.sh $FASTQ $OUTDIR/fastQC/ $NUM_THREADS
wait

# Determine file structure for results
mkdir $FASTQ"/case" $FASTQ"/control"
mv `ls $FASTQ/*.fastq* | head -1` $FASTQ"/case"
mv `ls $FASTQ/*.fastq*| head -1` $FASTQ"/control"

# Aligning to genome and producing read counts
time /data/compbio/timeor/scripts/run_HISAT2_HTSeq.sh $FASTQ $OUTDIR _1.fastq.gz 10 $NUM_THREADS 0 genom
wait

echo "TIMEOR finished 1) load data, 2) quality control check, 3) alignment, 4) count matrix creation."
