#!/bin/bash
# get_fastq_files.sh
# Ashley Mae Conard
# Last Mod. 7/9/2019
# Resource: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129292

if [ $# -ne 2 ]; then
	echo $0: "Usage: ./get_fastq_files.sh 
			1) /FULL/PATH/TO/OUTPUT/FASTQ_DIR/ \(to place fastq.gz files\) 
			2) /PATH/TO/SraAccList.txt"
	exit 1
fi

# Assigning input parameters to variable names
OUTPUT_DIR=$1
SRA_LIST=$2

# Defining timestamp function
timestamp() {
	date +"%Y-%m-%d_%H-%M-%S"
}

COMMAND_SCRIPT=${OUTPUT_DIR}"/run_get_fastq_files.txt" 

# Removing command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

# Creating accession number list
while IFS=\= read var value; do
    accessionList+=($var)
done < ${2}

# Getting number of fastq files
num_files=${#accessionList[@]}

NUM_PROCESSORS=1
# Determining how many processors
#if [[ $num_files -gt  10 ]]; then
#	NUM_PROCESSORS=10
#else
#	NUM_PROCESSORS=$num_files
#fi

# Creating output directory if needed
mkdir -p ${OUTPUT_DIR}

# Getting raw .fastq files from NCBI GEO
#echo "Starting to download "${accessionList[@]}" at $(timestamp)".
echo "(time parallel -j ${NUM_PROCESSORS} --bar fastq-dump {1} --split-files -O ${OUTPUT_DIR} --gzip ::: ${accessionList[@]}) &" >> $COMMAND_SCRIPT  
#echo "Finished data download at $(timestamp)"

# Run command script
bash ${COMMAND_SCRIPT}
