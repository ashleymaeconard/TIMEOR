#!/bin/bash
# run_HTSeq.sh
# Ashley Mae Conard
# Last Mod: July 11, 2019
# Runs HTSeq for all out.sorted.bam files in a given directory

# Check to make sure input is correct
if [ $# -ne 3 ]; then
	echo $0: "Usage: ./run_HTSeq.sh 
				1) /FULL/PATH/TO/timeor/results/preprocess/alignment/ALIGNMENT_METHOD (hisat2 OR bowtie2)
				2) ORGANISM (dme, hsa, or mmu)
				3) APP_dir (/PATH/TO/TIMEOR/app/)" 
				
	exit 1
fi

# Assign input to variable names
INPUT_DIR=$1
ORGANISM=$2
APP_DIR=$3

# Printing needed directories
ALIGNMENT_DIR=$(echo `dirname $INPUT_DIR`)
LOAD_DIR=$(echo `dirname $ALIGNMENT_DIR`)
RESULTS_DIR=${LOAD_DIR}"/count_matrix/htseq/"
echo "RESULTS_DIR: $RESULTS_DIR"

# Getting number of .bam files files
num_files=$(echo `ls -1q ${INPUT_DIR}/*/*/*out.sorted.bam | wc -l`)
echo "Number of input files: " $num_files

# Determining how many processors and commands to run in one batch
#if [ $((num_files%2)) -eq 0 ]; then # number is even
#    if [[ $num_files -gt  10 ]]; then
#        NUM_PROCESSORS=2
#        NUM_COMMANDS=5
#    else
#        NUM_PROCESSORS=$(echo $((num_files/2)))
#        NUM_COMMANDS=$(echo $((num_files - NUM_PROCESSORS)))
#    fi
#else # number is odd
#    num_filesPlus=$((num_files+1))
#	if [[ $num_files -gt  10 ]]; then
#        NUM_PROCESSORS=2
#        NUM_COMMANDS=5
#    else
#        NUM_PROCESSORS=$(echo `$((num_filesPlus / 2))`)
#        NUM_COMMANDS=$(echo `$((num_filesPlus - NUM_PROCESSORS))`)
#    fi
#fi

NUM_PROCESSORS=1
NUM_COMMANDS=1

if [ -d "$RESULTS_DIR" ]; then
	echo "${RESULTS_DIR} already exists."
else
    mkdir $RESULTS_DIR
fi

# Creating command script for HTSeq
COMMAND_SCRIPT=${RESULTS_DIR}"/run_HTSeq.txt" 

# Removing command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi

# Getting genome    
genome_HTSeq=${APP_DIR}"/../genomes_info/$ORGANISM/genes.gtf"

# Setting the counter for the number of prcesses to run in a batch
START=0
i=$START

echo "Processing .fastq files with HTSeq on $NUM_PROCESSORS processors."
# Adding commands to a script.txt file
for dir in $INPUT_DIR/*/*
    do
    echo "Subdir:" ${dir}
    
    # Iterating through input directory
    if [ -d ${dir} ] ; then
        echo "Creating read counts from aligned reads.";
        for fastq in ${dir}/*out.sorted.bam
            do
                echo "fastq" ${fastq}
                # Incrementing counter
                ((i = i + 1))

                # If number of processors reached, add wait to form a batch that will finish, and then process the next batch
                if (( ${i}%${num_files}==0 )); then # was NUM_PROCESSORS
                    echo "wait" >> $COMMAND_SCRIPT
                fi 
                fileName=$(echo `basename $fastq`) # get filename from out.sorted.bam                                
                replicateFolder=$(echo `basename $(dirname $fastq)`)
                sampleFolder=$(echo `basename $(dirname $(dirname $fastq))`)
                echo "replicateFolder: " $replicateFolder
                echo "sampleFolder: " $sampleFolder
                
                # create folder for all HTSeq outputs per file by removing R1 and R2 (e.g. MTb8-8)
                input_folderName=${INPUT_DIR}/${sampleFolder}/${replicateFolder}
                folderName=${RESULTS_DIR}/${sampleFolder}/${replicateFolder}
                
                # Output directory for /HTSeq/sample_name
                mkdir -p ${folderName}
                
                # write HTSeq command to a file, followed by .bam creation and sorting
                echo "Adding ${fileName} to run_HTSeq.txt script"
                echo " "                    
                echo "(htseq-count -f bam -r pos -i gene_id $input_folderName/*out.sorted.bam $genome_HTSeq > $folderName/htseq_counts) &">> $COMMAND_SCRIPT
        done
    fi
done
echo "wait" >> $COMMAND_SCRIPT

# run command_script (.txt file saved in $RESULTS_DIR)
bash $COMMAND_SCRIPT

