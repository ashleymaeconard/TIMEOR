#!/bin/bash
# run_HISAT2.sh
# Ashley Mae Conard
# Last Mod. 7/11/2019
# Purpose: Runs HISAT2 in for all .fastq.gz files in a given directory

# Checking input arguments
if [ $# -ne 5 ]; then
	echo $0: "Usage: ./run_HISAT2.sh 
		1) /PATH/TO/NAMED_FASTQ_DIR/ (NAMED folders containing replicate fastq.gz files) 
		2) /PATH/TO/RESULTS_DIR/ (desired location to add or create results directory) 
		3) PAIRED_OR_NOT (YES is 1, NO is 0) 
		4) APP_DIR
		5) ORGANISM (dme, hsa, mmu)"
	exit 1
fi

# Arguments
INPUT_DIR=$1
RESULTS_DIR=$2
PAIRED_OR_NOT=$3
APP_DIR=$4
ORGANISM=$5

# Checking to see if results directory exists
if [ -d "$RESULTS_DIR" ]; then
	echo "${RESULTS_DIR} already exists."
else
    mkdir $RESULTS_DIR
fi

# Getting number of fastq files
num_files=$(echo `ls -1q ${INPUT_DIR}/*/*/*.fastq* | wc -l`)
echo "Number of input .fastq files: " $num_files

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

NUM_PROCESSORS=4 # -p flag for hisat2
NUM_COMMANDS=1 # number of hisat2 commands

# Checking to see if results/HISAT2/ directory exists
#if [ -d "$RESULTS_DIR/hisat2/" ]; then
#	echo "${RESULTS_DIR}/hisat2/ already exists."
#    OUTPUT_DIR=${RESULTS_DIR}"/hisat2/"
#else
#    OUTPUT_DIR=${RESULTS_DIR}"/hisat2/"
#    mkdir $OUTPUT_DIR
#fi

# Creating command script for HISAT2
COMMAND_SCRIPT=${RESULTS_DIR}"/run_HISAT2.txt" 

# Removing command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
	echo -e "REPLACING ${COMMAND_SCRIPT}\n"
	rm -rf ${COMMAND_SCRIPT}
fi
    
# Getting genome
genome="/srv/genomes_info/${ORGANISM}/genome_hisat2/genome"
                          
START=0
i=$START

echo "Processing .fastq files with HISAT2 on $NUM_PROCESSORS processors."

for dir in $INPUT_DIR/*/*
    do
    echo " "
    echo "Subdir:" ${dir}
    
    # Iterating through input directory
    if [ -d ${dir} ] ; then
        
        # Checking if RNA-seq data is paired or not
        # SINGLE END READS
        if [ "$PAIRED_OR_NOT" -ne "1" ]; then
           echo "Initiating single-end RNA-seq data processing.";
           for fastq in ${dir}/*
               do
                   echo "Getting unpaired .fastq for $fastq"
                   # Iterating through each R1 to determine when to add 'wait' to script
                   ((i = i + 1))

                    # If number of processors reached, add wait to form a batch that will finish, and then process the next batch
                    # if (( ${i}%${num_files}==0 )); then # had been NUM_COMMANDS
                    #     echo "wait" >> $COMMAND_SCRIPT
                    # fi 

                    # Generating the HISAT2 command
                    fileName=$(echo `basename $fastq`) # get filename from .fq.gz
                    replicateFolder=$(echo `basename $(dirname $fastq)`)
                    sampleFolder=$(echo `basename $(dirname $(dirname $fastq))`)
                    
                    # Creating folder for all HISAT2 outputs per fastq
                    folderName=${RESULTS_DIR}/${sampleFolder}/${replicateFolder}
                    
                    # Outputting directory for HISAT2/sample_name
                    mkdir -p ${folderName}
                    
                    # Writing HISAT2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} to run_HISAT2.txt script"
                    echo " "                    
                    #echo "hisat2 -p ${NUM_PROCESSORS} --dta -x $genome -U $fastq -t --un-gz $folderName/out_unalign.gz --al-gz $folderName/out_al_atleastonce.gz --novel-splicesite-outfile $folderName/novel_splicesite --summary-file $folderName/summaryfile.txt --met-file $folderName/met-file.txt -S $folderName/out.sam 2> $folderName/alignmentsummary.txt; samtools sort -O sam -T $folderName/out.sorted -o $folderName/out.sorted.sam $folderName/out.sam; samtools view -S -b $folderName/out.sorted.sam > $folderName/out.sorted.bam; rm -rf $folderName/out.sam $folderName/out.sorted.sam">> $COMMAND_SCRIPT
                    echo "/hisat2-0f01dc6397a/./hisat2 -p ${NUM_PROCESSORS} --dta -x $genome -U $fastq -t --novel-splicesite-outfile $folderName/novel_splicesite --summary-file $folderName/summaryfile.txt -S $folderName/out.sam 2> $folderName/alignmentsummary.txt; samtools sort -O sam -T $folderName/out.sorted -o $folderName/out.sorted.sam $folderName/out.sam; samtools view -S -b $folderName/out.sorted.sam > $folderName/out.sorted.bam; rm -rf $folderName/out.sam $folderName/out.sorted.sam">> $COMMAND_SCRIPT
                    echo "Outputting results to $folderName"
           done
        
        # PAIRED END READS
        else
            echo "Initiating paired-end RNA-seq data processing.";
            # Iterating through only the R1 replicates
            for R1 in ${dir}/*R1*
                do
                    echo "Getting paired .fastq for $R1"
                    # Iterating through each R1 to determine when to add 'wait' to script
                    ((i = i + 1))

                    # If number of processors reached, add wait to form a batch that will finish, and then process the next batch
                    # if (( ${i}%${num_files}==0 )); then # had been NUM_COMMANDS
                    #     echo "wait" >> $COMMAND_SCRIPT
                    # fi 

                    # Generating the HISAT2 command for R1 and R2
                    fileName=$(echo `basename $R1`) # get filename from .fq.gz
                    replicateFolder=$(echo `basename $(dirname $fastq)`)
                    sampleFolder=$(echo `basename $(dirname $(dirname $fastq))`)
                    
                    # Creating folder for all HISAT2 outputs per fastq
                    folderName=${RESULTS_DIR}/${sampleFolder}/${replicateFolder}

                    # Getting 2nd read pair
                    R2=${R1//R1/R2}

                    # Outputting directory for HISAT2/sample_name
                    mkdir -p ${folderName}

                    # Writing HISAT2 command to a file, followed by .bam creation and sorting
                    echo "Adding ${fileName} to run_HISAT2.txt script"
                    echo " "
                    #echo "hisat2 -p ${NUM_PROCESSORS} --dta -x $genome -1 $R1 -2 $R2 -t --un-conc-gz $folderName/out_unconc.gz --al-conc-gz $folderName/out_al_conc_atleastonce.gz --novel-splicesite-outfile $folderName/novel_splicesite --summary-file $folderName/summaryfile.txt --met-file $folderName/met-file.txt -S $folderName/out.sam 2> $folderName/alignmentsummary.txt; samtools sort -O sam -T $folderName/out.sorted -o $folderName/out.sorted.sam $folderName/out.sam; samtools view -S -b $folderName/out.sorted.sam > $folderName/out.sorted.bam; rm -rf $folderName/out.sam $folderName/out.sorted.sam">> $COMMAND_SCRIPT            
                    echo "/hisat2-0f01dc6397a/./hisat2 -p ${NUM_PROCESSORS} --dta -x $genome -1 $R1 -2 $R2 -t --novel-splicesite-outfile $folderName/novel_splicesite --summary-file $folderName/summaryfile.txt -S $folderName/out.sam 2> $folderName/alignmentsummary.txt; samtools sort -O sam -T $folderName/out.sorted -o $folderName/out.sorted.sam $folderName/out.sam; samtools view -S -b $folderName/out.sorted.sam > $folderName/out.sorted.bam; rm -rf $folderName/out.sam $folderName/out.sorted.sam">> $COMMAND_SCRIPT            
                    echo "Outputting results to $folderName"
            done
       fi
    fi
done
#echo "wait" >> $COMMAND_SCRIPT

# Running command_script (.txt file saved in $RESULTS_DIR)
bash ${COMMAND_SCRIPT}
