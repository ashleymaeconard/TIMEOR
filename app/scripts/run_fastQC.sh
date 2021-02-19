# run_fastQC.sh
# Ashley Mae Conard
# Last Mod. 10/15/2019
# Purpose: Runs FastQC in /fastqc folder for all .fastq.gz files in a given directory

# Check to make sure input is correct
if [ $# -ne 2 ]; then
	echo $0: "Usage: ./run_fastQC.sh 
			1) /FULL/PATH/TO/FASTQ_DIR/ \(containing fastq.gz files\) 
			2) /FULL/PATH/TO/OUTPUT_DIR/ \(e.g. timeor/results/load/fastqc/\)"
	exit 1
fi

# Getting number of fastq files
num_files=$(echo `ls -1q ${1}/*.fastq* | wc -l`)
 
# Determining how many processors
#if [[ $num_files -gt  10 ]]; then
#	NUM_PROCESSORS=10
#else
#	NUM_PROCESSORS=$num_files
#fi
NUM_PROCESSORS=1

# Checking to see if output directory exists
mkdir -p $2

# Get all fastq files in an array format
declare -a arr
for file in ${1}/*.fastq*; do arr=(${arr[@]} "$file"); done
echo "FastQC on $num_files files: ${arr[@]}"
echo ${2}/run_fastQC.txt
echo "time parallel -j $NUM_PROCESSORS fastqc {1} --outdir=${2} ::: ${arr[@]}" > ${2}/run_fastQC.txt

# Summarize all FastQC results with MultiQC
echo "eval '$(conda shell.bash hook)'" >> ${2}/run_fastQC.txt
echo "conda deactivate" >> ${2}/run_fastQC.txt
echo "conda activate timeor_py38" >> ${2}/run_fastQC.txt
echo "cd ${2}" >> ${2}/run_fastQC.txt
echo "multiqc ." >> ${2}/run_fastQC.txt
echo "conda deactivate" >> ${2}/run_fastQC.txt
echo "conda activate timeor_conda_env" >> ${2}/run_fastQC.txt
echo "wait" >> $COMMAND_SCRIPT
# Run command script
bash ${2}/run_fastQC.txt
