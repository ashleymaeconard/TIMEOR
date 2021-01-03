#!/bin/bash
# run_meme.sh
# Ashley Mae Conard
# Last Mod. 2/25/2020
# Runs MEME to identify motifs de novo for a set of genes

# Checking input arguments
if [ $# -ne 1 ]; then
	echo $0: "Usage: ./run_meme_indiv_cluster.sh 
                1) /PATH/TO/INPUT_DIR/ e.g. /timeor/results/primary/a_results/<NUM_CLUST>/ 
                (output is placed in /results/primary/CONDITION_results/clusters/<NUM_CLUST>/MEME/)"
	exit 1
fi

# Inputting arguments
INPUT_DIR=$1 

for i in $INPUT_DIR
    do 
    echo ${i}

   for geneList in $INPUT_DIR/MEME/*DNAseqs* # for each DNA sequence file within each cluster
   do
       echo $geneList

        # Setting output directory
       OUTPUT_DIR=$(dirname "${geneList}")
       # -objfun classic --revcomp
       (meme $geneList -dna -nmotifs 3 -maxsize 150000 -mod anr -oc $OUTPUT_DIR/ ) &	

   done
   
    # Waiting for parallel processes to finish before proceeding next experiment.
   for job in `jobs -p`
   do
       echo $job
       wait $job || let "FAIL+=1"
   done
   if [ "$FAIL" == "0" ];
   then
       echo "Yes, finished parallel processes."
   else
       echo "Did not finish parallel processes yet! ($FAIL)"
   fi
done