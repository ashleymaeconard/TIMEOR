#!/usr/bin/python
# htseq_merge.py
# Ashley Conard
# Last Mod. 7/17/2019
# Purpose: Merges all HTSeq files into one table of genes by conditions
#          Assumes structure of .../project/results/hisat2_htseq/SAMPLE_NAME/REPLICATE/htseq_counts

# Importing libraries
import sys, os, csv, glob
import pandas as pd

if len(sys.argv) != 2 :
    print("Usage: python htseq_merge.py /FULL/PATH/TO/htseq/ (with subfolders e.g. ../htseq/SAMPLE/REP/htseq_counts stop at htseq/)")
    sys.exit (1)    

# Arguments
DIR=sys.argv[1]
OUTPUT_FILE = "htseq_merged" # saved as .csv
OUTPUT_DIR = DIR

#Check and/or make output directory
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def main(argv):
    first=1
    # Iterate through all samples
    for sample in glob.glob(DIR+"*/"):
        replicate_id = 0 # setting replicate id to be set per sample
        rep_name=os.path.basename(os.path.dirname(sample)) # getting replicate name
        print("Sample: ", sample)
        
        # Iterate through all replicates
        for rep in glob.glob(sample+"/*"):
            replicate_id+=1
            file_loc = rep+"/htseq_counts"
            col_name = ".".join((rep_name,str(replicate_id))) # creating replicate name with id
            print("Replicate: ", file_loc)
            
            # Creating initial dataframe
            if first:
                merged_htseq_df = pd.read_csv(file_loc,  sep='\t', names=["gene_id", col_name]) # creating merged_htseq_df
                merged_htseq_df = merged_htseq_df[~merged_htseq_df.gene_id.str.contains("_")]
                first=0
                len_df = len(merged_htseq_df.index) # checking size of dataframe
                print("Size: ", len_df)
                
            # Adding to merged dataframe
            else:
                df = pd.read_csv(file_loc,  sep='\t', names=["gene_id", col_name])
                merged_htseq_df = pd.merge(merged_htseq_df, df, on='gene_id')
                print("Size of df to merge: ", df.shape[0])
                print("Adding rep. Size: ", merged_htseq_df.shape[0])
    
    # Check that dataframe length is maintained when adding more replicates htseq files
    if len_df!=len(merged_htseq_df.index):
        sys.exit("ERROR - merged dataframe size changed during merge step.")
    
    merged_htseq_df = merged_htseq_df.reindex(sorted(merged_htseq_df.columns), axis=1)
    cols = list(merged_htseq_df)
    cols.insert(0, cols.pop(cols.index('gene_id')))
    merged_htseq_df = merged_htseq_df.ix[:, cols]
    
    merged_htseq_df.to_csv(OUTPUT_DIR+"/merged_htseq.csv", index=False)
    print("Placing merged_htseq.csv in ", OUTPUT_DIR)
    
if __name__ == "__main__":
    main(sys.argv[1:])
