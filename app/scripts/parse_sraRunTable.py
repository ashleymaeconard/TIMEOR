# parse_sraRunTable.py
# Ashley Mae Conard
# Last Mod. 10/23/2019
# Purpose: 1) create timeor/ folder structure
#		 2) create SraAccList.txt
# 		 3) create metadatafile.csv

import sys, os, csv, glob, re
import pandas as pd
import natsort as ns

if len(sys.argv) != 3:
    print("Usage: python parse_metadata.py \n 1)/FULL/PATH/TO/DESIRED/OUTPUT_DIR/ \n 2)/FULL/PATH/TO/SraRunTable.txt")
    sys.exit (1)    

# Arguments
DIR = sys.argv[1]
SRTABLE = sys.argv[2]
OUTPUT_DIR = "/timeor/data/"

sra = "SraAccList.txt" # saved as .csv (even though says .txt)
metadata= "metadata.csv" # saved as .csv

import os
if not os.path.exists(DIR):
    os.makedirs(DIR)

def create_file_structure():
    # Creating output directory trees as needed
    # Creating data stage
    create_folder(DIR,"/timeor/")
    create_folder(DIR,"/timeor/data/")
    create_folder(DIR,"/timeor/data/fastq/")
    create_folder(DIR,"/timeor/data/counts/")

    # Creating results stage
    create_folder(DIR,"/timeor/results/")
        
    # Creating results preprocess stage
    create_folder(DIR,"/timeor/results/preprocess/")
    create_folder(DIR,"/timeor/results/preprocess/fastqc/")
    create_folder(DIR,"/timeor/results/preprocess/alignment/")
    create_folder(DIR,"/timeor/results/preprocess/count_matrix/")
    create_folder(DIR, "/timeor/results/preprocess/count_matrix/htseq/")
    
    # Creating results primary and secondary stage analyses    
    create_folder(DIR,"/timeor/results/analysis/")

def create_metadata_and_accession():
    print("Creating metadata file and SRA accession list.")

    df_to_parse = pd.read_csv(SRTABLE, sep=",")

    # Getting treatment column 
    df_treatm = df_to_parse.filter(regex=re.compile("treatment", re.IGNORECASE))
    if df_treatm.empty:
        sys.exit("One column must have 'treatment' in the column name.")

    # Getting time column
    df_time = df_to_parse.filter(regex=re.compile("time", re.IGNORECASE))
    if df_time.empty:
        sys.exit("One column must have 'time' in the column name.")

    # Getting replicate column
    df_rep = df_to_parse.filter(regex=re.compile("replicate", re.IGNORECASE))
    if df_rep.empty:
        sys.exit("One column must have 'replicate' in the column name.")

    # Getting Run ID column
    df_run = df_to_parse.filter(regex=re.compile("run", re.IGNORECASE))
    if df_run.empty:
        sys.exit("One column must have 'run' in the column name.")

    # Getting batch column
    df_batch = df_to_parse.filter(regex=re.compile("batch", re.IGNORECASE))
    if df_batch.empty:
        sys.exit("One column must have 'batch' within the column name.")
            
    # Concatenate and rename necessary metadata columns
    sraRunTable_df = pd.concat([df_treatm, df_time, df_rep, df_batch, df_run], axis=1)
    sraRunTable_df.columns = ['treatment','time','replicate','batch','Run']

    # Casting potentially numerical columns as string 
    sraRunTable_df['time'] = sraRunTable_df['time'].astype(str)
    sraRunTable_df['replicate'] = sraRunTable_df['replicate'].astype(str)
    sraRunTable_df['batch'] = sraRunTable_df['batch'].astype(str)
    
    # Creating the unique ID (treatment_time.replicate)
    sraRunTable_df['ID'] = sraRunTable_df.treatment.str.split().str.get(0).str.cat(sraRunTable_df.time.str.split().str.get(0), sep ="_").str.cat(sraRunTable_df.replicate.str.split().str.get(0), sep =".").str.cat(sraRunTable_df.batch.str.split().str.get(0), sep =".")
    
    # Change column order
    sraRunTable_df = sraRunTable_df.reindex(columns=['ID','treatment','time','replicate','batch','Run'])
    
    # Sort by time
    sraRunTable_df['time'] = pd.Categorical(sraRunTable_df['time'], ordered=True, categories= ns.natsorted(sraRunTable_df['time'].unique()))
    sraRunTable_df = sraRunTable_df.sort_values(by=['time', 'ID'])

    # Save SRA accession list and metadata files.
    sraRunTable_df.Run.to_csv(DIR+"/timeor/data/SraAccList.csv", index=False)
    sraRunTable_df.to_csv(DIR+"/timeor/data/metadata.csv", index=False)
    print("Saved metadata file and SRA accession list here:",DIR+"/timeor/data/SraAccList.csv")
    return(sraRunTable_df)

def create_folder(FULL_DIR, TO_ADD):
    if not os.path.exists(FULL_DIR+TO_ADD):
        os.makedirs(FULL_DIR+TO_ADD)

def main(argv):

    # Creating file structure
    create_file_structure()
    
    # Creating metadata and accession list files
    metadataDf = create_metadata_and_accession()

# Calling main function
if __name__ == "__main__":
    main(sys.argv[1:])
