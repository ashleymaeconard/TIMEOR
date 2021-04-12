# fastq_into_folders.py
# Ashley Mae Conard
# Last Mod. 11/23/2019
# Purpose: Creates sample and replicate subfolders and moves fastq.gz files into those subfolders.

# Importing libraries
import sys, os, csv, glob
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: python fastq_into_folders.sh /FULL/PATH/TO/FASTQ_FILES/ /FULL/PATH/TO/METADATA.txt")
    sys.exit (1)    

# Arguments
FASTQ_DIR = sys.argv[1]
METADATA = sys.argv[2]

def create_folder(FULL_DIR, TO_ADD):
    if not os.path.exists(FULL_DIR+TO_ADD):
        os.makedirs(FULL_DIR+TO_ADD)

def create_subfolders(metadata_df):
    print(metadata_df)
    for i, j in metadata_df.iterrows():
        samples, replicates = (j["ID"].split(".")[0], j["Run"])
        subdir = samples+"/"+replicates
        create_folder(FASTQ_DIR, subdir)
        move_fastqs(FASTQ_DIR, samples, replicates)

def move_fastqs(DIR, samps, reps):
    fastq_star = (DIR+"/*"+reps+"*")
    for i in glob.glob(fastq_star):
        fastq_file = os.path.basename(i)
        fastq_move_from = (DIR+"/"+fastq_file)
        fastq_move_to = (DIR+"/"+samps+"/"+reps+"/"+fastq_file)
        os.rename(fastq_move_from, fastq_move_to)

def main(argv):
    # Reading in metadata file
    metadataDf = pd.read_csv(METADATA, sep=",")
    #print(metadataDf)
    
    # Creating subfolders for sample and replicates
    create_subfolders(metadataDf)

    print("Moved raw .fastq data into sample and replicate folders here:")
    print(FASTQ_DIR)

# Call main function
if __name__ == "__main__":
    main(sys.argv[1:])