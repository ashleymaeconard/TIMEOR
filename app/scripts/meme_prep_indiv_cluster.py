# meme_prep.py
# Ashley Mae Conard
# Last Mod. 7/7/2019
# Purpose: Creates MEME input list of DNA sequences for a set of genes

# Checking for arguments
import sys
if len(sys.argv) != 6:
    print("Type: python meme_prep.py \n"+
                "\t 1) /PATH/TO/reformatted_genes_gtf.csv (e.g. ~/Desktop/timeor_example/reformatted_genes_gtf.csv)\n"+
                "\t 2) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/a_results/<NUM_CLUST>/)\n"+
                "\t 3) /PATH/TO/CHROM_FAs (e.g. ~/timeor-shiny-app-v1/genomes_info/dme/dm6.fa)\n"+
                "\t 4) TSS_only (set to 1 to run +-1kb from transcription start site, 0 otherwise)\n"+
                "\t 5) GENOME (dme, hsa, or mmu)\n")
    sys.exit (1)    

# Inputting arguments
GENES_GTF_DM6 = sys.argv[1] 
OUTPUT_DIR = sys.argv[2] 
CHROMS = sys.argv[3] 
TSS_only = int(sys.argv[4])
GENOME = sys.argv[5] 

# Importing libraries
import os, math, re, glob, subprocess
import pandas as pd

# Choosing chromosome numbers:
if GENOME =="dme":
    list_chroms = ["2L", "2R", "3L", "3R", "4", "X", "Y"]
elif GENOME == "mmu":
    list_chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"]
elif GENOME == "hsa":
    list_chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
else:
    sys.exit("ERROR: must input genome type mus, hse, or dme.")

def main(argv):
    # Processing .gtf file into usable format
    df1_genes_gtf = pd.read_csv(GENES_GTF_DM6)

    # Form .txt files for input into MEME for each gene_cluster within each case/control pair
    if TSS_only:
        print("Processing",OUTPUT_DIR," into gene_name, gene_id, chrom, TSS-1000, TSS+1000 files.\n")
    else:
        print("Processing",OUTPUT_DIR," into gene_name, gene_id, chrom, start-1000, end+1000 files.\n")
    col_names = ["gene_name","gene_id","chrom","start_chrom","end_chrom"] 

    # For each gene_list (cluster)
    for i in glob.glob(OUTPUT_DIR+"/*geneList*"):
        geneClusterID = i.split("/")
        geneClusterNum = geneClusterID[-1].strip(".csv").split("_")[-1] # gene list (i.e. cluster) number 
        print("Genelist:", geneClusterID[-1].strip(".csv"), "and gene number: ", geneClusterNum)
        
        # Getting all gene names and set column names
        df = pd.DataFrame(columns = col_names)
        df_geneList = pd.read_csv(i)
        df_geneList = df_geneList.rename(columns={"Unnamed: 0": "gene_name"})
        gene_list = df_geneList["gene_name"].tolist() # create list of all gene names in gene list (i.e. cluster)
        print("Number genes in this cluster:", len(gene_list))
        
        # For each gene in the gene list (cluster), create a df for each gene: "gene_name gene_id chrom start_min1000 end_min1000"
        counter=0
        for g in gene_list:# gene in gene_list's gene_name, id, chrom, start, end
            gene_info_from_gtf = df1_genes_gtf.loc[df1_genes_gtf['gene_name']==g, ["gene_name","gene_id","chrom",'start_chrom','end_chrom']]
            if gene_info_from_gtf.empty: # no gene_name in gtf file
                counter+=1
            elif TSS_only: # if you should only look around TSS (+-1KB from start of gene)
                gene_info_1 = gene_info_from_gtf.to_string(index=False)
                gene_info_2 = re.sub(r'.*\n', '', gene_info_1)
                gene_info_3 = list(filter(None, gene_info_2.split(' ')))
                df.loc[len(df)] = [gene_info_3[0],gene_info_3[1],gene_info_3[2], int(gene_info_3[3])-1000, int(gene_info_3[3])+1000]# "gene_name","gene_id","chrom",'start_chrom-1KB','START_CHROM+1KB'             
            else: # look at entire gene
                gene_info_1 = gene_info_from_gtf.to_string(index=False)
                gene_info_2 = re.sub(r'.*\n', '', gene_info_1)
                gene_info_3 = list(filter(None, gene_info_2.split(' ')))
                df.loc[len(df)] = [gene_info_3[0],gene_info_3[1],gene_info_3[2], int(gene_info_3[3])-1000, int(gene_info_3[4])+1000] # "gene_name","gene_id","chrom",'start_chrom-1KB','end_chrom+1KB'
        
        # Creating bed file for bedtools getfastq 
        bed_file = df[['chrom','start_chrom','end_chrom','gene_name']].copy()
        
        # Creating output directories
        outdir_meme_input = (OUTPUT_DIR+"/MEME/")
        if not os.path.exists(outdir_meme_input): 
            os.makedirs(outdir_meme_input)
            print("Created: ", outdir_meme_input, "\n")
        
        # Saving bed and DNA sequence fasta file
        if TSS_only:# (TSS -/+1KB)
            # Saving TSS bed file in MEME subdirectory for each cluster in each experiment
            bed_file.to_csv(outdir_meme_input+geneClusterID[-1].split(".",1)[0]+"_TSS.bed", sep='\t', index=False, header=False)
            # Getting TSS bed file location for bedtools call
            bed_file_loc = outdir_meme_input+geneClusterID[-1].split(".",1)[0]+"_TSS.bed"
            # Creating sequence file output to pass to run_meme.txt
            gene_seqs_out_file = outdir_meme_input+geneClusterID[-1].split(".",1)[0]+"_TSS_DNAseqs.fasta"
        else:
            # Saving bed file in MEME subdirectory for each cluster in each experiment
            bed_file.to_csv(outdir_meme_input+geneClusterID[-1].split(".",1)[0]+".bed", sep='\t', index=False, header=False)
            # Getting bed file location for bedtools call
            bed_file_loc = outdir_meme_input+geneClusterID[-1].split(".",1)[0]+".bed"
            # Creating sequence file output to pass to run_meme.txt
            gene_seqs_out_file = outdir_meme_input+geneClusterID[-1].split(".",1)[0]+"_DNAseqs.fasta"                
        
        # Calling bedtools to produce fasta file of DNA sequences for a given list of genes
        subprocess.call(["/bedtools2/bin/./bedtools", "getfasta", "-fi", CHROMS, "-bed", bed_file_loc, "-name", "-fo",gene_seqs_out_file])
        
if __name__ == "__main__": # in Python 3 we need not use this, instead just main() but this is more universal.
    main(sys.argv[1:])
