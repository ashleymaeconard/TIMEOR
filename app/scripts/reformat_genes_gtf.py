# reformat_genes_gtf.py
# Ashley Mae Conard
# Last Mod. 2/3/2020
# Purpose: Reformats genes.gtf and creates reformatted_genes_gtf.csv with cols: 
#           index	chrom	type_prot	start_chrom	end_chrom	gene_name	
#           gene_biotype	gene_id	transcript_name	transcript_id	tss_id 

# Checking input arguments
import sys
if len(sys.argv) != 3:
    print("Type: python reformat_genes_gtf.py \n"+
            "\t 1) /PATH/TO/genes.gtf (e.g. ~/timeor-shiny-app-v1/genomes_info/dm6/genes.gtf)\n"+ 
            "\t 2) /PATH/TO/OUTPUT_DIR/")
    sys.exit (1)

# Importing libraries
import os, math, re, glob
import pandas as pd

# Importing arguments
GENES_GTF = sys.argv[1] # e.g. "/data/compbio/aconard/BDGP6/genes.gtf"
OUTPUT_DIR = sys.argv[2]

# Splitting allInfo in genes.gtf file 
def create_cols(row, sw):
    string = row['allInfo']
    list_of_words = string.split()
    next_word = list_of_words[list_of_words.index(sw) + 1]   
    return(next_word)

# Reformatting genes.gtf file
def process_genes_gtf(genes_gtf):
    # Importing genes.gtf file 
    df_genes_gtf = pd.read_csv(genes_gtf, sep="\t",usecols=[0,2,3,4,8], names=['chrom','type_prot', 'start_chrom', 'end_chrom', 'allInfo'], engine='python')

    # Parsing GTF for relevant information (gene_name, gene_biotype, ..., tss_id)
    df_genes_gtf['allInfo'] = df_genes_gtf['allInfo'].map(lambda x: str(x)[:-1])
    df_genes_gtf.columns.str.replace(' ', '')
    list_cols =["gene_name", "gene_biotype", "gene_id"] #, "transcript_name", "transcript_id", "tss_id"] # gene_biotype = gene_type, NO TSS_ID
    print("Example of 'allInfo' column in df_genes_gtf['allInfo'] that I will parse into these columns: \n\n", 
        list_cols,"\n\n", df_genes_gtf["allInfo"][0])

    # Iterating through all gene information (from allInfo in df_genes_gtf of genes.gtf) to keep
    for search_word in list_cols:
        df_genes_gtf[search_word] = df_genes_gtf.apply(create_cols,args=[search_word],axis=1)
        df_genes_gtf[search_word] = df_genes_gtf[search_word].str.replace(';','')
        df_genes_gtf[search_word] = df_genes_gtf[search_word].str.replace('"', '')
        
    # Setting index to gene_id
    df_genes_gtf = df_genes_gtf.drop(columns=['allInfo'])
    df_genes_gtf = df_genes_gtf.drop_duplicates(subset="gene_id")
    df_genes_gtf = df_genes_gtf.reset_index()
    print("Number of unique genes based on FlyBase ID from ", GENES_GTF,": ", len(df_genes_gtf))
    return(df_genes_gtf)

def main(argv):
    # Processing .gtf file into usable format
    df1_genes_gtf = process_genes_gtf(GENES_GTF)
    df1_genes_gtf.to_csv(OUTPUT_DIR+"/reformatted_genes_gtf.csv", index=False)

if __name__ == "__main__": # in Python 3 we need not use this, instead just main() but this is more universal.
    main(sys.argv[1:])