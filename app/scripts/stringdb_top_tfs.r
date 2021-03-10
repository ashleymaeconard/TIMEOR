# stringdb_top_tfs.r
# Ashley Mae Conard
# Last Mod. 1/9/2020
# Purpose: Finds the network (if any) connecting genes, given a list of genes.

# Checking input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type: /usr/local/bin/Rscript stringdb.r 
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/) 
            2) ORGANISM_NCBI_ID (# NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
            3) TFs_NOT_CLUSTERS (change to 1 if TFs, 0 otherwise (i.e. go through clusters)
            4) /APP/DIR/ (/PATH/TO/TIMEOR_APP/app/)", call.=FALSE)
} else if (length(args) == 4) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/) 
            2) ORGANISM_NCBI_ID (# NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
            3) TFs_NOT_CLUSTERS (change to 1 if TFs, 0 otherwise (i.e. go through clusters)
            4) /APP/DIR/ (/PATH/TO/TIMEOR_APP/app/)")
}

# Inputting Arguments
IN_OUTPUT <- args[1]
NCBI_TAXO <- strtoi(args[2]) # NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
TFs_NOT_CLUSTERS <- as.integer(args[3])
APP_DIR <- args[4]

# Installing libraries
library("STRINGdb")
library(dplyr)
library(stringr)

add_to_info_table <- function(gene_inters, str_db, pr_table, pval){
    # Creating new columns "A_gene_name" and "B_gene_name"
    gene_inters[,"A_gene_name"] <- NA
    gene_inters[,"B_gene_name"] <- NA
    gene_inters[,"pubmed_ids"] <- NA
    gene_inters <- gene_inters %>% select("B_gene_name", everything()) # moving "B" column to the front of gene_inters df
    gene_inters <- gene_inters %>% select("A_gene_name", everything()) # moving "A" column to the front of gene_inters df
    
    for (row in 1:nrow(gene_inters)){
        # Adding names to each row in new columns "A_gene_name" and "B_gene_name"
        gene_id_B <- gene_inters[row, "to"]
        #print(paste("gene_id_B", gene_id_B, sep=" "))
        gene_inters[row, "B_gene_name"] <- pr_table$preferred_name[pr_table$protein_external_id == gene_id_B]
        #print(paste("gene_id_B",gene_inters))
        gene_id_A <- gene_inters[row, "from"]
        gene_inters[row, "A_gene_name"] <- pr_table$preferred_name[pr_table$protein_external_id == gene_id_A]
        #print(paste("gene_id_A",gene_inters))
        
        # Fetching the Pubmed articles that support each interaction and adds that in a new column "pubmed_ids"
        pubmed_arti_list <- as.list(str_db$get_pubmed_interaction(gene_id_A, gene_id_B))
        if(length(pubmed_arti_list!=0)){
            gene_inters[row, "pubmed_ids"] <- toString(list(pubmed_arti_list))
            #print(paste("pubmed_ids",gene_inters))
        } else{
            gene_inters[row, "pubmed_ids"] <- toString(list("NA"))
            #print(paste("pubmed_ids NA",gene_inters))
        } 
    }
    return(gene_inters)
}

create_network_n_info_tables <- function(st_db, df1, gl, outdir, prot_tab){ 
    # Mapping input genes to stringdb network
    genes_mapped <- st_db$map(df1, "gene", removeUnmappedRows = TRUE)
    hits <- genes_mapped$STRING_id    

    # Plotting and saving network
    www_subdir <- file.path(outdir, "www") 
    if (!dir.exists(www_subdir)){
        dir.create(www_subdir)
    }
    
    # Plotting and saving network
    pdf(paste(outdir, "stringdb_network.pdf", sep="/"))
    st_db$plot_network(hits)
    dev.off()

    # Save png of network to visualize in TIMEOR
    st_db$get_png(hits, file=paste(www_subdir,"stringdb_network.png",sep="/"))

    # Creating pubmed_ge=ne_lists (e.g. #pubmed_gene_lists <- c(string_db$mp("Cyt-c-p"), string_db$mp("Nle")))
    gene_list_ids <-c()
    for(i in gl){
        gene_list_ids <- append(gene_list_ids, st_db$mp(i))
    }

    # Fetching and saving network interactions and Pubmed information    
    gene_interactions <- st_db$get_interactions(gene_list_ids)

    # Fetching enrichment p-value
    pvalEnrich <- st_db$get_ppi_enrichment(gene_list_ids)$enrichment
    
    # Adding stringdb, protein id table, and enrichment p-value to table (in that order)
    info_table <- add_to_info_table(gene_interactions, st_db, prot_tab, pvalEnrich)
    
    return(info_table)
}

save_stringdb_network_n_table <- function(s_t, IO){
    s_t$pubmed_ids <-gsub('c','',s_t$pubmed_ids)
    s_t$pubmed_ids <-gsub('\\(','',s_t$pubmed_ids) 
    s_t$pubmed_ids <-gsub(')','',s_t$pubmed_ids)
    s_t$pubmed_ids <-gsub('list','',s_t$pubmed_ids)
    s_t$pubmed_ids <-gsub('"','',s_t$pubmed_ids)
    write.table(s_t, paste(IO, "stringdb_info_table.tsv", sep="/"), row.names=FALSE, sep="\t")
    print(paste("Saved STRINGdb network and info_tables to", IO, sep=" "))
}

main <- function(){  
    print("Stringdb")
    # Creating stringdb object by instantiating the STRINGdb reference class
    string_db <- STRINGdb$new(version="10", species=NCBI_TAXO, score_threshold=0, input_directory="") # default threshold is 400
    
    # Creating the table for protein IDs to map gene name to Pubmed ID
    protein_ID_table <- string_db$get_proteins()
    
    APP_DIR_redirect <- "/srv/"
    if(NCBI_TAXO == 7227){
        strdb_file_folder <- paste(APP_DIR_redirect,"/genomes_info/dme/",sep="")
        
        if(!file.exists(paste(strdb_file_folder,"7227__protein_aliases_tf.tsv.gz", sep=""))){   
            print("FILE DOES NOT EXIST")
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/protein_aliases/7227__protein_aliases_tf.tsv.gz', oD = strdb_file_folder)
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/protein_links/7227__protein_links.tsv.gz', oD = strdb_file_folder)
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/proteins/7227__proteins.tsv.gz', oD=strdb_file_folder)
        }
     } else if(NCBI_TAXO == 9606){
        strdb_file_folder <- paste(APP_DIR_redirect,"/genomes_info/hsa/",sep="")
        if(!file.exists(paste(strdb_file_folder,"9606__protein_aliases_tf.tsv.gz", sep=""))){
            print("FILE DOES NOT EXIST")
            print(paste(strdb_file_folder,"9606__protein_aliases_tf.tsv.gz"))
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/protein_aliases/9606__protein_aliases_tf.tsv.gz', oD = strdb_file_folder)
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/protein_links/9606__protein_links.tsv.gz', oD = strdb_file_folder)
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/proteins/9606__proteins.tsv.gz', oD=strdb_file_folder)
         }
     } else if(NCBI_TAXO == 10090){
        strdb_file_folder <- paste(APP_DIR_redirect,"/genomes_info/mmu/",sep="")
        if(!file.exists(paste(strdb_file_folder,"10090__protein_aliases_tf.tsv.gz", sep=""))){
            print("FILE DOES NOT EXIST")
            print(paste(strdb_file_folder,"10090__protein_aliases_tf.tsv.gz"))
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/protein_aliases/10090__protein_aliases_tf.tsv.gz', oD = strdb_file_folder)
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/protein_links/10090__protein_links.tsv.gz', oD = strdb_file_folder)
            downloadAbsentFile('http://string.uzh.ch/permanent/string/10/proteins/10090__proteins.tsv.gz', oD=strdb_file_folder)
         }
     } else{
        stop("Add necessary organism information.")
     }
    
    # Getting results folder name for user's experiment 
    exp <- basename(IN_OUTPUT)
    print(paste("Results folder name: ",exp, sep=" "))
    
    # Iterating through all experiments in list_exp
    if(!grepl("_results", exp, fixed=T)){
        stop("Please pass in a named results folder (e.g. test_results)")
    }else{
       
        if(TFs_NOT_CLUSTERS){ # Iterating through list of TFs across clusters
            # Getting tf_reg_network folder
            dir_tf_rn <- paste(IN_OUTPUT, "temporal_relations", sep="/")

            # Collecting input data
            df <- read.table(paste(dir_tf_rn,"string_db_input.csv",sep="/"), sep=",", header=FALSE, skip = 1, check.names=FALSE)[,1:2]
            colnames(df) = c("gene","change")
            print(paste("Processing", dir_tf_rn, sep=" "))
            print(head(df))

            # Creating a list of genes within each cluster    
            a <- dput(as.character(df[1][,1]))
            
            # Removing any ballgown ID and/or transcript ID if present at end of gene name
            if (grepl("_", a[1])){
                cat("Removing values after and including '_' at end of gene/transcript name\n")
                gene_list <- gsub("\\-R[A-Z][_|\\>].*","",a)
            } else if (grepl("-R[A-Z]\\>", a[1])) {
                cat("Removing transcript ID at end of gene name\n")
                gene_list <- gsub("\\-R[A-Z]\\>.*","",a)
            } else{
                gene_list <- a
            }

            # Creating STRINGdb network and info. tables then save
            string_table <- create_network_n_info_tables(string_db, df, gene_list, dir_tf_rn, protein_ID_table)
            save_stringdb_network_n_table(string_table, dir_tf_rn)
            print(paste("Saved STRINGdb network and info_tables to", dir_tf_rn, sep=" "))
            
        } else{ # Iterating through clusters

            # Getting cluster directories
            dir_clusts <- paste(IN_OUTPUT, "clusters", sep="/")
            
            # Finding ONLY filenames that are *geneList*.csv
            filenames <- Sys.glob(file.path(dir_clusts, "/*/", "*geneList*.csv"))
            
            # For each gene cluster within each condition
            for(dir_file in filenames){
                OUTPUT_DIR <-paste(dirname(dir_file),'',sep="/")

                # Collecting input data
                df <- read.table(dir_file, sep=",", header=FALSE, skip = 1, check.names=FALSE)[,1:2]
                colnames(df) = c("gene","change")
                print(paste("Processing", dir_file, sep=" "))
                print(head(df))

                # Creating a list of genes within each cluster    
                a <- dput(as.character(df[1][,1]))
                
                # Removing any ballgown ID and/or transcript ID if present at end of gene name
                if (grepl("_", a[1])){
                    cat("Removing values after and including '_' at end of gene/transcript name\n")
                    gene_list <- gsub("\\-R[A-Z][_|\\>].*","",a)
                } else if (grepl("-R[A-Z]\\>", a[1])) {
                    cat("Removing transcript ID at end of gene name\n")
                    gene_list <- gsub("\\-R[A-Z]\\>.*","",a)
                } else{
                    gene_list <- a
                }

                # Creating STRINGdb network and info. tables then save
                string_table <- create_network_n_info_tables(string_db, df, gene_list, OUTPUT_DIR, protein_ID_table)
                save_stringdb_network_n_table(string_table, OUTPUT_DIR)
                print("save_stringdb_network_n_table")
            }
        }
    }
}

if(!interactive()) {
    main()
}