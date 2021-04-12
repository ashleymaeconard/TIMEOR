# stringdb.r
# Ashley Mae Conard
# Last Mod. 1/9/2020
# Finds the network (if any) connecting genes, given a list of genes.

# Checking input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type: /usr/local/bin/Rscript stringdb.r 
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/a_results/clusters/1/) 
            2) ORGANISM_NCBI_ID
            3) APP_DIR") # NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
             
} else if (length(args) == 3) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/a_results/clusters/1/) 
            2) ORGANISM_NCBI_ID
            3) APP_DIR") # NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
}

# Inputting Arguments
IN_OUTPUT <- args[1]
NCBI_TAXO <- strtoi(args[2]) # NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
APP_LOC <- args[3] 
timeor_app_folder <- dirname(APP_LOC)

# Installing libraries
library(STRINGdb)
library(dplyr)
library(stringr)

add_to_info_table <- function(gene_inters, str_db, pr_table, pval){
    # Creating new columns "A_gene_name" and "B_gene_name"
    gene_inters[,"A_gene_name"] <- NA
    gene_inters[,"B_gene_name"] <- NA
    
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
    }
    return(gene_inters)
}

create_network_n_info_tables <- function(st_db, df1, gl, outdir, prot_tab, det_links){ 
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
    
    # Fetching and saving experimental network interactions   
    interacting_genes <- st_db$get_interactions(gene_list_ids)
    uniq_gene_interactions <- data.frame(unique(interacting_genes))
    gene_interactions_n_types <- data.frame(from=character(),
                 to=character(), 
                 neighborhood=integer(),
                 fusion=integer(),
                 cooccurrence=integer(),
                 coexpression=integer(),
                 experimental=integer(),
                 database=integer(),
                 textmining=integer(),
                 combined_score=integer(),
                 stringsAsFactors=FALSE) 

    for(u in rownames(uniq_gene_interactions)){
        gene_A <- uniq_gene_interactions[u, 1]
        gene_B <- uniq_gene_interactions[u, 2]

        row_to_add = det_links %>% filter(protein1==gene_A & protein2==gene_B & experimental > 0)
        print("row_to_add")
        print(row_to_add)
        if(nrow(row_to_add)!=0){
            gene_interactions_n_types[nrow(gene_interactions_n_types) + 1,] = 
                det_links %>% filter(protein1==gene_A & protein2==gene_B & experimental > 0)
            gene_interactions_n_types[nrow(gene_interactions_n_types),]$from <- gene_A
            gene_interactions_n_types[nrow(gene_interactions_n_types),]$to <- gene_B
        }
    }
    # Fetching enrichment p-value
    pvalEnrich <- st_db$get_ppi_enrichment(gene_list_ids)$enrichment
    
    # Adding stringdb, protein id table, and enrichment p-value to table (in that order)
    info_table <- add_to_info_table(gene_interactions_n_types, st_db, prot_tab, pvalEnrich)
    write(" ", stderr())
    write("PVALLLLL ENRICH", stderr())
    write(" ", stderr())
    write(pvalEnrich, stderr())
    
    return(info_table)
}

save_stringdb_network_n_table <- function(s_t, IO){
    write.table(s_t, paste(IO, "stringdb_info_table.tsv", sep="/"), row.names=FALSE, sep="\t")
    print(paste("Saved STRINGdb network and info_tables to", IO, sep=" "))
}


main <- function(){
    
    APP_DIR_redirect <- "/srv/"
    # Getting alias files if needed
    if(NCBI_TAXO == 7227){
        strdb_file_folder <- paste(APP_DIR_redirect,"/genomes_info/dme/",sep="")
        detailed_links <- read.table(paste(strdb_file_folder, "7227.protein.links.detailed.v11.0.txt", sep="/"), header = TRUE, sep = " ")

        if(!file.exists(paste(strdb_file_folder,"7227.protein.aliases.v11.0.txt.gz", sep="/"))){   
            cat("FILE DOES NOT EXIST", file=stderr())
            downloadAbsentFile('https://stringdb-static.org/download/protein.aliases.v11.0/7227.protein.aliases.v11.0.txt.gz', oD = strdb_file_folder)
            downloadAbsentFile('https://stringdb-static.org/download/protein.links.v11.0/7227.protein.links.v11.0.txt.gz', oD = strdb_file_folder)
            downloadAbsentFile('https://stringdb-static.org/download/protein.info.v11.0/7227.protein.info.v11.0.txt.gz', oD=strdb_file_folder)
        }
     } else if(NCBI_TAXO == 9606){
        strdb_file_folder <- paste(APP_DIR_redirect,"/genomes_info/hsa/",sep="")
        detailed_links <- read.table(paste(strdb_file_folder, "9606.protein.links.detailed.v11.0.txt", sep="/"), header = TRUE, sep = " ")

        if(!file.exists(paste(strdb_file_folder,"9606.protein.aliases.v11.0.txt.gz", sep="/"))){
            cat("FILE DOES NOT EXIST", file=stderr())
            downloadAbsentFile('https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz', oD = strdb_file_folder)
            downloadAbsentFile('https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz', oD = strdb_file_folder)
            downloadAbsentFile('https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz', oD=strdb_file_folder)
         }
     } else if(NCBI_TAXO == 10090){
        strdb_file_folder <- paste(APP_DIR_redirect,"/genomes_info/mmu/",sep="")
        detailed_links <- read.table(paste(strdb_file_folder, "10090.protein.links.detailed.v11.0.txt", sep="/"), header = TRUE, sep = " ")

        if(!file.exists(paste(strdb_file_folder,"10090.protein.aliases.v11.0.txt.gz", sep="/"))){
            cat("FILE DOES NOT EXIST", file=stderr())
            downloadAbsentFile('https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz', oD = strdb_file_folder)
            downloadAbsentFile('https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz', oD = strdb_file_folder)
            downloadAbsentFile('https://stringdb-static.org/download/protein.info.v11.0/10090.protein.info.v11.0.txt.gz', oD=strdb_file_folder)
         }
     } else{
        stop("Add necessary organism information.")
     }
    
    # Creating stringdb object by instantiating the STRINGdb reference class
    string_db <- STRINGdb$new(version="11", species=NCBI_TAXO, score_threshold=0, input_directory=strdb_file_folder) # the default threshold is 400
    print("stringdb")
    
    # Creating the table for protein IDs to map gene name to Pubmed ID
    protein_ID_table <- string_db$get_proteins()
    print("protein_ID_table")

    # Process genes in a given cluster 
    # Finding ONLY filenames that are *geneList*.csv
    filename <- Sys.glob(file.path(IN_OUTPUT, "/*geneList*.csv"))
    print("filename")

    # Collecting input data
    df <- read.table(filename, sep=",", header=FALSE, skip = 1, check.names=FALSE)[,1:2]
    colnames(df) = c("gene","change")
    print(paste("Processing", filename, sep=" "))
    print("collect input data")

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
    print("string_db")
    print(string_db)
    print("df")
    print(df) 
    print("gene_list")
    print(gene_list)
    print("IN_OUT")
    print(IN_OUTPUT)
    print("protein_ID_table")
    print(protein_ID_table)

    # Creating STRINGdb network and info. tables then save
    string_table <- create_network_n_info_tables(string_db, df, gene_list, IN_OUTPUT, protein_ID_table, detailed_links)
    save_stringdb_network_n_table(string_table, IN_OUTPUT)
    print("save_stringdb_network_n_table")
}

if(!interactive()) {
    main()
}