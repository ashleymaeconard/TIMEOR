# get_top_tfs.r
# Ashley Mae Conard
# Last Mod: 2/13/2020
# Purpose: Runs rcis_target in 3 steps* then TIMEOR ranks results to find TOP TFs sets, creating 5 table types**.
#   * Rcistarget steps are (1) motif enrichment analysis, (2) motif-TF annotation, and (3) selection of significant genes (taken from Rcistarget tutorial)
#   ** Table 1: Rcistarget interactive output table via Shiny
#      Table 2: Rcistarget output table
#      Table 3: ENCODExplorer outputs table of ENCODE accession numbers and other information for ChIP-seq data for TF where applicable
#      Table 4: Top TFs highlighted in each database/method via Rcistarget, and the mode, percent_concordance, and encode numbers where applicable
#      Table 5: All clusters top TFs and accession numbers (where applicable)
# References: https://yulab-smu.github.io/clusterProfiler-book/chapter12.html          

# CHECKING ARGUMENTS 
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type: /usr/local/bin/Rscript get_top_tfs.r 
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/analysis/NAME_results/) 
            2) ORGANISM (dme, hsa, mmu) 
            3) MIN_NORM_ENRICH_SCORE (recommend 3)
            4) TOP_X (recommend returning the top 4 TFs per database)
            5) PERCENT_CONCORDANCE (between methods not returning NA - recommend 40%)
            6) APP_DIR (e.g. ~/Desktop/TIMEOR/app/)
	    7) TF_CONFIDENCE (hi, low)", call.=FALSE)

} else if (length(args) == 7) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/) 
                  2) ORGANISM (dme (Drosophila melanogaster), hsa (Homo sapiens), mmu (Mus musculus)
                  3) MIN_NORM_ENRICH_SCORE (recommend 3)
                  4) TOP_X (recommend returning the top 4 TFs per database)
                  5) PERCENT_CONCORDANCE (between methods not returning NA - recommend 40%)
                  6) APP_DIR (e.g. ~/Desktop/TIMEOR/app/)
    		  7) TF_CONFIDENCE (hi, low)")
}

# ASSIGNING INPUT ARGUMENTS
IN_OUTPUT <- args[1] #e.g. ~/timeor/results/primary/
ORGANISM <- args[2] # dme (Drosophila melanogaster), hsa (Homo sapiens), mmu (Mus musculus)
options(digits=5)
MIN_NORM_ENRICH_SCORE <- as.integer(args[3]) # recommend 3
TOP_X <- as.integer(args[4]) # recommend 4
PERCENT_CONCORDANCE <- as.double(args[5]) # recommend 40
APP_DIR <- args[6]
CONF <-args[7]

# LOADING NECESSARY PACKAGES
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
library(RcisTarget)

# ASSIGNING ORGANISM LIBRARY AND ANNOTATIONS
if(ORGANISM=="dme"){
    data("motifAnnotations_dmel")
    ENCODE_ORG = "Drosophila melanogaster"
    motifRankings <- importRankings(paste(APP_DIR,"/../genomes_info/dme/dm6-5kb-upstream-full-tx-11species.mc8nr.feather",sep="/"))
    TFS_FILE <- "/../genomes_info/dme/tfs_all.csv"

} else if(ORGANISM=="hsa"){
    data("motifAnnotations_hgnc")
    ENCODE_ORG = "Homo sapiens"
    motifRankings <- importRankings(paste(APP_DIR,"/../genomes_info/hsa/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",sep="/"))
    TFS_FILE <- "/../genomes_info/hsa/tfs_all.csv"

} else if(ORGANISM=="mmu"){
    data("motifAnnotations_mgi")
    ENCODE_ORG = "Mus musculus"
    motifRankings <- importRankings(paste(APP_DIR,"/../genomes_info/mmu/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",sep="/"))
    TFS_FILE <- "/../genomes_info/mmu/tfs_all.csv"

} else{
    stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
}

# LOADING REST OF PACKAGES
library(doMC)
library(ENCODExplorer)
library(doRNG)
library(doMC)
library(data.table)
library(DT)
library(stringr)
library(dplyr)
library(purrr)

find_perturbed_tfs <- function(geneLi, all_clust_top_X_tfs, geneClustNum, tfs_tfas, perturbed_tfs_list){
    
    # Filling all_clust_top_X_tfs dataframe perturbed_TFs column for geneClustNum cluster with a list of all perturbed TFs 
    for(i in 1:length(geneLi$cluster)){
        if(any(tfs_tfas$HGNC_symbol==geneLi$cluster[i])){ # check if perturbed gene is a TF 
            perturbed_tfs_list <- union(perturbed_tfs_list, geneLi$cluster[i])
        }
    }
    all_clust_top_X_tfs$perturbed_TFs[geneClustNum] <- list(perturbed_tfs_list)
    all_clust_top_X_tfs$cluster <- seq.int(nrow(all_clust_top_X_tfs))
    return(all_clust_top_X_tfs)
}

report_n_collect_encode_results <- function(top_X_tf_per_db_cleaned, OUTDIR){
    # Gathering ENCODE numbers if exist
    top_X_tf_per_db_cleaned$encode_accession_num <- "NA"
    for(gene in 1:nrow(top_X_tf_per_db_cleaned)){
        query_res <- queryEncodeGeneric(organism=ENCODE_ORG, assay="ChIP-seq", file_format="bigWig", target=top_X_tf_per_db_cleaned$mode_TF[gene], fuzzy=T)
        top_X_tf_per_db_cleaned$encode_accession_num[gene] <- list(unique(query_res$accession))
        per_conc <- as.double(format(top_X_tf_per_db_cleaned$prct_concordance[gene], digits=3))
        
        # Saving Table 3: ENCODExplorer outputs table of ENCODE accession numbers and other information for ChIP-seq data for TF where applicable
        print(paste0("Saving Table 3 in: ",OUTDIR))
        print("ENCODExplorer outputs ENCODE accession numbers table and other information for ChIP-seq data for TF where applicable.")
        write.csv(query_res, paste(OUTDIR, paste("top_tfs/encode_TF",gene,"chip",top_X_tf_per_db_cleaned$mode_TF[gene],per_conc,"prctConcor.csv", sep="_"), sep="/"), row.names=FALSE)
    }
    return(top_X_tf_per_db_cleaned)
}

top_X_tf_per_db <- function(rcistarget_res){
    # Adding database column to rcistarget results table
    rcistarget_res$database <- word(rcistarget_res$motif, 1, sep="_") # splitting motif str on "_"

    # Creating data.frame of top 4 TFs per database (Table 2: RcisTarget results output table)
    dbs <- unique(rcistarget_res$database)
    top_X_tf_each_db <- data.frame(matrix(ncol = length(dbs), nrow = TOP_X))
    colnames(top_X_tf_each_db) <- dbs 

    # Delimeters to remove in TF_highConf column
    delimeters_to_remove <- c(";",".","-"," ","(",")")

    if(CONF=="low"){
    	# Adding top 4 TFs to dataframe from TF low_confidence
    	for(db in unique(rcistarget_res$database)){
        	top_X_tf_each_db[db] <- word(rcistarget_res[which(rcistarget_res$database %like% db)][1:TOP_X]$TF_lowConf, 1, sep=" ")
    	}

    # Adding top 4 TFs to dataframe per database
    }else if(CONF=="high"){
	for(db in unique(rcistarget_res$database)){
        	top_X_tf_each_db[db] <- word(rcistarget_res[which(rcistarget_res$database %like% db)][1:TOP_X]$TF_highConf, 1, sep=" ")
    	}
    } else{
    	stop("Enter 'hi' or 'low' to argument 7 to indicate TF confidence type.")
    }

    # Adding NA if empty cell in dataframe
    top_X_tf_each_db_cleaned <- as.data.frame(apply(top_X_tf_each_db, 2, function(x)gsub(delimeters_to_remove, '', x))) %>% mutate_all(na_if,"")
    
    # Converting factor dataframe to character
    top_X_tf_each_db_cleaned %>% map_if(is.factor, as.character) %>% as_data_frame -> top_X_tf_each_db_cleaned

    # Getting mode of each row
    t_top_X_tf_each_db_cleaned<- as.data.frame(t(top_X_tf_each_db_cleaned))
    #mode_top_X_tf_each_db <- calculate_mode(t_top_X_tf_each_db_cleaned$V1)
    mode_top_X_tf_each_db <- apply(top_X_tf_each_db_cleaned, 1, function(x) { tab <- table(x); names(tab)[which.max(tab)] } )
    #apply(df[ ,2:length(df)], 1, mfv)

    # Calculating percent concordance and only reporting top TFs if concensus more than 40% for those methods with output (i.e. do not return NA)
    divide = get("/")
    times = get("*")
    per_conc_list = c()
    
    for(row in 1:nrow(top_X_tf_each_db_cleaned)){
        # Find number of occurrences of the mode TF
        num_mode_tf_occur <- as.double(length(grep(paste("^",mode_top_X_tf_each_db[row],"$", sep=""), top_X_tf_each_db_cleaned[row,])))

        # Get total number of methods reporting a TF
        row_df <- unlist(asplit(top_X_tf_each_db_cleaned,1)[row],recursive = TRUE, use.names = TRUE)
        denom_num_method_report_tf <- as.double(length(row_df[!is.na(row_df)]))
        
        # Calculate percent concordance for top X TFs
        perc_concord <- times(divide(num_mode_tf_occur,denom_num_method_report_tf),100)
        per_conc_list <- c(per_conc_list, perc_concord)
    }
    top_X_tf_each_db_cleaned$mode_TF <- mode_top_X_tf_each_db
    top_X_tf_each_db_cleaned$prct_concordance <- per_conc_list
    
    return(top_X_tf_each_db_cleaned)
}

run_rcistarget <- function(geneL, motifR, OUTDIR){
    
    # Calculating motif enrichment
    if(ORGANISM=="dme"){
        motifEnrichmentTable_wGenes <- cisTarget(geneL, motifR, motifAnnot=motifAnnotations_dmel)

    } else if(ORGANISM=="hsa"){
        motifEnrichmentTable_wGenes <- cisTarget(geneL, motifR, motifAnnot=motifAnnotations_hgnc)

    }else if(ORGANISM=="mmu"){
        motifEnrichmentTable_wGenes <- cisTarget(geneL, motifR, motifAnnot=motifAnnotations_mgi)

    } else{
        stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
    }
    motifs_AUC <- calcAUC(geneL, motifR, nCores=1)

    # # AUC histogram of motif scores
    # auc <- getAUC(motifs_AUC)["RELA_targets",]
    # hist(auc, main="RELA_targets", xlab="AUC histogram",
    #     breaks=100, col="#ff000050", border="darkred")
    # nes3 <- (3*sd(auc)) + mean(auc)
    # abline(v=nes3, col="red")

    # Selecting significant motis and annotate to TFs
    if(ORGANISM=="dme"){
        motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=MIN_NORM_ENRICH_SCORE, motifAnnot=motifAnnotations_dmel)

    } else if(ORGANISM=="hsa"){
        motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=MIN_NORM_ENRICH_SCORE, motifAnnot=motifAnnotations_hgnc)

    }else if(ORGANISM=="mmu"){
        motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=MIN_NORM_ENRICH_SCORE, motifAnnot=motifAnnotations_mgi)

    } else{
        stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
    }

    motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
        rankings=motifR,
        method="iCisTarget", 
        geneSets=geneL)
    dim(motifEnrichmentTable_wGenes)    
    motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
    resultsSubset <- motifEnrichmentTable_wGenes_wLogo

    # Creating datatable
    DT_resultsSubset <- datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE],
        escape = FALSE, # To show the logo
        filter="top", options=list(pageLength=20))
    datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], escape = FALSE,filter="top")
    # Save Shiny datatable and create output directory if does not exist
    print("Saving Table 1: Rcistarget interactive output table via Shiny.")
    dir.create(file.path(OUTDIR, "top_tfs"), showWarnings = FALSE)
    saveWidget(DT_resultsSubset, paste(OUTDIR,'top_tfs/interactive_rcistarget_results.html',sep="/"))
    return(resultsSubset)
}

main <- function(){
    
    # Getting all transcription factors for input organism
    if(ORGANISM=="dme"){
        is_tf  <- fread(paste(APP_DIR, TFS_FILE, sep="/"), sep=",", select=c("HGNC_symbol"))
    } else if(ORGANISM=="hsa"){
        all_tfs  <- fread(paste(APP_DIR, TFS_FILE, sep="/"), sep=",", select=c("HGNC_symbol","Is_TF"))
        is_tf <- all_tfs[all_tfs$Is_TF!="No"] #"Unlikely to be sequence specific TF"
    } else{
        stop("Need to download mmu (Mus musculus) TF table")
    }

    # Getting results folder name for user's experiment 
    exp <- basename(IN_OUTPUT)
    print(paste("Results folder name: ",exp, sep=" "))
    
    # Iterating through all experiments in list_exp
    if(!grepl("_results", exp, fixed=T)){
        stop("Please pass in a named results folder (e.g. test_results)")
    }else{
        # Creating perturbed TF list to fill across all clusters
        pert_g_l <- c()

        # Creating factor_binding folder
        dir_tf_rn <- paste(IN_OUTPUT, "factor_binding", sep="/")
        dir.create(dir_tf_rn, showWarnings=FALSE)
        # Finding the cluster directories for each experiment (in list_exp)
        dir_clusts <- paste(IN_OUTPUT, "clusters", sep="/")
        # Finding only filenames that are *geneList*.csv or *geneName*.csv to get list of genes to calc. enrichment
        filenames <- Sys.glob(file.path(dir_clusts, "/*/", "*geneList*.csv"))
        # Create Table 5: All clusters top perturbed and putative TFs with encode accession numbers (where applicable)
        x <- 1:TOP_X
        col_names <- c("perturbed_TFs")
        col_names <- union(col_names, paste("top_TF",x, sep="_"))
        col_names <- union(col_names, paste("ENCODE_accession_num_TF",x,sep="_"))
        all_clusters_top_X_tfs <- data.frame(lapply(data.frame(matrix(ncol = (TOP_X*2)+1, nrow = length(filenames))), as.character),stringsAsFactors=FALSE) # length(filenames) is number of clusters
        colnames(all_clusters_top_X_tfs) <- col_names
        print(filenames)
        # For all gene clusters within each experiment
        for(dir_file in filenames){
            geneNumClust=(sub(".*_([^.]+)\\.csv.*", "\\1", dir_file))
            OUTPUT_DIR <-paste(dirname(dir_file),'',sep="/")
            # Creating gene lists
            geneLists <- list(cluster=read.table(dir_file, sep=",",stringsAsFactors=FALSE, skip=1)[,1])
            # Running Rcistarget
            rcistarget_results_table <- run_rcistarget(geneLists, motifRankings, OUTPUT_DIR)
            # Saving Table 2: RcisTarget results output table
            print("Saving Table 2: RcisTarget results output table")
            write.csv(rcistarget_results_table, paste(OUTPUT_DIR, "top_tfs/rcistarget_results.csv", sep="/"), row.names=FALSE)
            # Creating and saving Table 2 with RcisTarget output, the mode, percent_concordance, and encode numbers where applicable
            final_TF_list <- top_X_tf_per_db(rcistarget_results_table)
            
            # Adding ENCODE numbers where applicable
            print("Adding ENCODE numbers where applicable to Table 2: RcisTarget results output table.")
            top_X_n_encode <-report_n_collect_encode_results(as.data.frame(final_TF_list), OUTPUT_DIR)

            # Saving Table 4: Top TFs highlighted in each database/method via Rcistarget, and the mode, percent_concordance, and encode numbers where applicable
            print("Saving Table 4: Top TFs highlighted in each database/method via Rcistarget, and the mode, percent_concordance, and encode numbers where applicable.")
            is.na(top_X_n_encode) <- top_X_n_encode == "motif"
            fwrite(top_X_n_encode, file = paste(OUTPUT_DIR, "top_tfs/top_TFs_per_method_n_encode_nums.csv", sep="/"))

            # Filling Table 5: All clusters top perturbed and putative TFs with encode accession numbers (where applicable)
            for(i in 1:nrow(top_X_n_encode)){
                # Only add those top TFs that pass the percent concordance threshold set by the user
                if(as.double(top_X_n_encode$prct_concordance[i]) >= as.double(PERCENT_CONCORDANCE)){
                    # Add top X TF
                    col_tf <- paste("top_TF",i,sep="_")
                    all_clusters_top_X_tfs[[col_tf]][as.integer(geneNumClust)] <- top_X_n_encode$mode_TF[i]
                    # Add top X TF ENCODE IDs
                    col_enc <- paste("ENCODE_accession_num_TF",i,sep="_")
                    all_clusters_top_X_tfs[[col_enc]][as.integer(geneNumClust)] <- top_X_n_encode$encode_accession_num[i]
                }
            }
           
            # Adding perturbed genes to Table 5: All clusters top perturbed and putative TFs with encode accession numbers (where applicable)
            all_clusters_top_X_tfs <- find_perturbed_tfs(geneLists, all_clusters_top_X_tfs, as.integer(geneNumClust), is_tf, pert_g_l)
            print(all_clusters_top_X_tfs)
        }

        # Saving Table 5: All clusters top perturbed and putative TFs with encode accession numbers (where applicable)
        print("Saving Table 5: All clusters top perturbed and putative TFs with encode accession numbers (where applicable)")
        is.na(all_clusters_top_X_tfs) <- all_clusters_top_X_tfs == "NULL"
        fwrite(all_clusters_top_X_tfs, file = paste(dir_tf_rn, "perturbed_putative_tfs_n_encode.csv", sep="/"))
    }
}

if(!interactive()) {
    print("Creating 5 table types.")
    main()
}
