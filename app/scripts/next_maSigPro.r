# next_MaSigPro.r
# Ashley Mae Conard
# Last Mod. 6/13/2020
# Purpose: Run nextMaSigPro to identify differentially expressed genes from time series data.

# Input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    write("Type:next_MaSigPro.r /FULL/PATH/TO/METDATA_FILE/ /FULL/PATH/TO/COUNT_MATRIX_FILE/ (normalized and corrected) /FULL/PATH/TO/OUTPUTDIR/ RESULT_FOLDER_NAME (e.g. insulin_stim) PVAL_THRESH ORGANISM", call.=FALSE, stderr())
} else if (length(args) == 6) {
    cat("Passed in:", args,"\n")
} else{
    write("Pass in 6 arguments. Type:next_MaSigPro.r /FULL/PATH/TO/METDATA_FILE/ /FULL/PATH/TO/COUNT_MATRIX_FILE/ (normalized and corrected) /FULL/PATH/TO/OUTPUTDIR/ RESULT_FOLDER_NAME (e.g. insulin_stim) PVAL_THRESH ORGANISM", stderr())
}

METADATA <- args[1] 
COUNTDATA <- args[2]
OUTDIR <- args[3]
RES_FOLDER <- args[4]
PVAL_THRESH <- as.numeric(args[5])
ORGANISM <- args[6]

# Import libraries
library(maSigPro)
library(MASS)
library(dplyr)
library(tidyverse)
library(tibble)
library(data.table)

# Assigning organism library
if(ORGANISM=="dme"){
    ORG_DB="org.Dm.eg.db"
    source("./scripts/geneID_converter.r", local=TRUE)
} else if(ORGANISM=="hse"){
    ORGANISM="hsa"
    require(biomaRt)
    ORG_DB="org.Hs.eg.db"
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
}else if(ORGANISM=="mus"){
    ORGANISM="mmu"
    require(biomaRt)
    ORG_DB="org.Mm.eg.db"
    mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
} else{
    write("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)", stderr())
}
library(ORG_DB, character.only = TRUE) # organism database library

# Import metaData
metaData <- read.csv(file=METADATA,row.names=1)
setnames(metaData, old=c("condition","time", "batch"), new=c("Group", "Time", "Batch"))

# Convert Group (i.e. condition) column to numeric for processing
iter <- 1
tmp <- as.character(metaData[,1][[1]]) # first row's group value
intGroup <- c()
print("first")
print(tmp)
for(i in metaData$Group){
  if(tmp!=i){
    iter <- iter + 1
    tmp <- i
  }
  intGroup <- c(intGroup,iter)
}
print("intGroup")
metaData$Group <- intGroup
print(metaData)
    
# Import 
ct = read.csv(file=COUNTDATA, sep=",",row.names=1)
ct<-na.omit(ct)
#write(head(ct), stderr())

# Defining the regression model
design <- make.design.matrix(metaData, degree = 1)

# Finding significant genes
write("design", stderr())
fit <- p.vector(ct, metaData)#,counts=TRUE) # "By default maSigPro corrects this p-value for multiple comparisons by applying the linear step-up (B-H) false discovery rate (FDR) procedure (Benjamini and Hochberg, 1995)."
write("fit", stderr())
# For edification:
  # fit$i # "returns the number of significant genes"
  # fit$alfa # "gives p-value at the Q false discovery control level"
  # fit$SELEC # "is a matrix with significant genes and their expression values"

# Finding significant differences
tstep <- T.fit(fit, step.method = "backward", alfa = as.numeric(PVAL_THRESH))
write("tstep", stderr())

# Obtaining lists of significant genes
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "each")

# Get number of DEGs
gene_names <- as.list(sigs$summary['independ'])
write("Number of differentially expressed genes: ", stderr())
write(length(gene_names[[1]]), stderr())

# Creating temporary dataframe
fit_select <- as.data.frame(fit$SELEC)
de_genes_nmsp <- subset(fit_select, rownames(fit_select) %in% sigs$summary['Time'][[1]])

# Initializing expression dataframe
de_expr_df <- data.frame(rownames(de_genes_nmsp))
colnames(de_expr_df) <- "ID"
rownames(de_expr_df) <- de_expr_df$"ID"
de_expr_df = de_expr_df[!(names(de_expr_df) %in% "ID")]

# Extract sample list (without replicate name)
samples.rep_list <- as.list(colnames(de_genes_nmsp))
l = list()
i <- 1
while(i<= length(as.list(colnames(de_genes_nmsp)))) { # (s in as.list(colnames(tmp_df))){
    #print(gsub("\\..*","",samples.rep_list[i]))
    l[[i]] <- gsub("\\..*","",samples.rep_list[i])
    i <- i + 1
}
samp_names <- unique(l)

# Loop through all columns and append means to new dataframe
for(s in samp_names){

    # Get replicate columns
    rep_cols <- de_genes_nmsp[, grepl(s, colnames(de_genes_nmsp))]

    # Get mean of replicate columns
    mean_reps <- subset(data.frame(ID=rep_cols, Means=rowMeans(rep_cols)),select=c("Means"))
    colnames(mean_reps) <- s

    # Merge dataframes
    merged_dfs <- merge(de_expr_df, mean_reps, by=0)
    rownames(merged_dfs) <- merged_dfs$"Row.names"
    merged_dfs= merged_dfs[!(names(merged_dfs) %in% "Row.names")]

    # Reset for next iteration (next sample)
    de_expr_df <- merged_dfs
    merged_dfs <- data.frame()
}

# Calculate Z-score for each: (x - mean(x)) / sd(x)
tmp_clustmap_zs <-as.data.frame(t(apply(de_expr_df, 1, scale)))
colnames(tmp_clustmap_zs) <- samp_names
#write(tmp_clustmap_zs, stderr())

# Get gene IDS to convert
gene_id_list <- sub('\\.[0-9]*$', '', rownames(tmp_clustmap_zs))

if(ORGANISM=="dme"){
  # Set database to work with
  gene_ID_database <- toTable(org.Dm.egFLYBASE)
  gene_ID_database_name <- "flybase"

  # Add gene name to resulting dataframe
  ids.type  <- gene_ID_database_name
  ids <- gene_id_list
  res <- as.vector(get.symbolIDsDm(ids,ids.type))
  tmp_clustmap_zs['gene_name'] <- res
  tmp_clustmap_zs$gene_id <- rownames(tmp_clustmap_zs)
  clustermap_zscore <- tmp_clustmap_zs %>% dplyr::select(gene_name, gene_id, everything())

} else if(ORGANISM=="hsa" || ORGANISM=="mmu"){
  bm_ids <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = gene_id_list, mart = mart)
  gene_tmp_clustmap_zs <- merge(tmp_clustmap_zs, bm_ids, by.x=0, by.y="ensembl_gene_id")
  gene_id_name_clustmap_zs <- gene_tmp_clustmap_zs %>% plyr::rename(c(Row.names = "gene_id", external_gene_name = "gene_name"))
  clustermap_zscore <- gene_id_name_clustmap_zs %>% dplyr::select(gene_name, gene_id, everything())
 
} else{
    write("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)", stderr())
}

# Check if NA in gene_name, copy gene_id in its place
if (NA %in% clustermap_zscore$gene_name){
  clustermap_zscore$gene_name <- ifelse(is.na(clustermap_zscore$gene_name), clustermap_zscore$gene_id, clustermap_zscore$gene_name)
}

# Save df
FULL_RES_DIR <- paste(RES_FOLDER,"results",sep="_")
dir.create(file.path(OUTDIR,FULL_RES_DIR), showWarnings = FALSE, recursive = TRUE)
FULL_OUTDIR <- paste(OUTDIR,FULL_RES_DIR,"nextmasigpro",sep="/")
if (!dir.exists(FULL_OUTDIR)){
  dir.create(FULL_OUTDIR, recursive = TRUE)
} else {
  write("NextMaSigPro directory already exists.", stderr())
}
outputLoc <- paste("NextMaSigPro results are here: ", FULL_OUTDIR)
write(outputLoc, stderr())
write.csv(clustermap_zscore,paste(FULL_OUTDIR,paste('nextMaSigPro_clustermapInput_padj',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE)