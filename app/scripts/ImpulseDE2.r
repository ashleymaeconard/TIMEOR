# Ashley Conard
# ImpulseDE2.r
# Last Mod. 12/12/2019
# Purpose: Run ImpulseDE2 to identify differentially expressed genes from time series data.

# Checking input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type:ImpulseDE2.r /FULL/PATH/TO/METDATA_FILE/ /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) /FULL/PATH/TO/OUTPUTDIR/ RES_FOLDER (e.g. insulin_stim) PVAL_THRESH", call.=FALSE)
} else if (length(args) == 5) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 5 arguments. Type:ImpulseDE2.r /FULL/PATH/TO/METDATA_FILE/ /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) /FULL/PATH/TO/OUTPUTDIR/ RES_FOLDER (e.g. insulin_stim) PVAL_THRESH")
}

# Input arguments
METADATA <- args[1] 
COUNTDATA <- args[2] 
OUTDIR <- args[3] 
RES_FOLDER <- args[4]
PVAL_THRESH <- as.numeric(args[5])

# Set up output directories
dir.create(file.path(dirname(OUTDIR), basename(OUTDIR)), recursive = TRUE)
CONDIT_DIR <- paste(RES_FOLDER,"results",sep="_")
dir.create(file.path(OUTDIR,CONDIT_DIR), showWarnings = FALSE, recursive = TRUE)
FULL_OUTDIR <- paste(OUTDIR,CONDIT_DIR,"impulsede2",sep="/")
if (!dir.exists(FULL_OUTDIR)){
  dir.create(FULL_OUTDIR, recursive = TRUE)
} else {
  print("ImpulseDE2 directory already exists.")
}

cat("Output dir for ImpulseDE2 output and clustermap input: ", FULL_OUTDIR,"\n")

# Import libraries
library("ImpulseDE2")
library("org.Dm.eg.db")
library(annotate)
source("./scripts/geneID_converter.r", local=TRUE)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(data.table)
library(zoo)

# Set database to work with
gene_ID_database <- toTable(org.Dm.egFLYBASE)
gene_ID_database_name <- "flybase"

# Import metadata
metaData <- read.csv(file=METADATA)
setnames(metaData, old=c("ID","condition", "time", "batch"), new=c("Sample","Condition","Time","Batch"))
tail(metaData)

# Import uncorrected data
ct = read.csv(file= COUNTDATA, sep=",",row.names=1)
ct1 <- as.matrix(floor(ct))

# Run ImpulseDE2
objectImpulseDE2 <- runImpulseDE2(
  matCountData    = ct1, 
  dfAnnotation    = metaData,
  boolCaseCtrl    = FALSE,
  vecConfounders  = "Batch",
  boolIdentifyTransients = TRUE,
  scaNProc        = 30)

# Sort results by adjusted p-value
df <- objectImpulseDE2$dfImpulseDE2Results
dfIDE2_sortedpadj <- df[order(df$padj),] 
df_IDE2_filtpadj <- dfIDE2_sortedpadj[which(dfIDE2_sortedpadj$padj < PVAL_THRESH),]

# Add gene_names to df_IDE2_filtpadj
ids.type  <- gene_ID_database_name
id <- rownames(df_IDE2_filtpadj)
df_IDE2_filtpadj['gene_id']  <- id
df_IDE2_filtpadj$gene_name <- as.vector(get.symbolIDsDm(id,ids.type))
df_IDE2_filtpadj_sym <- as.data.frame(df_IDE2_filtpadj) %>% dplyr::select(gene_name, gene_id, everything()) 
to_drop <- c("Gene")
df_IDE2_filtpadj_sym_front<- df_IDE2_filtpadj_sym[ , !(names(df_IDE2_filtpadj_sym) %in% to_drop)]

# Create heatmap object
lsHeatmaps <- plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           = "case",
  boolIdentifyTransients = FALSE,
  scaQThres              = PVAL_THRESH) # FDR corrected p-val of significance

# Create gene trajectory dataframe
ids <- lsHeatmaps$lsvecGeneGroups[[1]]
res <- as.vector(get.symbolIDsDm(ids,ids.type))
traj_df4 <- as.data.frame(lsHeatmaps$complexHeatmapRaw@matrix)
traj_df4['gene_name'] <- res 
traj_df4['gene_id'] <- ids
traj_df3 <- traj_df4 %>% select(gene_id, gene_name, everything()) # move names to front
traj_df2 <- data.frame(t(apply(traj_df3,1,na.locf)), check.names=FALSE)
names(traj_df2) <- colnames(traj_df3)
traj_df <- as.data.frame(traj_df2) %>% select(gene_name, everything())

# Check if NA in gene_name, copy gene_id in its place
if (NA %in% traj_df$gene_name){
  print("NA in clustermap dataframe gene_name")
  traj_df$gene_name <- ifelse(is.na(traj_df$gene_name), traj_df$gene_id, traj_df$gene_name)
}
if ('NA' %in% traj_df$gene_name){
  print("NA in clustermap dataframe gene_name")
  traj_df$gene_name <- ifelse(is.na(traj_df$gene_name), traj_df$gene_id, traj_df$gene_name)
} 

# Output all spreadsheet results with adjusted p-values and z-score trajectory plots
write.csv(df_IDE2_filtpadj_sym_front,paste(FULL_OUTDIR,paste('impulsede2_output_padj',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # full data for all genes
write.csv(traj_df,paste(FULL_OUTDIR,paste('impulsede2_clustermapInput_padj',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
