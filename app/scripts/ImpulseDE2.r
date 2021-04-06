# Ashley Conard
# ImpulseDE2.r
# Last Mod. 12/12/2019
# Purpose: Run ImpulseDE2 to identify differentially expressed genes from time series data.

# Checking input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type:ImpulseDE2.r /FULL/PATH/TO/METDATA_FILE/ /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) /FULL/PATH/TO/OUTPUTDIR/ RES_FOLDER (e.g. insulin_stim) PVAL_THRESH ORGANISM", call.=FALSE)
} else if (length(args) == 6) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 6 arguments. Type:ImpulseDE2.r /FULL/PATH/TO/METDATA_FILE/ /FULL/PATH/TO/COUNT_MATRIX_FILE/ (not normalized and corrected) /FULL/PATH/TO/OUTPUTDIR/ RES_FOLDER (e.g. insulin_stim) PVAL_THRESH ORGANISM")
}

# Input arguments
METADATA <- args[1] 
COUNTDATA <- args[2] 
OUTDIR <- args[3] 
RES_FOLDER <- args[4]
PVAL_THRESH <- as.numeric(args[5])
ORGANISM <- args[6]

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
library(annotate)
source("./scripts/geneID_converter.r", local=TRUE)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(data.table)
library(zoo)

# Assigning organism library
if(ORGANISM=="dme"){
    ORG_DB="org.Dm.eg.db"
} else if(ORGANISM=="hsa"){
    require(biomaRt)
    ORG_DB="org.Hs.eg.db"
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
}else if(ORGANISM=="mmu"){
    require(biomaRt)
    ORG_DB="org.Mm.eg.db"
    mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
} else{
    stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
}
library(ORG_DB, character.only = TRUE) # organism database library

# Import metadata
metaData <- read.csv(file=METADATA)
setnames(metaData, old=c("ID","condition", "time", "batch"), new=c("Sample","Condition","Time","Batch"))
head(metaData)

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
gene_id_list <- sub('\\.[0-9]*$', '', rownames(df_IDE2_filtpadj))
df_IDE2_filtpadj['gene_id']  <- gene_id_list

if(ORGANISM=="dme"){
  # Set database to work with
  gene_ID_database <- toTable(org.Dm.egFLYBASE)
  ids.type <- "flybase"  
  
  df_IDE2_filtpadj$gene_name <- as.vector(get.symbolIDsDm(gene_id_list,ids.type))
  df_IDE2_filtpadj_sym <- as.data.frame(df_IDE2_filtpadj) %>% dplyr::select(gene_name, gene_id, everything()) 
  to_drop <- c("Gene")
  df_IDE2_filtpadj_sym_front<- df_IDE2_filtpadj_sym[ , !(names(df_IDE2_filtpadj_sym) %in% to_drop)]

} else if(ORGANISM=="hsa" || ORGANISM=="mmu"){
  bm_ids <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = gene_id_list, mart = mart)
  df_IDE2_filtpadj_tmp <- as.data.frame(merge(df_IDE2_filtpadj, bm_ids, by.x="gene_id", by.y="ensembl_gene_id"))
  df_IDE2_filtpadj_tmp_gene_name <- df_IDE2_filtpadj_tmp%>%plyr::rename(c(external_gene_name ="gene_name"))
  df_IDE2_filtpadj_sym <- df_IDE2_filtpadj_tmp_gene_name %>% dplyr::select(gene_name, gene_id, everything())
  to_drop <- c("Gene")
  df_IDE2_filtpadj_sym_front<- df_IDE2_filtpadj_sym[ , !(names(df_IDE2_filtpadj_sym) %in% to_drop)]
  
} else{
    stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
}

# Create heatmap object
lsHeatmaps <- plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           = "case",
  boolIdentifyTransients = FALSE,
  scaQThres              = PVAL_THRESH) # FDR corrected p-val of significance

# Create gene trajectory dataframe
ids <- lsHeatmaps$lsvecGeneGroups[[1]] # gene_ids from heatmap 
# Create result gene_name column for gene_ids (not sorted!)
if(ORGANISM=="dme"){
  res <- as.vector(get.symbolIDsDm(ids,ids.type))
  traj_df4_tmp2 <- as.data.frame(lsHeatmaps$complexHeatmapRaw@matrix)
  traj_df4_tmp2['gene_name'] <- res
  traj_df4_tmp2['gene_id'] <- ids
   
} else if(ORGANISM=="hsa" || ORGANISM=="mmu"){
  bm_ids <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = ids, mart = mart)
  gene_id_name <- merge(data.frame(ids),data.frame(bm_ids), by.x="ids", by.y="ensembl_gene_id", sort = F) # do not sort rows with sort=F
  res <- as.vector(gene_id_name$external_gene_name)
  traj_df4 <- as.data.frame(lsHeatmaps$complexHeatmapRaw@matrix)
  traj_df4['gene_id'] <- ids
  traj_df4_tmp <- merge(traj_df4, gene_id_name, by.x="gene_id", by.y = "ids")
  traj_df4_tmp2 <- traj_df4_tmp%>%plyr::rename(c(external_gene_name="gene_name")) 

} else{
    stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
}

traj_df3 <- traj_df4_tmp2 %>% dplyr::select(gene_id, gene_name, everything()) # move names to front
traj_df2 <- data.frame(t(apply(traj_df3,1,na.locf)), check.names=FALSE)
names(traj_df2) <- colnames(traj_df3)
traj_df <- as.data.frame(traj_df2) %>% dplyr::select(gene_name, everything())

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
