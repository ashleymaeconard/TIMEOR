# Ashley Conard
# venns_intervene.r
# Last Mod. 12/12/2019

# Input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Pass in 5 arguments. Type: Rscript venns.r 
    1) /FULL/PATH/TO/IMPULSEDE2_OUTPUT_FILE/ 
    2) /FULL/PATH/TO/DESEQ2_OUTPUT_FILE/ 
    3) /FULL/PATH/TO/NEXTMASIGPRO_OUTPUT_FILE/
    4) /FULL/PATH/TO/PAST_STUDY_LIST_GENES 
    5) PAST_STUDY_NAME 
    6) /FULL/PATH/TO/OUTPUTDIR", call.=FALSE)
} else if (length(args) == 6) {
  cat("Passed in:", args,"\n")
} else{
  stop("Pass in 5 arguments. Type: Rscript venns.r 
    1) /FULL/PATH/TO/IMPULSEDE2_OUTPUT_FILE/ 
    2) /FULL/PATH/TO/DESEQ2_OUTPUT_FILE/ 
    3) /FULL/PATH/TO/NEXTMASIGPRO_OUTPUT_FILE/
    4) /FULL/PATH/TO/PAST_STUDY_LIST_GENES 
    5) PAST_STUDY_NAME 
    6) /FULL/PATH/TO/OUTPUTDIR")
}

# Parameters
IDE2 <- args[1] # ImpulseDE2 output file
DES2 <- args[2] # DESeq2 output file
NMSP <- args[3] # nextMaSigPro output file
PS <- args[4] # .txt file will be imported as a list where each gene is separated by a new line
PS_NAME <- args[5] # previous study name
OUTPUTDIR <- args[6] # output dir

# Create output directory if needed
dir.create(file.path(file.path(dirname(OUTPUTDIR),basename(OUTPUTDIR)))) 

# Load ImpulseDE2 and create gene list
IDE2_gene_df <- read.table(file = IDE2, sep = ",", header=TRUE)[c("gene_id")]
write.table(IDE2_gene_df, file=paste(OUTPUTDIR,"ImpulseDE2_gene_list.csv", sep="/"),  quote=FALSE, row.names=FALSE, col.names=FALSE)

# Load DESeq2 and create gene list
DES2_gene_df <- read.table(file = DES2, sep = ",", header=TRUE)[c("gene_id")]
write.table(DES2_gene_df, file=paste(OUTPUTDIR,"DESeq2_gene_list.csv", sep="/"),  quote=FALSE, row.names=FALSE, col.names=FALSE) 

# Load NextMaSigPro and create gene list
NMSP_gene_df <- read.table(file = NMSP, sep = ",", header=TRUE)[c("gene_id")]
write.table(NMSP_gene_df, file=paste(OUTPUTDIR,"nextMaSigPro_gene_list.csv", sep="/"),  quote=FALSE, row.names=FALSE, col.names=FALSE) 

# Check for empty string for previous study
if (PS == 'NA'){
  print("No previous study to compare.")
  
  # 3 Venn diagram overlap 
  command <- paste("intervene venn","--type", "list", "-i", paste(OUTPUTDIR,"ImpulseDE2_gene_list.csv", sep="/"), paste(OUTPUTDIR,"nextMaSigPro_gene_list.csv", sep="/"), paste(OUTPUTDIR,"DESeq2_gene_list.csv", sep="/"), "--figtype png", "--names=ImpulseDE2,Next_maSigPro,DESeq2" ,"--save-overlaps -o",OUTPUTDIR,sep=" ")
  system(command)

} else{ 
  print(paste("Previous study ",PS," to compare.", sep=""))
  
  # 4 Venn diagram overlap 
  command <- paste("intervene venn","--type", "list", "-i", paste(OUTPUTDIR,"ImpulseDE2_gene_list.csv", sep="/"), paste(OUTPUTDIR,"nextMaSigPro_gene_list.csv", sep="/"), paste(OUTPUTDIR,"DESeq2_gene_list.csv", sep="/"), PS, "--figtype png", paste("--names=ImpulseDE2,Next_maSigPro,DESeq2,", PS_NAME, sep=""), "--save-overlaps -o", OUTPUTDIR,sep=" ")
  system(command)

}
