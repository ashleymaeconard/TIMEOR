library(shiny, shinydashboard, DT)
library(ggplot2)
library(plotly)
library(devtools)
library(htmlwidgets)
library(dplyr)
library(Harman)
library(BiocManager)

# source("utils.R", local=T)
#################### Modules ###################

# metadataUI <- function(){
#   fileInput(
#     "metadataFile",
#     label = NULL,
#     accept = c(
#       "text/csv",
#       "text/comma-separated-values,text/plain",
#       ".txt"
#     )
#   )
# }


# import_metadata <- function(input, output, session, local_results_folder){
#   print("hi")
#   # Read in SraRunTable or Metadata input file from .csv
#   read_csv_sra_or_metadata_file <- reactive({
#     if(input$runDemo){
# 
#       # Set results folder to be demo data
#       file1 <- paste(local_results_folder(),"/timeor/data/metadata.csv",sep="")
#     }else {
# 
#       # Set input folder to be user defined
#       file1 <- input$metadataFile$datapath
#     }
#     if(is.null(file1)){return()}
#     read.csv(file = file1, sep = ",",header = TRUE)
#   })
# }
# 
# # Check metadata file is in correct format and then save as data frame
# df_sra_meta <- function(){
# 
#   # Read .csv (SraRunTable or Metadata file) and cast as dataframe
#   sra_metadata_df <- as.data.frame(read_csv_sra_or_metadata_file())
#   print(sra_metadata_df)
# 
#   # If metadata used
#   if(input$sra_or_meta == TRUE){
#     print("Metadata toggle")
# 
#     # Check format of metadata file
#     if(!(any(grepl("ID", colnames(sra_metadata_df), ignore.case = TRUE)) &
#          any(grepl("condition", colnames(sra_metadata_df), ignore.case = TRUE)) &
#          any(grepl("time", colnames(sra_metadata_df), ignore.case = TRUE)) &
#          any(grepl("batch", colnames(sra_metadata_df), ignore.case = TRUE)))){
#       shinyalert("Please fix 1) toggle 'Input file' to 'Metadata' in question 4 above.
#                  2) Make sure metadata file includes these columns:'id', 'condition', 'time', batch'.")
#     } else{
#       return(sra_metadata_df)
#     }
# 
#     # If SraRunTable used
#   }else{
#     print("SraRunTable toggle")
# 
#     # Check format of SraRunTable file
#     if(!(any(grepl("treatment", colnames(sra_metadata_df), ignore.case = TRUE)) &
#          any(grepl("time", colnames(sra_metadata_df), ignore.case = TRUE)) &
#          any(grepl("replicate", colnames(sra_metadata_df), ignore.case = TRUE)) &
#          any(grepl("batch", colnames(sra_metadata_df), ignore.case = TRUE)) &
#          any(grepl("Run", colnames(sra_metadata_df), ignore.case = TRUE)))){
#       shinyalert("Please fix 1) toggle 'Input file' to 'SraRunTable' in question 4 above.
#                  2) Make sure SraRunTable file includes these columns:'run', 'treatment', 'replicate', 'time', batch'.")
# 
#     } else{
#       parsed_sra_to_metadata_df <- parse_sra(sra_metadata_df)
#       return(parsed_sra_df)
#     }
# 
#   }
# }

#################### Individual Functions ###################
# Bookmarking user location
enableBookmarking(store = "server")

# Run upper quartile normalization 
# Two methods of upper quartile determination are presented here, one 
# which is based on the total counts, and one based on the expressed 
# counts (i.e. excluding values with no expression).
upperNormalization <- function(cD){
  # Filter low counts with mean less than 5
  cD <- cD[rowMeans(cD) > 5,];
  cD.upquantileAll <- apply(cD, 2, function(x){quantile(x, 0.75)});
  cD.upquantileExpressed <- apply(cD, 2, function(x){quantile(x[x>0], 0.75)});
  cD.norm <- t(t(cD) / cD.upquantileAll);
  return(cD.norm)
}

# Run Harman Correction
harmanCorrection <- function(nD, mD){
  harman_corr <- harman(nD, expt = mD$time, batch=mD$batch, limit=0.99)
  return(harman_corr)
}










