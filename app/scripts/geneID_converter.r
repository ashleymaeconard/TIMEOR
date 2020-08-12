# geneID_converter.r
# Ashley Conard
# Last Modified: Aug. 9, 2019
# Resource: function get.symbolIDsDm TAKEN FROM https://www.researchgate.net/publication/308990864_R_function_getsymbolIDsDm_to_convert_Uniprot_Flybase_etc_to_gene_symbol


# args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#   stop("Type: /usr/local/bin/Rscript clusterProfiler.r /PATH/TO/INPUT_OUTPUT_DIR/ (Pass in Input/Output Directory - will be the same - e.g. /Users/ashleymaeconard/Desktop/RESULTS/Feb4_min5_clusters/) OVERLAP_NOT_EXPERIMENTS (set to 1 to run overlap comparison, 0 otherwise)", call.=FALSE)
# } else if (length(args) == 2) {
#     cat("Passed in:", args,"\n")
# } else{
#     stop("Pass in 1) Input/Output Directory (e.g. /Users/ashleymaeconard/Desktop/RESULTS/Feb4_min5_clusters/) 2) OVERLAP_NOT_EXPERIMENTS (set to 1 to run overlap comparison, 0 otherwise)")
# }

get.symbolIDsDm <- function(id,id.type){
  
  # ************************************************
  # get.symbolIDs function programmed by 
  # Benjamin Tovar | February 25, 2014
  # http://tata-box-blog.blogspot.de/2014/02/convert-ensembl-unigene-uniprot-and.html
  
  # modified to get.symbolIDsDm function to work with Drosophila melanogaster by
  # Christoph Metzendorf | October 11, 2016
  
  # INSTRUCTIONS
  # id = vector of the original IDs, for example: 
  # 	c("FBgn0040373","FBgn0040372","FBgn0261446" )
  # id.type = type of the original IDs, in the example is ensembl
  
  # NOTE: only refseq, ensembl, uniprot, unigene and flybase are supported,
  # this function depends on the Bioconductor package org.Dm.eg.db
  # that can be installed:
  
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("org.Dm.eg.db")
  
  # ************************************************
  # # USAGE EXAMPlE: ENSEMBL
  # 	require(org.Dm.eg.db)
  # 	ensembl <- toTable(org.Dm.egENSEMBL)
  # 	id  <- ensembl[1:100,2]
  # 	id.type  <- "ensembl"
  # 	res <- get.symbolIDsDm(id,id.type)
  
  # # USAGE EXAMPlE: UNIPROT
  #   require(org.Dm.eg.db)
  # 	uniprot <- toTable(org.Dm.egUNIPROT)
  # 	id  <- uniprot[1:100,2]
  # 	id.type  <- "uniprot"
  # 	res <- get.symbolIDsDm(id,id.type)
  
  # # USAGE EXAMPlE: REFSEQ
  # require(org.Dm.eg.db)
  # refseq.id <- toTable(org.Dm.egREFSEQ)
  # id  <- refseq.id[1:100,2]
  # id.type  <- "refseq"
  # res <- get.symbolIDsDm(id,id.type)	
  
  # # USAGE EXAMPlE: UNIGENE
  # require(org.Dm.eg.db)
  # unigene <- toTable(org.Dm.egUNIGENE)
  # id  <- unigene[1:100,2]
  # id.type  <- "unigene"
  # res <- get.symbolIDsDm(id,id.type)	
 
  # # USAGE EXAMPlE: FLYBASE
  # require(org.Dm.eg.db)
  # flybase <- toTable(org.Dm.egFLYBASE)
  # id  <- flybase[1:100,2]
  # id.type  <- "flybase"
  # res <- get.symbolIDsDm(id,id.type)
  
  
  # LOAD THE ANNOTATION INFORMATION
  cat("Note: Running function get.symbolIDsDm from Christoph Metzendorf.")
  cat("1) Loading annotation library",date(),"\n")   	
  require(org.Dm.eg.db)
  cat("2) Loading annotation for symbol IDs",date(),"\n")   	    
  symbol <- toTable(org.Dm.egSYMBOL)
  # IF THE ORIGINALS IDS = ENSEMBL
  if(id.type=="ensembl"){
    cat("3) Loading annotation for ensembl IDs",date(),"\n")   	
    annotation <- toTable(org.Dm.egENSEMBL)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }
  # IF THE ORIGINALS IDS = UNIPROT
  if(id.type=="uniprot"){
    cat("3) Loading annotation for uniprot IDs",date(),"\n")   
    annotation <- toTable(org.Dm.egUNIPROT)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }
  # IF THE ORIGINALS IDS = REFSEQ
  if(id.type=="refseq"){
    cat("3) Loading annotation for refseq IDs",date(),"\n")   
    annotation <- toTable(org.Dm.egREFSEQ)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }
  # IF THE ORIGINALS IDS = UNIGENE
  if(id.type=="unigene"){
    cat("3) Loading annotation for unigene IDs",date(),"\n")   
    annotation <- toTable(org.Dm.egUNIGENE)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }	
    if(id.type=="flybase"){
    cat("3) Loading annotation for flybase IDs",date(),"\n")   
    annotation <- toTable(org.Dm.egFLYBASE)
    # extract the indexes of the original database 
    index <- as.numeric(sapply(id, function(x) which(annotation[,2]==x)[1]))
    # extract the gene_ids
    gene_id.index <- as.numeric(annotation[index,1])
    # parse the indexes back to the symbol IDs 
    index <- as.numeric(sapply(gene_id.index, function(x) which(symbol[,1]==x)))
    # extract the IDs
    symbolIDs <- symbol[index,2]
    # export the output
    return(symbolIDs)
  }	
  cat("** ERROR: DATABASE TYPE NOT SELECTED | TRY: ensembl OR unigene OR uniprot OR refseq **",date(),"\n") 
  return(NA)
}