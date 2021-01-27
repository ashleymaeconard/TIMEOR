# DESeq2.r
# Ashley Mae Conard
# Last Mod. 12/12/2019
# Purpose: Run DESeq2 to identify differentially expressed genes from time series data.
# Resources: http://seqanswers.com/forums/showthread.php?t=64039
#            http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#            https://lashlock.github.io/compbio/R_presentation.html
#            http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html
#            http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2
#            http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments    
#            https://support.bioconductor.org/p/101613/

# Importing libraries
library("DESeq2") # make sure 1.22
library(dplyr)
source("./scripts/geneID_converter.r", local=TRUE)
library(data.table)

run_DESeq2 <- function(METDATA, COUNTDATA, OUTPUTDIR, CONDITION, BATCH_EFFECT, TC, PVAL_THRESH, CONTROL, ORGANISM){
    
    # Assigning organism library
    if(ORGANISM=="dme"){
        ORG_DB="org.Dm.eg.db"
    } else if(ORGANISM=="hsa"){
        ORG_DB="org.Hs.eg.db"
    }else if(ORGANISM=="mmu"){
        ORG_DB="org.Mm.eg.db" 
    } else{
        stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
    }
    library(ORG_DB, character.only = TRUE) # organism database library
    
    # Set database to work with
    gene_ID_database <- toTable(ORG_DB)
    if(ORGANISM=="dme"){
        gene_ID_database_name <- "flybase"
    } else if(ORGANISM!="dme"){
        gene_ID_database_name <- "NA"
    } else{
        stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
    }

    # Create or set output directory
    CONDIT_DIR <- paste(CONDITION,"results",sep="_")
    FULL_OUTDIR <- paste(OUTPUTDIR,CONDIT_DIR,"deseq2",sep="/")
    if (!dir.exists(CONDIT_DIR)){
        dir.create(CONDIT_DIR)
    } else {
        print("Results directory already exists.")
    }
    if (!dir.exists(FULL_OUTDIR)){
        dir.create(FULL_OUTDIR)
    } else {
        print("DESeq2 directory already exists.")
    }
    cat("Output directory: ", FULL_OUTDIR)

    # Load data and metadata
    # Import count data
    countData1 <- read.csv(COUNTDATA, header=TRUE, sep=",")
    col_name = "ID"
    
    # Import metadata
    metaData <- read.csv(METDATA, header=TRUE, sep=",")
    col_for_index="ID"

    # Format metadata
    keep1 <- as.vector(metaData[[col_for_index]]) # conditions to keep
    rownames(metaData) <- metaData$ID # set rownames to be condition labels
    drop<-c(col_for_index)
    metaData = metaData[,!(names(metaData) %in% drop)]
    metaData$time <- as.factor(metaData$time) # use factor for Time (which is all integers)
    metaData$condition <- as.factor(metaData$condition)
    metaData$batch <- as.factor(metaData$batch)
    metaData <- metaData[ order(row.names(metaData)), ]
    
    # Checking to see if case vs. control or just case or control
    if (dim(table(metaData$condition)) == 1){ 
        CASEvsCONT <- 0  
    } else{
        CASEvsCONT <- 1 # yes this is case vs. control of some sort
        
        # Checking case vs. control (i.e. do all timepoints have a control or is it one control for all other timepoints)
        freq_case_contr = data.frame(table(metaData$condition))
    }
    print(head(metaData))

    # Make sure all Flybase ID names are unique
    n_occur <- data.frame(table(countData1[col_name]))
    non_unique = n_occur[n_occur$Freq > 1,]
    if(length(non_unique) == 0){
        print("ERROR - COL. NAMES (FLYBASE) NOT ALL UNIQUE")
        print(nrow(countData1))
        print(non_unique)
        print(countData1[countData1[col_name] %in% n_occur$Var1[n_occur$Freq > 1],])
        quit()
    }
    rownames(countData1) <- countData1[[col_name]]# use the Flybase name as unique ID
    countData1 <- countData1[, !(names(countData1) %in% col_name)] # remove Flybase name from columns

    # Subsetting merged_htseq based on metadata file
    ids_to_keep <- c(rownames(metaData))
    print("ids_to_keep")
    print(ids_to_keep)
    col.num <- which(colnames(countData1) %in% ids_to_keep)
    countData <- countData1[,col.num]
    countData <- countData[ , order(names(countData))]
    print(head(countData))

    # Check that all row names for metaData match column names for countData
    if(!all(rownames(metaData) %in% colnames(countData))){
        write("ERROR: row names for metaData do not match column names for countData.", stderr())
    }

    # Create DESeq2 object and run DESeq2
    # LRT tests whether the terms removed in reduced model explain a significant amount of variation in the input data. It is recommended for time series - http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
    # LRT is recommended without betaPrior - https://support.bioconductor.org/p/87224/
    
    # 3 variables
    if(BATCH_EFFECT & TC & CASEvsCONT){ # batch effect, timecourse, and case vs. control
        write("Batch effect, timecourse, and case vs. control", stderr())
        if(!any(freq_case_contr$Freq >= floor(dim(metaData)[1]/2))){
            write("Batch effect, timecourse with timepoint t compared to t-1", stderr())
            dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch + condition + time + condition:time) # examine effect of condition over time
            
            # Set reference level for case control comparison
            dds$condition <- relevel(dds$condition, ref = CONTROL)
           
            # Run DESeq2
            dds <- DESeq(dds, test="LRT", reduced= ~ batch + condition + time)
        } else{ # Reference: http://52.71.54.154/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html#time
            write("Batch effect and timecourse with control as 1st timepoint", stderr())
            dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch + time)# examine effect of time
            
            # Set reference level for case control comparison
            dds$condition <- relevel(dds$condition, ref = CONTROL)

            # Run DESeq2            
            dds <- DESeq(dds, test="LRT", reduced= ~ batch)
        }
    # 2 variables
    } else if(!BATCH_EFFECT & TC & CASEvsCONT){ # timecourse and condition (i.e. case vs. control)
        write("Timecourse and case vs. control", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ condition + time + condition:time)# examine effect of condition over time
        dds <- DESeq(dds, test="LRT", reduced= ~ condition + time)
    } else if(BATCH_EFFECT & !TC & CASEvsCONT){ # batch effect and condition
        write("Batch effect and case vs. control", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch + condition) # examine effect of condition
        dds <- DESeq(dds, test="Wald", betaPrior=TRUE, reduced= ~ batch) ## add Wald test here
    } else if(BATCH_EFFECT & TC & !CASEvsCONT){ # batch effect and timecourse
        write("Batch effect and timecourse", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch + time)# examine effect of time
        dds <- DESeq(dds, test="LRT", reduced= ~ batch)
    # 1 variable
    } else if(!BATCH_EFFECT & !TC & CASEvsCONT){ # condition (normally what DESeq2 is used for - case vs. control, no time)
        write("Case vs. control", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ condition)
        dds <- DESeq(dds, test="Wald", betaPrior=TRUE, reduced= ~1)
    } else if(!BATCH_EFFECT & TC & !CASEvsCONT){ # condition (normally what DESeq2 is used for - case vs. control, no time)
        write("Timecourse", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ time)
        dds <- DESeq(dds, test="LRT", reduced= ~1)
    } else if(BATCH_EFFECT & !TC & !CASEvsCONT){ # condition (normally what DESeq2 is used for - case vs. control, no time)
        write("Batch effect", stderr())
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ batch)
        dds <- DESeq(dds, test="LRT", reduced= ~1)
    } else{
        write("ERROR: TIMEOR processes datasets with columns 'time', 'batch', and 'condition'.", stderr())
    }

    # Filtering based on adjusted p-value (i.e. padj)
    res_no_padj <- results(dds)
    res <- res_no_padj[which(res_no_padj$padj < PVAL_THRESH),]
    
    # Performing apeglm shrinkage transformation
    # Resource: https://rdrr.io/bioc/DESeq2/man/lfcShrink.html
    resApe_no_padj <- lfcShrink(dds, type="apeglm")
    resApe <- resApe_no_padj[which(resApe_no_padj$padj < PVAL_THRESH),]

    # Saving MA plots before and after shrinkage
    pdf(paste(FULL_OUTDIR,paste('ma_plot_noShrinkage_padj',toString(PVAL_THRESH),'.pdf',sep=""),sep="/"))
    plotMA(res, ylim=c(-4,4), cex=.8)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()
    pdf(paste(FULL_OUTDIR,paste('ma_plot_apeShrinkage_padj',toString(PVAL_THRESH),'.pdf',sep=""), sep="/"))
    plotMA(resApe, ylim=c(-4,4), cex=.8)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    dev.off()

    # Adding gene symbol and placing it in the front for no and apeglm shrinkage matricies
    ids.type  <- gene_ID_database_name
    ids <- rownames(res)
    idsN <- rownames(resApe)
    
    res['gene_id'] <- rownames(res)
    res$gene_name <- as.vector(get.symbolIDsDm(ids,ids.type))
    res_sym_front <- as.data.frame(res) %>% dplyr::select(gene_name, gene_id, everything())    
    
    resApe['gene_id'] <- rownames(resApe)
    resApe$gene_name <- as.vector(get.symbolIDsDm(idsN,ids.type))
    resN_sym_front <- as.data.frame(resApe) %>% dplyr::select(gene_name, gene_id, everything())    

    # Creating clustermap inputs for no and apeglm shrinkage
    betasTC <- coef(dds)
    colnames(betasTC)
    topGenes <- which(res_sym_front$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    topGenesN <- which(resN_sym_front$padj < PVAL_THRESH, arr.ind = FALSE) # get indicies for all results
    
    # Generating clustermap no shrinkage matrix (NOTE 1,2,3 removed the batch effect columns)
    batch_cols <- seq(1,length(unique(metaData$time))-1)
    mat <- betasTC[topGenes, -(batch_cols)]
    write("batch_cols", stderr())
    write(batch_cols, stderr())
    
    # Generating clustermap apeglm shrinkage matrix (NOTE 1,2,3 removed the batch effect columns)
    batch_cols <- seq(1,length(unique(metaData$time)))
    matN <- betasTC[topGenesN, -batch_cols]

    # Adding gene symbol and placing it in the front for clustermap input matricies (no shrinkage)
    df_mat <- as.data.frame(mat)
    idsM <- rownames(df_mat)
    df_mat['gene_id'] <- idsM
    df_mat$gene_name <- as.vector(get.symbolIDsDm(idsM,ids.type))
    df_mat_sym_front <- as.data.frame(df_mat) %>% dplyr::select(gene_name, gene_id, everything()) 
    
    # Check if NA in gene_name, copy gene_id in its place
    if (NA %in% df_mat_sym_front$gene_name){
        df_mat_sym_front$gene_name <- ifelse(is.na(df_mat_sym_front$gene_name), df_mat_sym_front$gene_id, df_mat_sym_front$gene_name)
    }

    # Adding gene symbol and placing it in the front for clustermap input matricies (apeglm shrinkage)
    df_matN <- as.data.frame(matN)
    idsMN <- rownames(df_matN)
    df_matN['gene_id'] <- idsMN
    df_matN$gene_name <- as.vector(get.symbolIDsDm(idsMN,ids.type))
    df_matN_sym_front <- as.data.frame(df_matN) %>% dplyr::select(gene_name, gene_id, everything()) 

    # Check if NA in gene_name, copy gene_id in its place
    if (NA %in% df_matN_sym_front$gene_name){
        df_matN_sym_front$gene_name <- ifelse(is.na(df_mat_sym_front$gene_name), df_mat_sym_front$gene_id, df_mat_sym_front$gene_name)
    }

    # Saving sorted (by padj) results
    resSort <- res_sym_front[order(res_sym_front$padj),]
    resSortApeShr <- resN_sym_front[order(resN_sym_front$padj),]
    
    write.csv(as.data.frame(resSort), file=paste(FULL_OUTDIR,paste("deseq2_output_noShrinkage_padj",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(as.data.frame(resSortApeShr), file=paste(FULL_OUTDIR,paste("deseq2_output_ApeShrinkage_padj",toString(PVAL_THRESH),'.csv',sep=""), sep="/"), row.names=FALSE, quote=FALSE)
    write.csv(na.omit(df_mat_sym_front),paste(FULL_OUTDIR,paste('deseq2_noShrinkage_clustermapInput_padj',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
    write.csv(na.omit(df_matN_sym_front),paste(FULL_OUTDIR,paste('deseq2_ApeShrinkage_clustermapInput_padj',toString(PVAL_THRESH),'.csv',sep=""),sep='/'), row.names=FALSE, quote=FALSE) # clustermap input
}  
