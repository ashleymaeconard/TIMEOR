# Ashley Mae Conard
# clusterProfiler_pathview_indiv_cluster.r
# Last Modified: Jan 15, 2019
# Runs clusterProfiler and pathview to find and visualize gene enrichment within each cluster.
# References: https://yulab-smu.github.io/clusterProfiler-book/chapter12.html          

# CHECKING ARGUMENTS 
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type: /usr/local/bin/Rscript clusterProfiler.r 
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/<EXPERIM_results/clusters/1/) 
            2) EXPERIMENT_NAME (e.g. test_results, write 'test') 
            3) SEPARATE_TIMEPOINTS (set to 1 to run GO for each timepoint separately, 0 otherwise) 
            4) ORGANISM (dme, hsa, mmu) 
            5) ADJ_PVAL (recommend 0.05)", call.=FALSE)
} else if (length(args) == 5) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/<EXPERIM_results/clusters/1/) 
                  2) EXPERIMENT_NAME (e.g. test_results, write 'test') 
                  3) SEPARATE_TIMEPOINTS (set to 1 to run GO for each timepoint separately, 0 otherwise) 
                  4) ORGANISM (dme (Drosophila melanogaster), hsa (Homo sapiens), mmu (Mus musculus)
                  5) ADJ_PVAL (recommend 0.05))")
}

# ASSIGNING INPUT ARGUMENTS
IN_OUTPUT <- args[1] #e.g. ~/timeor/results/primary/
EXPERIMENT_NAME <- args[2] # without "_results" added
SEP_TPS <- args[3] # 0 or 1 
ORGANISM <- args[4] # dme (Drosophila melanogaster), hsa (Homo sapiens), mmu (Mus musculus)
options(digits=5)
ADJ_PVAL <- as.double(args[5]) # recommend 0.05

# ASSIGNING ORGANISM LIBRARY
if(ORGANISM=="dme"){
    ORG_DB="org.Dm.eg.db"
} else if(ORGANISM=="hsa"){
    ORG_DB="org.Hs.eg.db"
}else if(ORGANISM=="mmu"){
    ORG_DB="org.Mm.eg.db"
} else{
    stop("Please enter dme (Drosophila melanogaster), hsa (Homo sapiens), or mmu (Mus musculus)")
}

# LOADING PACKAGES
library(biomaRt)
library(devtools)
#install_github("GuangchuangYu/bitr")
#library(bitr)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(ORG_DB, character.only = TRUE) # organism database library
library("pathview")
library(cowplot)

write("Loaded Organism: ", stderr())
write(ORG_DB, stderr())

assess_enrichment <- function(geneENS, type_enr){

    # GO over-representation (enrichement) test 
    ego <- enrichGO(gene        = geneENS,
                    OrgDb         = ORG_DB,
                    keyType       = "ENSEMBL",
                    ont           = type_enr,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = ADJ_PVAL,
                    qvalueCutoff  = ADJ_PVAL)                
    ego <- setReadable(ego, OrgDb = ORG_DB)
    return(ego)
}

mapNameToPathwayID <- function(organism, enriched_name) { # function alterned from https://biobeat.wordpress.com/category/r/  
  # Getting the pathway URL for the input organism
  KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
  pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, organism, sep="")
  
  # Setting up the table to fill with the enriched_name, if any
  pathway_id_name <- data.frame()
  
  # Iterating through pathway list to find enriched_term and return the pathway ID
  for (line in readLines(pathway_list_REST_url)) {
    tmp <- strsplit(line, "\t")[[1]]
    pathway_id <- strsplit(tmp[1], organism)[[1]][2]
    pathway_name <- tmp[2]
    pathway_name <- strsplit(pathway_name, "\\s+-\\s+")[[1]][1]
    
    # If the enriched_name is part of the pathway name, get the pathway ID
    if(grepl(enriched_name, pathway_name, ignore.case=TRUE)){
        pathway_id_name[pathway_name, 1] = pathway_id
    }
  }
  # Returning the pathway name and ID
  pathway_id_name
}

plotPathway <- function(gl, exper, pathwayName, subdi, num_g_clust){
    print("Assessing pathway enrichment.")
    # Getting pathway ID
    pathwayID <- mapNameToPathwayID(ORGANISM,pathwayName)
    print(paste("Fetched pathway ID:", pathwayID, sep=" "))
    
    if(!length(pathwayID)){
        print(paste("Pathway ID is empty for", pathwayName, sep=" "))
    } else {

        # Generate www folder within each cluster subdir for png rendering in R Shiny
        wwwsvg_subdirec <- file.path(subdi, "www") 
        if (!dir.exists(wwwsvg_subdirec)){
            dir.create(wwwsvg_subdirec)
        } else {
            cat(wwwsvg_subdirec, "www/ subdirectory for *.png multi-timepoint pathway plot exists.")
        }

        # Plotting pathway (no significance limit because searching based on pvalue thresholded clusterProfiler GO terms)
        pathview(gene.data = gl,
                            pathway.id = paste(ORGANISM,pathwayID[,1],sep=""),
                            species = ORGANISM, 
                            out.suffix = paste(num_g_clust,"_pathview",sep="") )
        
        # Moving automatic output from Pathview to correct output folder 
        pathview_higlighted_genes_file_name <- paste(paste(ORGANISM,pathwayID[,1],sep=""),paste(num_g_clust,"pathview.multi.png",sep="_"),sep=".")
        file.rename(paste(getwd(),pathview_higlighted_genes_file_name,sep="/" ),  
                    paste(wwwsvg_subdirec, pathview_higlighted_genes_file_name, sep="/"))

        pathview_all_genes_file_name <- paste(ORGANISM,pathwayID[,1],".png",sep="")
        file.rename(paste(getwd(),pathview_all_genes_file_name,sep="/" ),  
                    paste(subdi, pathview_all_genes_file_name, sep="/"))

        pathview_xml_file_name <- paste(ORGANISM,pathwayID[,1],".xml",sep="")
        file.rename(paste(getwd(),pathview_xml_file_name,sep="/" ),  
                    paste(subdi, pathview_xml_file_name, sep="/"))
    }
}

generate_plots <- function(egoo, gl, ex, type_enrich, outdir, num_gene_cluster){
    # Checking if the enrichment object is empty, and if so, exit function and script
    if(nrow(egoo)>0){

        # Generating subfolder within specific cluster 
        subdirec <- file.path(outdir, type_enrich) 
        if (!dir.exists(subdirec)){
            dir.create(subdirec)
        } else {
            cat(subdirec, "subdirectory exists.")
        }
        # Generate www folder within each subfolder for dotplot rendering in R Shiny
        wwwsvg_subdirec <- file.path(subdirec, "www") 
        if (!dir.exists(wwwsvg_subdirec)){
            dir.create(wwwsvg_subdirec)
        } else {
            cat(wwwsvg_subdirec, "www/ subdirectory for *.svg dotplot exists.")
        }

        # Plotting GOrilla plot
        write.csv(egoo, paste(subdirec, paste(ex,"cluster_", num_gene_cluster, "clustProf",type_enrich, ".csv", sep="_"), sep="/"))
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_dag",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_plotGOgraph <- plotGOgraph(egoo)
        print(plt_plotGOgraph)
        dev.off()

        # Plotting dot plot
        svg(paste(paste(wwwsvg_subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_dotplot",paste(type_enrich,".svg",sep=""), sep="_"), sep="/")))
        plt_dotplot <- dotplot(egoo)
        print(plt_dotplot)
        dev.off()

        # Plotting clusterProfiler relationship graph plot
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_emaplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_emaplot <- emapplot(egoo) 
        print(plt_emaplot)
        dev.off()

        # Barplot 
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_barplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
        plt_bar<- barplot(egoo, showCategory=10)
        print(plt_bar)
        dev.off()

        # Gene concept network
        pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_conceptplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 30, height = 15)
        p2 <- cnetplot(egoo, categorySize="pvalue", foldChange=gl)
        plt_cnet <- cnetplot(egoo, foldChange=gl, categorySize="pvalue", colorEdge = TRUE) 
        plt_cnet_cir <- cnetplot(egoo, foldChange=gl, circular = TRUE, categorySize="pvalue", colorEdge = TRUE)
        plt_concept <- plot_grid(plt_cnet, plt_cnet_cir, ncol=2)
        print(plt_concept)
        dev.off()

        # Plotting clusterProfiler phylogeny plot (if there is more than 1 enriched GO term)
        if(dim(egoo)[1]>1){
            pdf(paste(paste(subdirec, paste(ex, "geneCluster",num_gene_cluster, "clustProf_goplot",paste(type_enrich,".pdf",sep=""), sep="_"), sep="/")), width = 10, height = 10)
            plt_goplot <- goplot(egoo)
            print(plt_goplot)
            dev.off()
        }
    } else{
        print(paste(type_enrich," is empty", sep=""))
    }
}

main <- function(){
    # Finding only filenames that are *geneList*.csv or *geneName*.csv to get list of genes to calc. enrichment
    if (SEP_TPS == 0){
        dir_file <- Sys.glob(file.path(IN_OUTPUT, "*geneList*.csv"))
    } else {
        dir_file <- Sys.glob(file.path(IN_OUTPUT, "*geneName*.csv"))
    }

    if (! length(dir_file)>1){ # check to make sure only one geneList file per cluster
        geneNumClust=(sub(".*_([^.]+)\\.csv.*", "\\1", dir_file))
        OUTPUT_DIR <-paste(dirname(dir_file),'',sep="/")
        
	# Generating ranked gene list for gene set enrichment analysis (GSEA)
        f <- read.csv(dir_file, sep=",", header = TRUE)
        geneList<-f[,2]
        names(geneList)<- as.character(f[,1])
        geneList<-sort(geneList, decreasing=TRUE)
        
        # Reading the gene file for each cluster    
        a <- dput(as.character(f[1][,1]))
        
        # Removing ballgown ID and/or transcript ID if present at end of gene name
        if (grepl("_", a[1])){
            cat("Removing values after and including '_' at end of gene/transcript name\n")
            xs <- gsub("\\-R[A-Z][_|\\>].*","",a)
        } else if (grepl("-R[A-Z]\\>", a[1])) {
            cat("Removing transcript ID at end of gene name\n")
            xs <- gsub("\\-R[A-Z]\\>.*","",a)
        } else{
            xs <- unique(a)
        }
        
        # Converting gene names to ensemble and entrezids
        gene = bitr(geneID=xs, fromType="SYMBOL", toType=c("ENSEMBL","ENTREZID"), OrgDb=ORG_DB)

        # Creating ENTREZID matrix of gene trajectories for pathview.multi
        gene_n <- gene[,-1]
        rownames(gene_n) <- gene[,1]
        # Moving gene symbol to row names
        df_trajs<-f[,-1]
        rownames(df_trajs)<-f[,1]
        # Fetching gene ENTREZID and creating dataframe with ENTREZID as rowname
        df_ens <- merge(df_trajs, gene_n, by=0, all=TRUE)
        df_ens_n <- na.omit(df_ens)
        df_ens_n2 <- df_ens_n[,0:(ncol(df_ens_n)-2)]
        rownames(df_ens_n2) <- df_ens_n[,ncol(df_ens_n)]
        # Dropping gene symbol
        drops <- c("Row.names")
        df_trajs_ens <- df_ens_n2[ , !(names(df_ens_n2) %in% drops)]

        # Writing an intermediate file with gene names
        type_enrichment <- c("BP","MF","CC") # Biological Process (BP), Molecular Function (MF), Cellular Component (CC)

        # Plotting GO and pathway enrichment per type (BP, MF, CC)
        pathway_found=0 # flip to 1 when the most enriched pathway is identified
        for(types in type_enrichment){
            print(paste("Processing", geneNumClust, sep=" "))
            print(paste("Determining",types,"enrichment", sep=" "))
            ego <- assess_enrichment(gene$ENSEMBL, types)
            #ego3 <- assess_GSEA(geneList,type)
            
            # Iterating through Biological Processes, Molecular Function, and Cellular Component to find and plot top enriched pathway
            if(types=="BP"){ 
                for(enr_term in ego[,2]){
                    print(paste("enriched term", enr_term, sep=" "))
                    id_pathway <- plotPathway(df_trajs_ens, EXPERIMENT_NAME, enr_term, OUTPUT_DIR, geneNumClust)
                    if(!id_pathway==TRUE){ # if no pathway was identified (if id_pathway is empty)
                        next # search for the next most enriched term
                    } else{ # if enriched pathway identified, stop looking and exit out of loop
                        print(paste("Pathway found", id_pathway, sep=" "))
                        pathway_found=1
                        break
                    }
                }
            } else if (types=="MF" && pathway_found==0) { # if no enriched pathway was identified for BP
                for(enr_term in ego[,2]){
                    id_pathway <- plotPathway(gene$ENTREZID, EXPERIMENT_NAME, enr_term, OUTPUT_DIR, geneNumClust)
                    if(!id_pathway==TRUE){ # if no pathway was identified (if id_pathway is empty)
                        next # search for the next most enriched term
                    } else{ # if enriched pathway identified, stop looking and exit out of loop
                        print(paste("Pathway found", id_pathway, sep=" "))
                        pathway_found=1
                        break
                    }
                }   
            } else if (types=="CC" && pathway_found==0) { # if no enriched pathway was identified for BP or MF 
                for(enr_term in ego[,2]){
                    id_pathway <- plotPathway(gene$ENTREZID, EXPERIMENT_NAME, enr_term, OUTPUT_DIR, geneNumClust)
                    if(!id_pathway==TRUE){ # if no pathway was identified (if id_pathway is empty)
                        next # search for the next most enriched term
                    } else{ # if enriched pathway identified, stop looking and exit out of loop
                        print(paste("Pathway found", id_pathway, sep=" "))
                        pathway_found=1
                        break
                    }
                }
            } else if (pathway_found==1) {
                print("Pathway already enrichment found.")
            } else{
                print("No pathway enrichement found.")
            }
            print(paste("Plotting",types,"enrichment for", EXPERIMENT_NAME,sep=" "))
            generate_plots(ego, geneList, EXPERIMENT_NAME, types, OUTPUT_DIR, geneNumClust)
        }
        
    }
}

if(!interactive()) {
    main()
}
