# Ashley Conard
# clustermap.r
# Last Mod. 12/12/2019
# Resources: https://declara.com/content/lgAyPkg4
#           If issues installing packages: # install.packages("vegan", repos = "https://cran.rstudio.com", dependencies = TRUE)

# Input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Pass in 7 arguments: Rscript clustermap.r \n 1)/FULL/PATH/TO/RESULTS_DIR/ \n 2)/FULL/PATH/TO/HEATMAP_INPUT_FILE \n 3)DESIRED_OUTPUT_NAME (e.g. insulin_test) \n 4)CLOSE_TIMEPOINT (1: close timepoint, 0: distant timepoint) \n 5)USER_CHOOSE_CLUST_NUMBER (if 0 TIMEOR finds optimal, else type desired number of clusters) \n 6)DIST_METHOD (e.g. euclidean) \n 7)CLUSTER_METHOD (e.g. ward.D2)", call.=FALSE)
} else if (length(args) == 7) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 7 arguments: Rscript clustermap.r \n 1)/FULL/PATH/TO/RESULTS_DIR/ \n 2)/FULL/PATH/TO/HEATMAP_INPUT_FILE \n 3)DESIRED_OUTPUT_NAME (e.g. insulin_test) \n 4)CLOSE_TIMEPOINT (1: close timepoint, 0: distant timepoint) \n 5)USER_CHOOSE_CLUST_NUMBER (if 0 TIMEOR finds optimal, else type desired number of clusters) \n 6)DIST_METHOD (e.g. euclidean) \n 7)CLUSTER_METHOD (e.g. ward.D2)")
}

# Parameters
RESULTS_DIR <- args[1] # e.g. ../timeor/results/
HEATMAP_INPUT_FILE <- args[2] # e.g. ../timeor/results/primary/DE/insulin_stim_zs.csv
LIST_EXP <- args[3] # insulin
CT <- as.numeric(args[4]) # 1 if close timepoint, 0 otherwise
USER_CHOOSE_CLUST_NUMBER <- as.numeric(args[5]) # change this to some non-negative integer
DIST_METHOD = args[6] #"euclidean","maximum", "canberra", "manhattan", "binary", "minkowski"
HCLUSTER_METHOD = args[7] #"ward.D", "ward.D2", "single", "complete", "average", mcquitty, "median", "centroid"

# Arguments that could be parameters if desired
NUM_CLUSTERS_TO_TEST = 15

# Libraries
# library(data.table)
# source("geneID_converter.r") # find in scripts folder!
# library(plotly)
# library(htmlwidgets)
# library(heatmaply)


# Libraries
library(data.table)
library("org.Dm.eg.db")
source("./scripts/geneID_converter.r") # find in scripts folder!
library(plotly)
library(htmlwidgets)
library(heatmaply)
library(ggplot2)
library(gplots)
library(devtools)
library(dplyr)
library(knitr)
library(gclus)
library(fpc)
library(dendextend)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)}
if (!require("vegan")) {
  install.packages("vegan")
  library(vegan)}

# Parameters
#produceClusterMap <- function(RESULTS_DIR, HEATMAP_INPUT_FILE, LIST_EXP, CT, USER_CHOOSE_CLUST_NUMBER, DIST_METHOD, HCLUSTER_METHOD){

# Arguments that could be parameters if desired
NUM_CLUSTERS_TO_TEST = 15
USER_CHOOSE_CLUST_NUMBER = as.numeric(USER_CHOOSE_CLUST_NUMBER)

# Get mode between list of numbers (used to find optimal number clusters from 3 clustering methods, not including Elbow)
getmode <- function(v) {
uniqv <- unique(v)
uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Add the gene_name (i.e. gene symbol) to the input dataframe
convertInput <- function(dir, gene_ID_db_name){
# Read csv
d <- read.csv(dir, sep=",", header = TRUE)
# If gene_name (i.e. gene symbol) already exists, return dataframe with rows indexed by gene_name
if(any("gene_name"==colnames(d))){
  return(d <- read.csv(dir, sep=",", header = TRUE, row.names="gene_name"))

  # Otherwise, add gene_name to dataframe and return dataframe with rows indexed by gene_name
} else{
  ids <- rownames(d)
  ids.type  <- gene_ID_db_name
  gene_names <- as.vector(get.symbolIDsDm(ids,ids.type)) # from geneID_converter.r
  d1 <- cbind(gene_names, d)
  setDT(d1, keep.rownames = TRUE)[]
  colnames(d1)[1] = "gene_id"
  colnames(d1)[2] = "gene_name"

  # If some gene names are NA, insert the FlybaseID
  if(length(d1$gene_name[is.na(d1)]) != 0){ # check if there are any NA in gene_name
    # Make sure gene_id and gene_name columns are of type chr 
    d1 %>% mutate_if(is.factor, as.character) -> d1
    d1 %>% mutate_if(is.list, as.character) -> d1
    d1$gene_name[is.na(d1$gene_name)] <-as.character(d1$gene_id[is.na(d1$gene_name)])
  }
  return(d1)
}
}

# Plot likely number of clusters and choose mode between 3 clustering methods
find_plot_num_clusters <- function(a, experiment, NUM_DESIRED_CLUSTERS, CHOOSE_NUM_CLUSTERS) {
# If CHOOSE_NUM_CLUSTERS == 1, return number of clusters and change to 0 for next time, 
# IF CHOOSE_NUM_CLUSTERS == 0, use NUM_DESIRED_CLUSTERS as consensus and save the actual clusters (number based on CHOOSE_NUM_CLUSTERS)

# Make the gene_names the rownames
#a <- data.frame(a[,-1], row.names=a[,1]) 
print(head(a))
# User defines the number of clusters
if (USER_CHOOSE_CLUST_NUMBER >= 1 & CHOOSE_NUM_CLUSTERS == 1){
  cat("User set the number of clusters to", USER_CHOOSE_CLUST_NUMBER,"\n")
  optim_num_clust = USER_CHOOSE_CLUST_NUMBER
  return(optim_num_clust)

  # Determine the number of clusters and plot clustering methods' results
} else if (USER_CHOOSE_CLUST_NUMBER==0 & CHOOSE_NUM_CLUSTERS==1) { 

  # Create output file name to display clustering methods
  if (CT){
    ending='tc'
  } else{
    ending='l2FC'
  }
  output_filename <- paste(experiment,'_',ending,'_3_methods.pdf',sep='')
  pdf(paste(RESULTS_DIR, paste(experiment,"_results", sep=""), 'clusters',output_filename ,sep="/")) 

  # Plot Elbow
  tot_wss = c()
  for(i in 1:NUM_DESIRED_CLUSTERS)
  {
    cl <- kmeans(a,centers = i)
    tot_wss[i] <- cl$tot.withinss
  }           

  par(mfrow=c(2,2))
  plot(x=1:NUM_DESIRED_CLUSTERS,
      y=tot_wss,
      type = "b",
      xlab = "Number of clusters",
      ylab = "Within groups sum of squares")  
  title(paste("Elbow number clusters: ", experiment, sep=' '))  

  # Choosing number clusters from mode of 3 methods: partition, silouhette, calinsky criterion  

  # Partition around medoids
  pamk.best <- pamk(a) 
  plot(pam(a, pamk.best$nc))
  cat("Optimal partition around medoids gives ",pamk.best$nc," clusters.\n")

  # Silouhette
  asw <- numeric(NUM_DESIRED_CLUSTERS) 
  for(k in 3:NUM_DESIRED_CLUSTERS) {
    asw[[k]] <- pam (a, k) $ silinfo $ avg.width 
  }    
  k.best <- which.max(asw) 
  cat("Optimal cluster number with Silouhette gives",k.best," clusters.\n")

  # Calinsky criterion
  fit <- cascadeKM(scale(a, center = TRUE, scale = TRUE),1, NUM_DESIRED_CLUSTERS, iter = 500)
  plot(fit, sortg=TRUE, grpts.plot = TRUE)
  calinski.best <- as.numeric(which.max(fit$results[2,]))
  cat("Optimal cluster number with Calinski criterion gives", calinski.best," clusters.\n")
  dev.off()

  # Return optimal number of clusters from 3 clustering methods
  v <- c(pamk.best$nc, k.best, calinski.best)
  optim_num_clust <- getmode(v)
  return(optim_num_clust) # return the number of consensus # optimal clusters 

  # Return the gene groups (clusters) themselves
} else if (CHOOSE_NUM_CLUSTERS==0) 

  # Get distance between all genes in a
  d = dist(a)

# Cluster with gene and cluster level parameters
hr <- hclust(d, method = HCLUSTER_METHOD)

# Assign each gene to cluster. Specifically, mycl is a vector of integers which assigns each point (i.e. gene) to a cluster.
mycl <- cutree(hr, NUM_DESIRED_CLUSTERS)

# Partition all genes into their assigned cluster defined by mycl.
clusters <- lapply(unique(mycl), function(grp){
  a[which(mycl==grp),]
})

clust_details = list(clus=clusters,ordering=hr)

# Return clusters of genes
return(clust_details) 
}  

# Create or set output directory
if (!dir.exists(RESULTS_DIR)){
dir.create(RESULTS_DIR)
} else {
print("Directory already exists.")
}
cat("Output directory: ", RESULTS_DIR)

# Assigning database
gene_ID_database <- toTable(org.Dm.egFLYBASE)
gene_ID_database_name <- "flybase"

# Plotly's clustering method is the same as hcluster's method
hclust_heatmaply_method = HCLUSTER_METHOD

# Heatmap color scheme
color_scheme = ggplot2::scale_fill_gradient2(low = "midnightblue", high = "darkred") 

# Begin processing
print("\nCalculating # clusters and producing cluster plots and interactive heatmaps using log2FC.")
cat("List of experiments from input: ", LIST_EXP)

for(experiment in LIST_EXP){
cat("\n","\n", experiment, "\n")

# Create output dir for each experiment
outdir <- (paste(RESULTS_DIR, paste(experiment,"_results", sep=""),sep="/"))
out_subdir <- paste(outdir,'clusters',sep="/")

# Folder for /NAME_results/
if (!dir.exists(outdir)){
  dir.create(outdir)
  cat("Created:", outdir)
}

# Folder for /clusters/
if (!dir.exists(out_subdir)){
  dir.create(out_subdir)
  cat("Created:", out_subdir)
}

# Add gene_name (i.e. gene symbol) to input matrix
x <- convertInput(HEATMAP_INPUT_FILE, gene_ID_database_name)           
x_m <- data.frame(x[,-1])# remove Flybase ID column
#x_m <- data.frame(x[,-1], row.names=x[,1])

# Find optimal number of clusters from mode of 3 clustering methods
if(USER_CHOOSE_CLUST_NUMBER==0){
  CHOOSE_NUM_CLUST = 1 # call find_plot_num_clusters to get consensus on number of clusters 
  opti_num_clust <- find_plot_num_clusters(x_m, experiment, NUM_CLUSTERS_TO_TEST, CHOOSE_NUM_CLUST)
  cat("Consensus on optimal number cluster to use is: ", opti_num_clust,"\n")
} else{
  cat("User chose number of clusters to be: ",USER_CHOOSE_CLUST_NUMBER)
  opti_num_clust=USER_CHOOSE_CLUST_NUMBER
}
cat("opti", head(opti_num_clust))
# Use opti_num_clust to get actual clusters
CHOOSE_NUM_CLUST = 0 
cluster_details <- find_plot_num_clusters(x_m, experiment, opti_num_clust, CHOOSE_NUM_CLUST)
cat("\nTotal number of genes: ", dim(x_m),"\n")

# Saving clusters to output files
for(i in 1:opti_num_clust){ 

  # Create output dir for each experiment
  outdir_cluster <- (paste(RESULTS_DIR,paste(experiment,"_results", sep=""), 'clusters',i,sep="/")) 
  if (!dir.exists(outdir_cluster)){
    dir.create(outdir_cluster)
    cat("\nCreated:", outdir_cluster, opti_num_clust)
  }

  # if close timepoint
  if (CT){
    output_file_ending <- '_cluster_tc_geneList_'
  # else distant timepoint
  } else{
    output_file_ending <- '_cluster_l2FC_geneList_'
  }
  output_filename <- paste(experiment,output_file_ending,i,".csv",sep='')

  #  Save the divided input matrix into cluster specific folders.
  write.csv(cluster_details$clus[[i]], file = paste(RESULTS_DIR,paste(experiment,"_results", sep=""),'/clusters/', i, output_filename,sep='/'))
  # Print statement highlights to user that all genes have been accounted for in a cluster. 
  cat("\n","Cluster" , i, "has (dim of matrix) genes", dim(cluster_details$clus[[i]])) 
}

# Remove cluster subfolders that are left over
dirs <- list.dirs(path = out_subdir)
dirs <- dirs[2:length(dirs)]
for(d in dirs){
    if(!grepl("heatmaply", d, fixed = TRUE)){
        if(strtoi(basename(d), base = 0L) > strtoi(opti_num_clust, base = 0L)){
            unlink(d, recursive=TRUE)
        }
    }
}

# Plotting interactive clustermaps (heatmaply)
if (CT){
  ending='tc'
  title=paste("Clustermap of", experiment, "DEGs Across Time",sep=" ")
} else{
  ending='l2FC'
  title=paste("Log_2 Fold Change Clustermap of", experiment, "of DEGs Across Time",sep=" ")
}
output_filename <- paste(experiment,"_heatmaply_",ending,".html",sep="")
outdir_withfile <- (paste(RESULTS_DIR,paste(experiment,"_results", sep=""), '/clusters/', output_filename, sep="/"))

cat("\nGene distance methd: ", DIST_METHOD)
cat("\nCluster distance method: ", hclust_heatmaply_method)

# Getting colors of dendrogram to map to cluster folder labels (positive integers)
dend <- as.dendrogram(cluster_details$ordering)
dend <- color_branches(dend, opti_num_clust)
colors <- get_leaves_branches_col(dend)
colors_uniq <- rev(unique(colors))

# Create or replace cluster_color_to_number_map file
cluster_num_n_colors <- paste(dirname(outdir_cluster),'cluster_color_to_number_map.txt', sep="/")

if(file.exists(cluster_num_n_colors)){
    file.remove(cluster_num_n_colors)
}
write.table(colors_uniq, cluster_num_n_colors, append=T, col.names=F, row.names=F, quote=F, sep=",")

p <- heatmaply(x_m, Rowv=dend, dist_method = DIST_METHOD, hclust_method = hclust_heatmaply_method, xlab = "Stages", ylab = "Genes", 
              main = title, seriate="none",
              dendrogram = "row", k_row = opti_num_clust, margins = c(40, 130), 
              scale_fill_gradient_fun = color_scheme)
htmlwidgets::saveWidget(p, outdir_withfile)    
}
