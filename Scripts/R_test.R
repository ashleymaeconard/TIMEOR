###################################
# Test R script to load packages
###################################

# Test loads in R 
install.packages("pacman")
library(pacman)
p_load(c("openssl", "httr", "rvest", "xml2")) # works!!
p_load(c("UpSetR, corrplot","autoplotly","BiocManager","broom","cluster","d3heatmap","data.table","devtools","doMC","doRNG","dplyr","DT",
                   "europepmc","factoextra","farver","fastmatch","fpc","ggforce","ggfortify","ggplot2","ggplotify","ggridges","GlobalOptions","graphlayouts","gridGraphics","heatmaply","markdown",
                   "plotly","png","polyclip","promises","reticulate","rmarkdown","rvcheck","shiny","shinyalert","shinycssloaders","shinydashboard","shinydashboardPlus","shinyjs","shinyLP",
                   "shinyWidgets","stringr","tibble","tidygraph","tidyr","tidyverse","triebeard","tweenr","urltools","vegan","zoo"))  # works! 

library("clusterProfiler")
library("Harman")
library("ImpulseDE2") # failed 
library("maSigPro")
library("ENCODExplorer")
library("enrichplot")
library("fgsea")
library("IRanges")
library("pathview")  # failed - glmnet
library("RcisTarget") # failed - glmnet 
library("STRINGdb") # failed - glmnet 
library("topGO") # failed - glmnet 
library("future") 
library("heatmaply")
library("org.Dm.eg.db") # failed - glmnet 
library("org.Hs.eg.db")
library("org.Mm.eg.db") # failed - glmnet 
library("DESeq2")
library("BiocParallel")
library("biomaRt") 
library("DOSE")
