#######################################
# R script for Ashley's Dockerfile
#######################################

# Install necessary packages in R from CRAN and Bioconductor 

# This works
install.packages(c("openssl", "httr", "rvest", "xml2"), dependencies=TRUE, repos='http://cran.rstudio.com/')

# This works 
install.packages("UpSetR", dependencies=TRUE, repos='http://cran.rstudio.com/')

# This works 
install.packages(c("corrplot","autoplotly","BiocManager","broom", "cluster"), dependencies=TRUE, repos='http://cran.rstudio.com/')

#########################################################################################
# d3heatmap doesn't work because it was removed from CRAN - we use heatmaply instead
#########################################################################################
#install.packages("d3heatmap", dependencies=TRUE, repos='http://cran.rstudio.com/')

# This works
install.packages(c("data.table","devtools","doMC","doRNG","dplyr","DT"), dependencies=TRUE, repos='http://cran.rstudio.com/')

# This works 
install.packages(c("europepmc","factoextra","farver","fastmatch","fpc","ggforce","ggfortify","ggplot2"), dependencies=TRUE, repos='http://cran.rstudio.com/')
# This works 
install.packages(c("ggplotify","ggridges","GlobalOptions","graphlayouts","gridGraphics","markdown"), dependencies=TRUE, repos='http://cran.rstudio.com/')
# This works 
install.packages(c("plotly","png","polyclip","promises","reticulate","rmarkdown","rvcheck","shiny","shinyalert"), dependencies=TRUE, repos='http://cran.rstudio.com/')
# This works 
install.packages(c("shinycssloaders","shinydashboard","shinydashboardPlus","shinyjs","shinyLP"), dependencies=TRUE, repos='http://cran.rstudio.com/')
# This works 
install.packages(c("shinyWidgets","stringr","tibble","tidygraph","tidyr","tidyverse","triebeard","tweenr"), 
                 dependencies=TRUE, repos='http://cran.rstudio.com/')
# This works
install.packages(c("urltools","vegan","zoo"), dependencies=TRUE, repos='http://cran.rstudio.com/')
# This works 
BiocManager::install("clusterProfiler", dependencies=TRUE)
# This works 
BiocManager::install("Harman", dependencies=TRUE)
# Don't think you need this 
#install.packages("glmnet", dependencies=TRUE, repos='http://cran.rstudio.com/')
# This works 
BiocManager::install("ImpulseDE2", dependencies=TRUE) 
# This works
BiocManager::install("maSigPro", dependencies=TRUE)
# This works
BiocManager::install("ENCODExplorer", dependencies=TRUE)
# This works 
BiocManager::install("enrichplot", dependencies=TRUE)
# All code below works 
BiocManager::install("fgsea", dependencies=TRUE)
BiocManager::install("IRanges", dependencies=TRUE)
BiocManager::install("pathview", dependencies=TRUE) 
BiocManager::install("RcisTarget", dependencies=TRUE) 
BiocManager::install("STRINGdb", dependencies=TRUE) 
BiocManager::install("topGO", dependencies=TRUE) 
install.packages("future", dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages("heatmaply", dependencies=TRUE, repos='http://cran.rstudio.com/')
BiocManager::install("org.Dm.eg.db", dependencies=TRUE) 
BiocManager::install("org.Hs.eg.db", dependencies=TRUE)
BiocManager::install("org.Mm.eg.db", dependencies=TRUE) 
BiocManager::install("DESeq2", dependencies=TRUE)
BiocManager::install("BiocParallel", dependencies=TRUE)
BiocManager::install("biomaRt", dependencies=TRUE) 
BiocManager::install("DOSE", dependencies=TRUE)