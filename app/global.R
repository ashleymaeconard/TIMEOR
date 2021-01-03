# global.R
# Ashley Mae Conard
# Last Mod. 6 June 2020

library(shiny, shinydashboard, DT)
library(ggplot2)
library(plotly)
library(devtools)
library(htmlwidgets)
library(dplyr)
library(Harman)
library(BiocManager)

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










