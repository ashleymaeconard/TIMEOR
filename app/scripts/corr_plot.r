# corr_plot.r
# Ashley Mae Conard
# Last Mod. 11/11/2019
# Purpose: Creates plotly correlation plot between samples using pearson

# Input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
    cat("Passed in:", args,"\n")
} else{
    stop("Usage: Rscript corr_plot.r /FULL/PATH/TO/geneXsample.csv /FULL/PATH/TO/OUTPUT_DIR/")
}
GENExSAMP <- args[1]
OUTPUT_DIR <- args[2]

# Make directory if does not exist
dir.create(file.path(OUTPUT_DIR), recursive=TRUE)

# Importing libraries
library(heatmaply)
library(ggplot2)
library(htmlwidgets)

# Reformatting 
geneXsamp <- read.csv(GENExSAMP, sep=",", header = TRUE)
rownames(geneXsamp)<-geneXsamp[,1]
geneXsamp[,1] <- NULL
geneXsamp <- geneXsamp[which(rowSums(geneXsamp) > 0),]

# Plot and save correlation plot
output_filename <- paste("corr_samples.html",sep="")
outdir_withfile <- (paste(OUTPUT_DIR, output_filename, sep="/"))

p <- heatmaply(cor(geneXsamp, method = c("spearman")), 
          xlab = "Experiments", ylab = "Experiments", 
          main = paste("Spearman Correlation Between Experiments"),
          margins = c(40, 40),
          limits = c(-1,1), 
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"))
htmlwidgets::saveWidget(p, outdir_withfile)
