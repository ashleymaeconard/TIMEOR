#!/bin/bash
# create_avg_prof.sh
# Ashley Mae Conard
# Last Mod: Dec 29, 2019
# Resources: https://www.biostars.org/p/42844/
#            http://seqanswers.com/forums/archive/index.php/t-45817.html (remove chr from file and remove header)

# Checking to make sure input is correct
if [ $# -ne 5 ]; then
	echo $0: "Usage: ./create_avg_profs.sh
            1) /PATH/TO/ENCFF.bigWig 
            2) TF_NAME 
            3) /PATH/TO/RESULTS/ANALYSIS/NAME_results/ 
            4) /PATH/TO/NAME_results/factor_binding/
            5) /PATH/TO/CLUST_NUM_TO_COLOR_FILE.txt\n"
	exit 1
fi

# Inputting arguments
INPUT_BIGWIG=$1
TF=$2
CLUSTERS_FOLDER=$3
OUTPUT_DIR=$4
CLUST_NUM_2_COLOR=$5

# Importing colors as list
IFS=$'\r\n' GLOBIGNORE='*' command eval  'LIST_COLORS=($(cat ${CLUST_NUM_2_COLOR}))'

# Creating list of all .bed files of genes for each cluster
LIST_BEDS=()
for i in ${CLUSTERS_FOLDER}/clusters/
do 
    for geneList in ${CLUSTERS_FOLDER}/clusters/*/MEME/*geneList*.bed # for each gene cluster
    do
        LIST_BEDS+=(${geneList})
    done
done

# Computing matrix
computeMatrix scale-regions -S ${INPUT_BIGWIG} -R ${LIST_BEDS[@]} --binSize 250 --beforeRegionStartLength 1000 --afterRegionStartLength 1000 --regionBodyLength 5000 -o ${OUTPUT_DIR}/matrix.genes.clusters.mat.gz --skipZeros --smartLabels --sortRegions descend

# Plotting heatmap
plotHeatmap -m ${OUTPUT_DIR}/matrix.genes.clusters.mat.gz -out ${OUTPUT_DIR}/heatmap_genes.clusters.${TF}.png --heatmapHeight 25 --heatmapWidth 15 --labelRotation 45 --missingDataColor red --whatToShow="heatmap and colorbar"

# Plotting average profile - with standard error (used colorspace::rainbow_hcl(15) to determine colors used in heatmaply)
plotProfile -m ${OUTPUT_DIR}/matrix.genes.clusters.mat.gz -out ${OUTPUT_DIR}/avg_profile.${TF}.png --plotType=se --labelRotation 45 --plotHeight 15 --plotWidth 15 --colors ${LIST_COLORS[@]} -T="${TF} Average Profile"

echo "Saved heatmap.clusters.${TF}.png and average avg_profile.${TF}.png in ${OUTPUT_DIR}"
rm ${OUTPUT_DIR}"/matrix.genes.clusters.mat.gz"
