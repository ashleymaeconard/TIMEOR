---
title: "Web Server"
output:
  html_document:
    css: github.css
---

Quick Start: 3 Steps
=======
 TIMEOR accepts 2 input types: (1) raw .fastq files and SraRunTable [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/real_data_subset/timeor/data/SraRunTable.csv) or a (2) RNA-seq time-series read count matrix [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/countMatrix.csv) and metadata file [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/metadata.csv).


1. Visit https://timeor.brown.edu.
2. For (1) in 'Example Data' (side-bar) under 'Load raw data' click the 'SraRunTable & .fastq files' button. This will guide you through the 'Process Raw Data' tab demo. **Follow pop-ups and fill in grey boxes**. See [Run TIMEOR](#run-timeor-from-raw-data-starting-from-fastq-time-series-rna-seq) below for walk-through.
3. Next, for (2) in 'Example Data' (side-bar) under 'Load count matrix' click the 'Metadata & read count file' button. This will guide you through the rest of the full method demo. **Follow pop-ups and fill in grey boxes**. See [Run TIMEOR](#run-timeor-using-simulated-data-starting-from-read-count-matrix) below for full application walk-through.


Website
=======

### Computational Biology Core at Brown Univeristy and DRSC/TRiP Functional Genomics Resources at Harvard Medical School Partnership
TIMEOR is available online at https://timeor.brown.edu.

Run TIMEOR
===================

### Two ways to input data:

1.   Import **SraRunTable from GEO**\* where TIMEOR will process raw data
    through retrieving .fastq files, quality control, alignment, and
    read count matrix creation. Read [this section](#run-timeor-from-raw-data-starting-from-fastq-time-series-rna-seq) below.

2.   Import **metadata file\*\* and count matrix \*\*\*** (skipping raw
    data retrieval, quality control, alignment, and read count matrix
    creation) and proceeding straight to normalization and correction. 
    Read [this section](#run-timeor-using-simulated-data-starting-from-read-count-matrix) below.

Then simply follow the prompts. Fill out the **grey** boxes to begin
interacting with each stage and tab. 

#### Input file types:

  \* **SraRunTable from GEO** follow instructions in TIMEOR first tab
    (“Process Raw Data”)

   \*\* **metadata file** requires at least these columns.
    -   *ID, condition, time, batch*
        -   *ID*: a unique identifier (ID) for the user
            (e.g. case\_1min\_rep1)
        -   *condition*: one word description (e.g. case, control)
        -   *time*: numerical values e.g. (0, 20, 40)
        -   *batch*: string description of batch (e.g. b1, b2, b3)

  \*\*\* **count matrix** : rows should be unique gene identifiers
    (e.g. Flybase, Ensembl or Entrez IDs) and columns should be the IDs
    from metadata file.

## Run TIMEOR from Raw Data: Starting from .fastq Time-Series RNA-seq
This tutorial uses a subset of real data used in the TIMEOR publication to
take the user through TIMEOR's "Process Raw Data" tab. You will first see this pop-up. Please read.
There are 4 steps.

<p>
     
</p>
<center>
<img src="T000.png" style="width:95.0%" />
</center>
<p>
     
</p>

1.  In the far-left side-bar click on “Example Data” and then
    under “Load raw data” click on "SraRunTable & raw .fastq files".

<p>
     
</p>
<center>
<img src="T00.png" style="width:95.0%" />
</center>
<p>
     
</p>

2.  Follow the pop-up prompt to explore the default settings to questions 1-6 to set the adaptive 
    default parameters, and then click the "Run" button to begin retrieving the raw data (SRR8843738
    and SRR8843750), performing quality control, and aligning the reads using HISAT2 and Bowtie2.

<p>
     
</p>
<center>
<img src="T01.png" style="width:95.0%" />
</center>
<p>
     
</p>



<p>
     
</p>
<center>
<img src="T02.png" style="width:95.0%" />
</center>
<p>
     
</p>

3.  Once the data have been retrieved and quality control has finished. You can view a summary under
    the "Quality Control" panel on the right. Interactive results can be downloaded and viewed from MultiQC.

<p>
     
</p>
<center>
<img src="T03.png" style="width:95.0%" />
</center>
<p>
     
</p>
 
<p>
     
</p>
<center>
<img src="T04.png" style="width:95.0%" />
</center>
<p>
     
</p>   

4.  You can explore the alignment results between both methods in the "Alignment Quality" panel. Note that
    HISAT2 is splice-site aware). You can choose the method and then click "Generate count matrix" to have 
    TIMEOR generate the read count matrix for the next tab "Load Count Matrix". You will see a green check mark when finished. Follow pop-up. You have completed this section of the tutorial for loading raw time-series RNA-seq data.

<p>
     
</p>
<center>
<img src="T05.png" style="width:95.0%" />
</center>
<p>
     
</p> 

<p>
     
</p>
<center>
<img src="T050.png" style="width:95.0%" />
</center>
<p>
     
</p> 

## Run TIMEOR Using Simulated Data: Starting from Read Count Matrix

This tutorial uses simulaated data and takes the user through TIMEOR’s full
functionality beginning from a read count matrix (genes x sample/time). **NOTE**: figures with two panels are the same page,
just split. There are 20 steps.

**The user can begin this tutorial before *or* after following ["Run TIMEOR from Raw Data: Starting from .fastq Time-Series RNA-seq"](#run-timeor-from-raw-data-starting-from-fastq-time-series-rna-seq).**

1.  In the far-left side-bar click on “Example Data” and then
    under “Load simulated data” click on "Metadata & read count file".

<p>
     
</p>
<center>
<img src="T1.png" style="width:95.0%" />
</center>
<p>
     
</p>

2.  Follow the pop-up to explore results on the Pre-processing
    tabs "Process Raw Data" and "Process Count Matrix".

<p>
 
</p>
<center>
<img src="T2.png" style="width:95.0%" />
</center>
<p>
 
</p>
<p>
 
</p>
<center>
<img src="T3.png" style="width:95.0%" />
</center>
<p>
 
</p>

3.  On the Normalize and Correct Data tab, choose from normalization and
    correction methods and click “Run” to view result. You have completed pre-processing. Follow pop-up.

<p>
 
</p>
<center>
<img src="T4.png" style="width:95.0%" />
</center>
<p>
 
</p>

<p>
 
</p>
<center>
<img src="T40.png" style="width:95.0%" />
</center>
<p>
 
</p>

4.  Proceed to Primary Analysis and click “Run”.

5.  At the bottom right you will see a pop-up to click “Render
    Venn Diagram” in the top right to compare differential expression
    results between three methods (ImpulseDE2, Next maSigPro, and
    DESeq2) and choose which method results to proceed with. 
    You will then see a pop-up saying that you have completed Primary Analysis. 
    Feel free to move on with TIMEOR's default parameters, **or explore Primary Analysis options (see next step)**.

<p>
 
</p>
<center>
<img src="T50.png" style="width:95.0%" />
</center>
<p>
 
</p>

<p>
 
</p>
<center>
<img src="T51.png" style="width:95.0%" />
</center>
<p>
 
</p>

6.  See 3-way Venn diagram in white box. You can compare a previous study list of genes with the three differential expression
    results by first downloading [prev\_study.txt](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/prev_study.txt). Then upload this file using the Browse button. See 4-way Venn diagram in blue box. 

<p>
 
</p>
<center>
<img src="T52.png" style="width:95.0%" />
</center>
<p>
 
</p>

7.  Examine differential expression method results in the bottom row.
    Toggle under “Display Desired Differential Expression Method
    Results” between ImpulseDE2, Next maSigPro, and DESeq2 on the left,
    and the interactive clustermap with automated clustering will
    display the differentially expressed gene trajectories for the
    chosen method. You can then toggle under “Cluster Gene Expression Trajectories” to choose the
    number of clusters desired. 
    For this demonstration, we chose "ImpulseDE2" differentially expressed gene output. 
    We also chose "automatic" clustering of gene trajectories. On new data the user can choose these two parameters.
    **NOTE** ImpulseDE2 is chosen because
    it has the largest differential expressed gene overlap with the
    previous study and other methods.

<p>
 
</p>
<center>
<img src="T5.png" style="width:95.0%" />
</center>
<p>
 
</p>
8. As said in pop-up, proceed to Secondary Analysis tab in side-bar.

9. Clusters are labeled in ascending order from 1 for top-most cluster. Under Gene Expression Trajectory Clusters choose cluster 1, 2, or 3
    in the dropdown. On the right under “Chosen Cluster Gene Set” you
    will see the genes in that cluster. **Genes are the same color as
    the gene trajectory cluster to which they belong**.

10. Once you have chosen which genes set to test for enrichment, click
    the “Analyze” toggle to “ON”.

11.  Wait to view any enriched gene ontology (GO) terms (Molecular
    Function, Biological Process, or Cellular Component), pathway,
    network, and/or motif analysis. **NOTE** you may download the
    interactive motif results for viewing.

12.  Toggle the “Analyze” button to “OFF” to choose another gene set, and
    repeat steps 9-12.

<p>
 
</p>
<center>
<img src="T6.png" style="width:95.0%" />
</center>
<p>
 
</p>

13. Proceed to the "Factor Binding" tab to view the *observed* (that is, differentially expressed genes in data) and top
    predicted transcription factors in each gene cluster (under
    “Observed and Top Predicted Transcription Factors”). At least 40% of the transcription factor prediction methods must agree on their top predicted transcription factors, otherwise cell is left blank. Row names are gene expression trajectory clusters. **NOTE** for this demonstration the threshold is reduced.
    
<p>
 
</p>
<center>
<img src="T70.png" style="width:95.0%" />
</center>
<p>
 
</p>

14.  In that same table on the far right you will see ENCODE IDs indicating
    published ChIP-seq data for the predicted transcription factors. For this tutorial, can either see an example provided with "pho". 
    You may also either download these read-depth normalized .bigWig files here ([ENCFF467OWR](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/ENCFF467OWR.bigWig), [ENCFF609FCZ](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/ENCFF609FCZ.bigWig), [ENCFF346CDA](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/ENCFF346CDA.bigWig)) or follow the prompts (step 16) in the grey box under “Upload
    .bigWig Files”. Any .bigWig files from protein-DNA data are accepted. 
    
15. If you are interested, click on the “+” under
    “See each method's predicted transcription factors:” to
    see the ranked lists of transcription factors and motifs by method.
    Blanks indicate an enriched motif is not assigned to a transcription factor region (to see motifs click 'Download interactive cluster motif result'). Search for a method (e.g. transfac in blue box), enrichment score, etc. Row names are top 1 - 4 transcription factors. 

<p>
 
</p>
<center>
<img src="T7.png" style="width:95.0%" />
</center>
<p>
 
</p>

16.  Under “Average Profiles Across Each Gene Expression Trajectory
    Cluster”, in the second box type “stat92e”, upload
    ENCFF467OWR.bigWig, and click “Go”. In the third box type
    “CG7786”, upload ENCFF346CDA.bigWig, and click “Go”. You will see 3
    average profile distributions in each plot and a heatmap with 3 partitions (one for gene trajectory each cluster). 
    **NOTE** each color of the distribution matches the corresponding gene trajectory cluster.

<p>
 
</p>
<center>
<img src="T8.png" style="width:85.0%" />
</center>
<p>
 
</p>

17.  Proceed to the last tab to view the temporal relations between
    transcription factors. On the first row you are reminded of the
    observed and predicted transcription factors to bind each gene
    cluster. 
  <p>
 
</p>
<center>
<img src="T9.png" style="width:95.0%" />
</center>
<p>
 
</p> 

18. Specifically, zooming into the second row to the left you will see “Transcription
    Factor Network" of both obiserved and predicted transcription
    factors. On the right (”Temporal Relations Between Observed and
    Predicted Transcription Factors") you will see a table highlighting
    the temporal relations between transcription factors. **TIMEOR identified the transcription factor (TF) gene regulatory network (GRN).**
    These temporal relationships are identified and represented in the legend
    (far-right). <span style="color:#3F88DE"> **Blue arrow/highlight:** predicted TF to observed TF, experimentally determined interaction</span>. <span style="color:#D6678D"> **Pink arrow/highlight:** observed TF to observed TF, experimentally determined interaction</span>. <span style="color:#F7C144"> **Yellow arrow/highlight:** observed TF to observed TF, predicted interaction</span>. <span style="color:#5B8179"> **Green arrow/highlight:** predicted TF to observed TF, predicted interaction</span>. Network displayed in table format in app to enhance flexibility of GRN visualization.

<p>
 
</p>
<center>
<img src="T90.png" style="width:95.0%" />
</center>
<p>
 
</p>

19.  On the third row (“Network Customization: move and add desired genes
    to describe temporal relation”) the user can use this information to
    create a customized network to temporally relate all transcription
    factors and other genes. Do so by clicking “Search” and then
    “Multiple proteins”.
    
<p>
 
</p>
<center>
<img src="T10.png" style="width:95.0%" />
</center>
<p>
 
</p>

20. Your results folder can be downloaded on the far-left side-bar under "Download Results Folder". **NOTE** The original simulated data and results can be downloaded [here](https://github.com/ashleymaeconard/TIMEOR/tree/master/demos/simulated_data). 

<p>
 
</p>
<center>
<img src="T11.png" style="width:95.0%" />
</center>
<p>
 
</p>

Details
=======

### Real Data Subset in Tutorial

The original temporal RNA-seq data analyzed in our paper comes from Zirin et al., 2019). In this tutorial SRR8843750 and SRR8843738 are analyzed to demonstrate the "Process Raw Data" tab in which raw RNA-seq data are retrieved, quality checked, aligned (with HISAT2 and Bowtie2), and converted to a read count matrix. The real data subset folder (which TIMEOR automatically generates) can be downloaded [here](https://github.com/ashleymaeconard/TIMEOR/tree/master/demos/real_data_subset).

### Simulated Data in Tutorial

The original simulated data folder can be downloaded [here](https://github.com/ashleymaeconard/TIMEOR/tree/master/demos/simulated_data).

#### Secondary Analysis: Factor Binding

To get the top 4 TFs a 25% concensus threshold was used, with a
normalized enrichement score threshold of 3.

Command used:
`Rscript get_top_tfs.r /PATH/TO/simulated_results/ dme 3 4 25 /PATH/TO/TIMEOR/`

The following bigWig files were collected:

-   ENCFF467OWR (read-depth normalized signal between both replicates)
    within dataset ENCSR240ADR for Stat92E

-   ENCFF609FCZ (read-depth normalized signal between both replicates)
    within dataset ENCSR681YMA for pho

-   ENCFF346CDA (read-depth normalized signal between all three replicates)
    within dataset ENCSR776AVR for CG7786

### Real Data in Publication

The results presented in TIMEOR's publication can be downloaded in TIMEOR's automatically generated folders [here](https://github.com/ashleymaeconard/TIMEOR/tree/master/demos/real_data_full).
