<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/timeor_logo.png" height="150">

### Trajectory Inference and Mechanism Exploration with Omics in R

####  Author: Ashley Mae Conard

Motivation
==========

It's about time! [Click here](https://www.youtube.com/watch?v=bBodWvGPD6g&t=8s) for a quick video demonstration of the TIMEOR application, and [click here](https://www.youtube.com/watch?v=-BV9B3L1ymg) for a video guide through the tutorial below.

Analyzing time series differential gene expression and other multi-omics
data is computationally laborious and full of complex choices for the
various types of time series experiment analyses. The TIMEOR
web-application offers a solution as an interactive and adaptive R-Shiny
web interface to **support reproducibility and the best tools for a given
experimental design**. 

TIMEOR can take in either raw .fastq files or a raw
pre-computed count matrix (genes by samples), and after answering six
questions about your experiment, **performs all analysis from quality
control and differential gene expression to gene trajectory clustering
and transcription factor gene regulatory network construction**. TIMEOR
also suggests and allows users to integrate ChIP-seq data from ENCODE to
help vialidate predicted transcription factors binding to gene clusters.
Through a suite of published and novel methods, TIMEOR guides the user
through three stages, suggesting adaptive default methods and comparing
results for multiple normalization, alignment, and differential
expression (DE) methods. Those stages are **Pre-processing (generating
count matrix, normalizing and correcting data), Primary Analysis (DE and
gene trajectory clustersing) and Secondary Analysis (enrichment, factor
binding, and temporal relations)**. 

TIMEOR’s **interactive data visualizations and publication-ready figures 
streamline the process of time series data analysis and assist the user 
to design follow-up experiments**. At any point in the analysis, the user 
can scroll down and click **“Save your place”** to return to analysis later. 
TIMEOR is available for *Homo sapiens, Mus musculus,* and *Drosophila melanogaster*.
The web server is **completely free**, and is **accessible** at https://timeor.brown.edu in **a partnership** between the [Computational Biology Core at Brown University](https://cbc.brown.edu/) and [Harvard Medical School's DRSC TRiP Core Facility](https://fgr.hms.harvard.edu/). It is **also available** through both a Conda environment and Docker.

Quick Start: 3 Steps
=======
 TIMEOR accepts 2 input types: (1) raw .fastq files and SraRunTable [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/real_data_subset/timeor/data/SraRunTable.csv) or a (2) RNA-seq time-series read count matrix [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/countMatrix.csv) and metadata file [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/metadata.csv).


1. Visit https://timeor.brown.edu.
2. For (1) in 'Example Data' (side-bar) under 'Load raw data' click the 'SraRunTable & .fastq files' button. This will guide you through the 'Process Raw Data' tab demo. **Follow pop-ups and fill in grey boxes**. See [Run TIMEOR](#run-timeor-from-raw-data-starting-from-fastq-time-series-rna-seq) below for walk-through.
3. Next, for (2) in 'Example Data' (side-bar) under 'Load count matrix' click the 'Metadata & read count file' button. This will guide you through the rest of the full method demo. **Follow pop-ups and fill in grey boxes**. See [Run TIMEOR](#run-timeor-using-simulated-data-starting-from-read-count-matrix) below for full application walk-through.


Paper and Citation
=======

Please join a number of labs already using TIMEOR. [Read our paper here](https://www.biorxiv.org/content/10.1101/2020.09.14.296418v1), and cite:

Conard, A. M., Goodman, N., Hu, Y., Perrimon, N., Singh, R., Lawrence, C., & Larschan, E. (2020). TIMEOR: a web-based tool to uncover temporal regulatory mechanisms from multi-omics data. bioRxiv.

Overview
========

The [TIMEOR software (web, Conda, and Docker)](https://github.com/ashleymaeconard/TIMEOR) gives users a flexible and intuitive web platform with which to
upload their time-series RNA-seq data and protein-DNA data, and step through the entire temporal differential gene
expression and gene dynamics analysis pipeline to generate gene regulatory networks. The application is organized into three
separate stages: Pre-processing, Primary Analysis, and Secondary Analysis. [Click here](https://www.youtube.com/watch?v=bBodWvGPD6g&t=8s) for a quick video demonstration of the TIMEOR webserver, and [click here](https://www.youtube.com/watch?v=-BV9B3L1ymg) for a video guide through the webserver tutorial (below in [Run TIMEOR](https://github.com/ashleymaeconard/TIMEOR#run-timeor).

<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/timeor_steps.png" style="width:95.0%" />
</center>
<p>
 
</p>

A.  ***Pre-processing*: Gather and configure time series RNA-seq data**.
    The user can choose to process raw data (.fastq files) using a GEO
    identifier, or upload a raw count matrix (genes by samples). TIMEOR
    then automatically chooses from several methods with which to
    perform quality control, alignment, and produce a count matrix. The
    user can then choose between several methods to normalize and
    correct the data.

B.  ***Primary Analysis*: Use methods to perform differential gene
    expression analysis and determine gene trajectory clusters**. TIMEOR
    provides two continuous and one categorical DE method for the user.
    Specifically, DE genes are determined using one or more of
    [ImpulseDE2](https://bioconductor.org/packages/release/bioc/html/ImpulseDE2.html),
    [Next
    maSigPro](https://www.bioconductor.org/packages/release/bioc/html/maSigPro.html)
    and/or
    [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html),
    depending on your answers to questions in the Pre-processing stage.
    The user can then compare (via Venn diagram) the DE results between
    methods with a previous study of their choice to determine which DE
    method results to use for downstream analysis. The user can toggle
    between methods to determine which results produce the most
    resonable results to pass on to the next Secondary Analysis. TIMEOR
    then automatically clusters and creates an interactive clustermap of
    the selected DE gene trajectories over time. The user can choose a
    different number of clusters if desired.

C.  ***Secondary Analysis*: Assess enrichment, factor binding, and
    temporal relations**. The user can analyze the gene trajectory
    clusters using three categories of analysis in different tabs:
    *Enrichment*, identifies the genes and gene types that are
    over-represented within each cluster; *Factor Binding*, predicts which
    TFs are post-transcriptionally influencing the expression of each
    gene cluster using motif and ChIP-seq data; and *Temporal Relations*,
    identifies transcription factor (TF) gene regulatory network (GRN). <span style="color:#3F88DE"> **Blue arrow:** predicted TF to observed TF, experimentally determined interaction</span>. <span style="color:#D6678D"> **Pink arrow:** observed TF to observed TF, experimentally determined interaction</span>. <span style="color:#F7C144"> **Yellow arrow:** observed TF to observed TF, predicted interaction</span>. <span style="color:#5B8179"> **Green arrow:** predicted TF to observed TF, predicted interaction</span>. Network displayed in table format in app to enhance flexibility of GRN visualization.


Website
=======

### Computational Biology Core at Brown Univeristy and DRSC/TRiP Functional Genomics Resources at Harvard Medical School Partnership
TIMEOR is available online at https://timeor.brown.edu.

Installation
============

After cloning the TIMEOR repo from GitHub
(<a href="https://github.com/ashleymaeconard/TIMEOR.git" class="uri">https://github.com/ashleymaeconard/TIMEOR.git</a>)

### Conda Environment

1. Install Anaconda2/Miniconda2 (version 4.8.3). 
2. Type `$eval "$(/PATH/TO/<CONDA_DIR>/bin/conda shell.bash hook)"`
3. Locate packages.yml in /TIMEOR/
4. Type `$conda env create -f packages.yml`
5. Type `$conda activate timeor_conda_env`
6. Type `$R`
7. In R
    1. Type `>library(shiny)`
    2. Type `>runApp("app/", launch.browser=F)`
    3. Go to browser URL from R output
        (e.g. <a href="http://127.0.0.1:3681" class="uri">http://127.0.0.1:3681</a>
        )

Note, if TIMEOR is running on a remote machine, you may access the website through the ssh gateway. For example on port 8888, `ssh -L 8888:127.0.0.1:8888 USERNAME@SERVER -t host=SERVER -L 8888:127.0.0.1:3681`.

### Docker

1. Install Docker (version 20.10.0 recommended)
2. Building Docker image in TIMEOR directory:
    1. `$docker build -t timeor_env .`
3. Running Docker:
    1. `$docker run -p 3838:3838 timeor_env`
    2. Shiny server will be running on port 3838. Thus, in a browser visit `localhost:3838`.

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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T000.png" style="width:95.0%" />
</center>
<p>
     
</p>

1.  In the far-left side-bar click on “Example Data” and then
    under “Load raw data” click on "SraRunTable & raw .fastq files".

<p>
     
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T00.png" style="width:95.0%" />
</center>
<p>
     
</p>

2.  Follow the pop-up prompt to explore the default settings to questions 1-6 to set the adaptive 
    default parameters, and then click the "Run" button to begin retrieving the raw data (SRR8843738
    and SRR8843750), performing quality control, and aligning the reads using HISAT2 and Bowtie2.

<p>
     
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T01.png" style="width:95.0%" />
</center>
<p>
     
</p>



<p>
     
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T02.png" style="width:95.0%" />
</center>
<p>
     
</p>

3.  Once the data have been retrieved and quality control has finished. You can view a summary under
    the "Quality Control" panel on the right. Interactive results can be downloaded and viewed from MultiQC.

<p>
     
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T03.png" style="width:95.0%" />
</center>
<p>
     
</p>
 
<p>
     
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T04.png" style="width:95.0%" />
</center>
<p>
     
</p>   

4.  You can explore the alignment results between both methods in the "Alignment Quality" panel. Note that
    HISAT2 is splice-site aware). You can choose the method and then click "Generate count matrix" to have 
    TIMEOR generate the read count matrix for the next tab "Load Count Matrix". You will see a green check mark when finished. Follow pop-up. You have completed this section of the tutorial for loading raw time-series RNA-seq data.

<p>
     
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T05.png" style="width:95.0%" />
</center>
<p>
     
</p> 

<p>
     
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T050.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T1.png" style="width:95.0%" />
</center>
<p>
     
</p>

2.  Follow the pop-up to explore results on the Pre-processing
    tabs "Process Raw Data" and "Process Count Matrix".

<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T2.png" style="width:95.0%" />
</center>
<p>
 
</p>
<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T3.png" style="width:95.0%" />
</center>
<p>
 
</p>

3.  On the Normalize and Correct Data tab, choose from normalization and
    correction methods and click “Run” to view result. You have completed pre-processing. Follow pop-up.

<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T4.png" style="width:95.0%" />
</center>
<p>
 
</p>

<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T40.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T50.png" style="width:95.0%" />
</center>
<p>
 
</p>

<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T51.png" style="width:95.0%" />
</center>
<p>
 
</p>

6.  See 3-way Venn diagram in white box. You can compare a previous study list of genes with the three differential expression
    results by first downloading [prev\_study.txt](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/prev_study.txt). Then upload this file using the Browse button. See 4-way Venn diagram in blue box. 

<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T52.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T5.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T6.png" style="width:95.0%" />
</center>
<p>
 
</p>

13. Proceed to the "Factor Binding" tab to view the *observed* (that is, differentially expressed genes in data) and top
    predicted transcription factors in each gene cluster (under
    “Observed and Top Predicted Transcription Factors”). At least 40% of the transcription factor prediction methods must agree on their top predicted transcription factors, otherwise cell is left blank. Row names are gene expression trajectory clusters. **NOTE** for this demonstration the threshold is reduced.
    
<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T70.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T7.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T8.png" style="width:85.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T9.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T90.png" style="width:95.0%" />
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
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T10.png" style="width:95.0%" />
</center>
<p>
 
</p>

20. Your results folder can be downloaded on the far-left side-bar under "Download Results Folder". **NOTE** The original simulated data and results can be downloaded [here](https://github.com/ashleymaeconard/TIMEOR/tree/master/demos/simulated_data). 

<p>
 
</p>
<center>
<img src="https://github.com/ashleymaeconard/TIMEOR/blob/master/app/www/T11.png" style="width:95.0%" />
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
