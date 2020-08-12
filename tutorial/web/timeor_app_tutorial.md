Introduction
------------

Analyzing time course differential gene expression of multi-omics data
is computationally expensive, laborious, and sometimes fallible,
particularly for those without knowledge of Linux or Unix. The TIMEOR
web-application is an interactive R-Shiny web interface that takes raw
.fastq files and performs all analysis from quality control and
differential gene expression to gene regulatory network reconstruction,
allowing users to investigate the altered regulatory gene networks from
an experiment. Through a suite of published and novel methods, the
TIMEOR application retrieves and loads raw .fastq data from NCBI to be
processed using adaptive default methods in response to six questions.
TIMEOR then compares several alignment and differential expression
methods automatically (both close and distant timepoint options) from
which the optimal set of interactive results, visualizations, and
publication-ready figures are generated to highlight the transcription
factor regulatory network. Additionally, TIMEOR outputs geneontology,
motif, and ChIP-seq results to help validate the regulatory network.
User outputs and data are stored in the cloud and accessible via a
private user key. The web server is completely free, and accessible at
[timeor.org](timeor.org).

The TIMEOR Application
----------------------

### Overview

The TIMEOR software gives users an intuitive web platform with which to
upload their raw data and step through the entire differential gene
expression analysis pipeline. The application is organized into three
separate pages, Load Stage, Primary Analysis, and Secondary Analysis.

<p>
 
</p>
<center>
<img src="timeor%20steps.png" style="width:95.0%" />
</center>
<p>
 
</p>
1.  **Load Stage: Gather and process time course RNA-seq data from
    <https://www.ncbi.nlm.nih.gov/geo/> corresponding to the desired
    data**. TIMEOR will retrieve and store the relevant raw .fastq
    datasets in your generated user directory. After answering six
    questions, each .fastq file is run through quality control, and a
    count matrix is generated, which you can then use to correlate your
    samples (i.e. replicates), as well as normalize, correct, and
    perform Principal Component Analysis (PCA) on your data.

2.  **Primary Analysis: Use methods to perform differential gene
    expression analysis and determine gene cluster dynamics**. Depending
    on your inputs from 1,
    [ImpulseDE2](https://bioconductor.org/packages/release/bioc/html/ImpulseDE2.html)
    and/or
    [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)
    is run to perform differential expression. You can then compare the
    results of these methods with interactive data visualizations and
    publication-ready figures.

3.  **Secondary Analysis: Use the results of 2 to do more with your
    data**. Use the differential expression results to uncover perturbed
    mechanisms in each experiment. Via gene ontology (GO) and pathway
    analysis, you can uncover the gene network which follow the same
    expression dynamics, perform de novo motif identification for each
    gene cluster, reconstruct the transcription factor regulatory
    network, and validate the network with ChiP-seq data provided
    through TIMEOR.

### Architecture

TIMEOR was built on a 2.3 GHz Intel Core i5 Processor running OS X.
[timeor.org](timeor.org) is hosted at the [DRSC/TRiP Functional Genomics
Resources](https://fgr.hms.harvard.edu/) at the Harvard Medical School.
The application is run on an instance of an Ubuntu EC2 quad-core virtual
machine (find out exactly from Perrimon) and user data is stored in S3
storage buckets (change). A unique file system is granted for each user
using an AWS S3 Bucket (change); the directories and any data stored in
them will persist for two weeks after the user originally launches the
instance. Dependencies are handled with a Docker container and
Anaconda’s environment management system. The software was built using
RStudio Version 1.2.5001 using the R-Shiny package.

<center>
<img src="TIMEOR%20full%20pipeline.png" style="width:95.0%" />
</center>
<p>
 
</p>
### Starting the application

To use the TIMEOR pipeline, visit [timeor.org](timeor.org).

First you will be prompted to enter in your email address. Shortly
after, you will be emailed a link to a file system containing the
intermediate outputs and results of your analysis. This directory will
be accessible over the course of two weeks so that you can pause your
progress and come back to it. Please note that closing the window will
eliminate all progress made on the application side; however, your
intermediate outputs will still persist in the directory sent to you.

If you are installing the application locally, make sure you have the
latest version of anaconda – Python 3.7 – installed on your machine. You
also should place the timeor\_env.yml file in the same directory from
which the application is run. Then, in your terminal, navigate to this
directory and run:
<p>
 
</p>
`source ~/anaconda/bin/activate`, then `conda activate timeor_env`

<p>
 
</p>
### Load Stage

#### Process Data

To load your data, make sure you are on the tab “Process Data”
<p>
 
</p>
<center>
<img src="process_data.png" style="width:95.0%" />
</center>
<p>
 
</p>
**Search and Retrieve**. If you know the GEO IDs associated with the
datasets you wish to analyze, upload them in a text file as shown in the
image below. To locate your GEO IDs for your study of interest:

1.  Go to GEO <https://www.ncbi.nlm.nih.gov/geo/> to search for the
    dataset of interest, e.g. “time series RNA-seq Drosophila S2R+” (A)
2.  Click on the SRA number in “Relations”.
3.  In the top right hand corner select “Send to” then “Run Selector”
    then click “Go”.
    <p>
     
    </p>

<center>
<img src="load_data.png" alt="Example NCBI screen." style="width:95.0%" />
</center>
<p>
 
</p>
1.  In the new window in the “Select” box under “Download” click
    “metadata”.
    <p>
     
    </p>

<center>
<img src="sraRunSelector.png" alt="Example SRA Run Selector." style="width:95.0%" />
</center>
<p>
 
</p>
1.  If personal data, upload metadata. Otherwise, upload the
    resulting“SraRunTable.txt” as the metadata file.
    <p>
     
    </p>

<center>
<img src="step5.png" style="width:95.0%" />
</center>
<p>
 
</p>
1.  Upload the SRR\_Acc\_List.txt or any text file containing the SRA
    number(s).
    <p>
     
    </p>

Next, answer six questions located on the right, from which the optimal
settings for your differential gene expression analysis will be
configured.
<p>
 
</p>
<center>
<img src="process_data.png" style="width:95.0%" />
</center>
<p>
 
</p>
Once you click the “Run” button, the following will happen to your data:

1.  The metadata file uploaded will be parsed and the necessary results
    file structure for further analysis will be created.

2.  The FASTQ files associated with the GEO IDs you uploaded will be
    loaded.

3.  Each replicae will be aligned to the genome using at least one of
    two alignment methods,
    [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and
    [HISAT2](http://daehwankimlab.github.io/hisat2/). If you answer yes
    to question 4, the comparison of the two methods’ results, and the
    method with the highest median overall percent, are displayed on the
    front page.

4.  A gene name by replicate (i.e. sample) matrix detailing the genes
    and their expression levels for each time point will be produced and
    stored in your directory.

#### Load Count Matrix

Once the above steps are completed, you will be ready to upload your
gene count table for further analysis. Navigate to the top of the page
and click “Load Count Matrix”
<p>
 
</p>
<center>
<img src="load_button.png" style="width:95.0%" />
</center>
<p>
 
</p>
Upload your gene count matrix and wait for the matrix to be displayed
before moving further.
<p>
 
</p>
<center>
<img src="count.png" style="width:95.0%" />
</center>
<p>
 
</p>
Once the above plots and count matrix are generated, you may interact
with them as you wish. Your unnormalized and uncorrected replicates will
be displayed using both interactive principal component analysis (PCA)
and correlation plots. To correlate your replicates (i.e. samples), you
have the option to choose between Pearson and Spearman correlation
coefficients and can choose the one that best suits the distribution of
the gene counts.

#### Normalization and Correction

Upon seeing the initial interactive PCA results, you may want to
normalize and correct your data. To do so, please click on the “PCA
After Normalization and Correction” tab. Your screen should now resemble
the image below. TIMEOR currently supports the
[Harman](https://www.bioconductor.org/packages/devel/bioc/vignettes/Harman/inst/doc/IntroductionToHarman.html)
method for correction and [Trimmed Mean of
M-Values](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)
and
[Upper-quartile](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94)
methods for normalization. To run these methods on your data, select the
ones that best suit your experiment and click the “Run” button. After a
few seconds, the interactive normalized and corrected PCA plots will
display on your screen. Now these data can be used to determine the
temporally differentially expressed genes.
<p>
 
</p>
<center>
<img src="norm_and_correct.png" style="width:95.0%" />
</center>
<p>
 
</p>
### Primary Analysis

The goal of Primary Analysis is to compare multiple close and distant
time point differential expression methods in order to understand which
gene dynamis were either altered by an experiment or changed through
time. Once you are ready to perform this step of the pipeline, navigate
to the button that says “Primary Analysis” on the left.

<p>
 
</p>
<center>
<img src="prim_button.png" style="width:95.0%" />
</center>
<p>
 
</p>
In the top left window, enter in two pieces of information.

1.  The name of the folder in which you would like your results to be
    stored. This will determine the names of the results files as well.
2.  The p-value threshold at which the differential expression methods
    stop considering expression to be statistically significant. Note, a
    p-value of 1 is used as time points have been truncated and there
    are fewer significant genes.

Now you are ready to run differential expression. Press “Go” — if the
answer to question 4 is "yes,
[ImpulseDE2](https://bioconductor.org/packages/release/bioc/html/ImpulseDE2.html)
and
[DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)
will both be run. The results are stored in the directory the user
selected above.
<p>
 
</p>
<center>
<img src="run_de.png" style="width:95.0%" />
</center>
<p>
 
</p>
If you would like to compare the differential expression results, TIMEOR
will generate a venn disagram that visualizes the number of
differentially expressed genes found between ImpulseDE2, DESeq2, and an
optional, user-defined previous study. You can generate this venn
diagram by clicking the “Render Venn Diagram” button. If you would like
to compare the results with a previous study, please upload a text file
containing the results and re-render the venn diagram.

<p>
 
</p>
<center>
<img src="venn_results.png" style="width:95.0%" />
</center>
<p>
 
</p>
You may select the method whose results are to be displayed once both
methods finish as a matrix below. By default, DESeq2’s results will be
displayed.

<p>
 
</p>
<center>
<img src="which_method.png" style="width:95.0%" />
</center>
<p>
 
</p>
Once both methods finish running, your screen should look something like
the image below. As shown below, a clustermap that clusters the genes
using a Euclidean distance and Wards.D2 method is generated. Users can
choose the number of clusters they would like to see, ranging from 0 to
15 clusters. In the former case, TIMEOR will choose the optimal number
of clusters using the mode of 3 measures: silhouette, K-means, and
Calinski Harabasz criterian.
<p>
 
</p>
<center>
<img src="de_results.png" style="width:95.0%" />
</center>
<p>
 
</p>
