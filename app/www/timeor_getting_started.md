---
title: "Welcome! Thank you for using TIMEOR."
output:
  html_document:
    css: github.css
---

## Quick Start: 3 Steps
 TIMEOR accepts 2 input types: (1) raw .fastq files and SraRunTable [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/real_data_subset/timeor/data/SraRunTable.csv) or a (2) RNA-seq time-series read count matrix [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/countMatrix.csv) and metadata file [(e.g. here)](https://github.com/ashleymaeconard/TIMEOR/blob/master/demos/simulated_data/timeor/data/metadata.csv).

1. Visit https://timeor.brown.edu.
2. For (1) in 'Example Data' (side-bar) under 'Load raw data' click the 'SraRunTable & .fastq files' button. This will guide you through the 'Process Raw Data' tab demo. **Follow pop-ups and fill in grey boxes**. See [Run TIMEOR](#run-timeor-from-raw-data-starting-from-.fastq-time-series-rna-seq) for walk-through.
3. Next, for (2) in 'Example Data' (side-bar) under 'Load count matrix' click the 'Metadata & read count file' button. This will guide you through the rest of the full method demo. **Follow pop-ups and fill in grey boxes**. See [Run TIMEOR](#run-timeor-using-simulated-data-starting-from-read-count-matrix) for full application walk-through.

## Important Points to Remember

- <span style="color:red">We strongly encourage the user to input a read count matrix when possible</span>, as the input file size limit is 10GB.
- For larger dataset processing, the user is encouraged to use our ready-to-use Docker image. Read steps [here](). If that is not possible, please feel free to contact us for specific space requirements. **We are happy to help**. 
- <span style="color:red">While TIMEOR analysis is running, simply make sure to revisit the page at least once an hour.</span>
- <span style="color:red">We strongly encourage the user to keep '5.Compare multiple methods' set to 'Yes' to see TIMEOR's full functionality.
- Please click each button just once.
- <span style="color:red"> Once analysis has begun, please proceed through TIMEOR sequentially</span>. **The user can visit previous tabs**, but proceed forward sequentially. Before beginning the analysis, the user can skim through each tab to see what is to come.
- **TIMEOR supports these types of time-series data** (note this is asked in Question 3 of "Determine Adaptive Default Methods"): 
  - control at 1st time point vs case (i.e. treatment) at subsequent time points
  - control **or** case over time
- Please compare one set of time-series experiments at a time.
  - TIMEOR v1 does **not** support time-series RNA-seq of these types:
      - control v.s. case 1 vs. case 2
  - In this case, simply perform the analysis **separately**:
      - control v.s. case 1
      - control v.s. case 2
- <span style="color:red"> Thank you for using TIMEOR! Please help us improve to better assist you. Please contact us with questions, ideas, and suggestions. If an error occurs with your data, please download the log file (far left) to check. When contacting us with questions, please send the *time*, the log file, and if possible *a screenshot* so we know where in TIMEOR you are.</span>

## Two ways to input data:

1.   Import **SraRunTable from GEO**\* where TIMEOR will process raw data
    through retrieving .fastq files, quality control, alignment, and
    read count matrix creation. Read first tab of TIMEOR (Getting Started) for information about 
    this input specification. Read [this section](https://timeor.brown.edu/app_direct/timeor/timeor_app_tutorial.html#run-timeor-from-raw-data-starting-from-.fastq-time-series-rna-seq) for information about how to process these data in TIMEOR. **We strongly encourage users to upload a read count matrix**, or process raw .fastq data via TIMEOR's interface locally using Docker ([see 5 steps here](https://timeor.brown.edu/app_direct/timeor/timeor_app_tutorial.html#local-installation)).

2.   Import **metadata file\*\* and count matrix \*\*\*** (skipping raw
    data retrieval, quality control, alignment, and read count matrix
    creation) and proceeding straight to normalization and correction. Read [this section](https://timeor.brown.edu/app_direct/timeor/timeor_app_tutorial.html#run-timeor-using-simulated-data-starting-from-read-count-matrix) for information about how to process these data in TIMEOR.

Then simply follow the prompts. Fill out the **grey** boxes to begin
interacting with each stage and tab. 

#### Input file types:

  - \* **SraRunTable from GEO** requires at least these columns (which will be reordered to produce the metadata file).
    - *treatment, time, Run, replicate, batch*
        - *treatment*: one word describing experiment
        - *time*: numerical values e.g. (0, 20, 40)
        - *replicate*: one word description of replicate

  - \*\* **metadata file** requires *at least* these columns.
    -   *ID, condition, time, replicate*
        -   *ID*: a unique identifier (ID) for the user
            (e.g. case1min\_rep1)
        -   *condition*: one word description (e.g. case, control)
        -   *time*: numerical values e.g. (0, 20, 40)
        -   *replicate*: one word description of replicate (e.g. b1, b2, b3)

  - \*\*\* **count matrix**  requires Ensembl or Flybase unique gene identifiers, and columns should be the IDs
    from metadata file, and in the same order as metadata file.

## Inputs Detailts: SraRunTable
- Please upload the SraRunTable.txt, which has comma delimiters.
- Please remove as many unneeded colums as you can. Some dataset SraRunTables have odd delimiters that are difficult to parse.
- Make sure that the resulting metadata file meets the requirements for the input metadata file [(link here)](https://timeor.brown.edu/app_direct/timeor/timeor_getting_started.html#input-data-read-count-matrix).
- For paired-end reads, the pairs must have "_NUM" to distinguish them (e.g. SRRXXX_1.fastq.gz, SRRXXX_2.fastq.gz).

## Input Data: .fastq Files
- The .fastq files are downloaded as .fastq.gz.

## Input Details: Metadata File
- Make sure there are an equal number of replicates for each sample.
- Label control as “control”.
- Have time point data in order where control or time point 1 is at the top and the last time point is at the bottom.
- have unique IDs follow these formats:
  - NAMETIME_REPLICATE_TIME 
- Sometimes when using editors such as Excel, odd delimiters specific to the user's machine are added at the end of lines. We advise users to check that these are not present.
- Please upload .csv files.

## Input Data: Read Count Matrix
- Pre-filter out any rows you are not interested to process (such as low count genes across all samples).
- Gene ID column should be named “ID” and populated with Ensembl IDs.
- Make sure columns are in the same order as the rows of metadata file.
- Please upload .csv files. 

## Suggestions for How to Answer Six "Determine Adaptive Default Methods" Questions
- Overall the user must select at least the organism, sequencing, and experiment type, then load metadata or SraRunTable.txt.

- Question 1 asks: "What type of organism?" The user can choose from fruit fly, human, or mouse.

- Question 2 asks: "What type of sequencing?" If the user is uploading <span style="color:red">a read count matrix, strongly encouraged</span>, the user can choose "not applicable".

- Question 3 asks: "What type of experiment?" There are two options - "case vs. control", and "just case or control" types of time-series that TIMEOR supports [(see this section)](https://timeor.brown.edu/app_direct/timeor/timeor_getting_started.html#important-points-to-remember).

- Question 4 asks: "What type of time-series?" There are three options - "close time point and long time series", "close time point and short time series", and "distant time point". Based on the user's understanding of the biological system, the user should decide whether the timepoints are considered close or far in time. This question is important to determine how to model differential gene expression (DE) trajectories over time. 
    * DESeq2 is a categorical DE method generally used to analyze timepoints separately. When time points are far apart this is a good option. TIMEOR uses DESeq2 if the user toggles to "distant time point". 
    * When we are interested to model gene trajectories (when time points are close), we assess the temporal dynamic expression between time point $t$ given $t-1$. In this context, it is advised a continuous DE method. TIMEOR uses ImpulseDE2 if the user toggles to "close time point ...". ImpulseDE2 employs an impulse model to determine differentally expressed genes. 
    * Importantly, it is <span style="color:red">**strongly**</span> advised to compare all three (ImpulseDE2, Next maSigPro, and DESeq2) DE methods' results (<span style="color:red">by keeping 'Yes' for Question 5</span>), especially when there are "close time points and short time-series". Recent studies such as Spies et al. 2019 show that DESeq2 performs well when determining differentially expressed genes when time-series is short. To compare all three these, keep 'Yes' as the answer for Question 5 (below).
    * *Please see our manuscript for more a more robust explanation and a series of citations for further reading.*

- Question 5 asks: "Compare multiple methods (alignent and differential expression)?" If this question is left to 'Yes' (which is **strongly** encouraged), TIMEOR will run all methods for the user to determine the best suited method. This is important because in many cases the categorical method DESeq2 which does not consider gene trajectories, still returns a robust set of differentially expressed genes. If this is set to 'No', TIMEOR will run for alignment (if applicable): HISAT2, and for DE: DESeq2 (if distant time points selected in Question 4), or ImpulseDE2 (if close time points selected in Question 4).

- Question 6 asks: "What is the maximum number of time steps over which one gene can influence the transcription of another gene?" This question prompts the user to tell TIMEOR the timespan over which one gene can *directly* influence another. Within this timespan all interactions are considered. It is advised to keep this value small if the time points are spaced out. Said differently, at each time point $t$ for a differentially expressed gene $g$, if Question 6's answer were 2, TIMEOR would be asking, what are potential interactions of $g$ with other TFs across $t+1$ and $t+2$.

## Method and question choice assistance

- "Normalize and Correct" tab: there are two normalization options - upper quartile and trimmed mean of M-values. It is advised to try both methods through TIMEOR's interactive interface because the influence of normalization differs depending on the RNA-seq data structure. 
    * There are several recent papers that discuss these differences such as Zyprych-Walczak et al. 2015, Pereira et al. 2018, and Abbas-Aghababazadeh et al. 2018.

- "Normalize and Correct" tab: there are two options for correlating samples/replicates using the Pearson or Spearman correlation. The choice of correlation method depends heavily on the assumptions the user wants to make about their data, and it is encouraged to try both in TIMEOR's interface. The user knows more about which samples/replicates (e.g. time points) should cluster together and how to identify outliers. 
    * Both correlation methods define the strength of the relationship between the samples/replicates. The Pearson correlation accounts for differences in the samples/replicates mean and standard deviation when defining the linear relationship. The Spearman correlation is actually a nonparametric measure that uses the rank values of the samples/replicates. 
    * Importantly, the more similar the expression profiles between samples/replicates, the higher the correlation coefficient will be.
    * Furthermore, the user is encouraged to remove any outliers (if needed) for further analysis.

- "Primary Analysis" stage: the user can choose to allow TIMEOR to automatically cluster the DE gene trajectories, or the user can choose the number of gene trajectory clusters. Importantly, finding the optimal solution to this hierarchical clustering problem is an NP-hard. Thus, user input is needed to assess a reasonable number of clusters for downstream analysis. To help, TIMEOR provides an automatic clustering option (PDF visible when folder is downloaded) which takes the mode between three unsupervised clustering methods (partition around medoids (Reynolds et al. 2006), Silhouette (Rousseeuw et al. 1987), and Calinski criterion (Calinkski et al. 1974)) to automatically return the number of gene trajectory clusters to the user. TIMEOR also provides an Elbow plot to show the user how the explained variation changes as a function of the number of clusters. The user can leverage this plot by picking the elbow of the curve. The user is encouraged to use the interactive clustermap and the clustering plots (available on download) to determine whether the automatic clustering option provides suitable clusters.

- "Primary Analysis" stage: NOTE, there is not a fold change cut-off for the DE gene trajectories, only an adjusted p-value cutoff. This allows the user to view significant differences in expression trajectories while the fold change might be smaller. This is useful to observe changes for genes including non-coding genes and genes involved in dosage compensation. 

- "Secondary Analysis: Factor Binding": the user is encouraged to "see each method's predicted transcription factors" and search for protein-DNA data (in .bigWig format) to view the binding profile of that transcription factor across each gene trajectory cluster.

- "Secondary Analysis: Temporal Relations": the user can add additional genes or transcription factors (potentially viewed on Factor Binding tab) to the final gene regulatory network (GRN) within STRINGdb. NOTE: TIMEOR only reports the TF GRN using the "Observed and Top Predicted Transcription Factors" and the user is **encouraged** to view the results from individual methods (on Factor Binding tab) when constructing the final GRN.

