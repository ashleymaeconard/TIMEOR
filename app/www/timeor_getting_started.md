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
- Please click each button just once.
- <span style="color:red"> Once analysis has begun, please proceed through TIMEOR sequentially</span>. **The user can visit previous tabs**, but proceed forward sequentially. Before beginning the analysis, the user can skim through each tab to see what is to come.
- **TIMEOR supports these types of time-series data**: 
  - control at 1st time point vs case (i.e. treatment) at subsequent time points
  - control **or** case over time
- Please compare one set of time-series experiments at a time.
  - TIMEOR v1 does **not** support time-series RNA-seq of these types:
      - control v.s. case 1 vs. case 2
  - In this case, simply perform the analysis **separately**:
      - control v.s. case 1
      - control v.s. case 2
- <span style="color:red"> Thank you for using TIMEOR! Please help us improve to better assist you. Please contact us with questions, ideas, and suggestions. If an error occurs, please download the log file (far left) to check. When contacting us with questions, please send the *time*, the log file, and if possible *a screenshot* so we know where in TIMEOR you are.</span>

## Two ways to input data:

1.   Import **SraRunTable from GEO**\* where TIMEOR will process raw data
    through retrieving .fastq files, quality control, alignment, and
    read count matrix creation. Read [this section]() for information about 
    this input specification. Read [this section](https://timeor.brown.edu/app_direct/timeor/timeor_app_tutorial.html#run-timeor-from-raw-data-starting-from-.fastq-time-series-rna-seq).

2.   Import **metadata file\*\* and count matrix \*\*\*** (skipping raw
    data retrieval, quality control, alignment, and read count matrix
    creation) and proceeding straight to normalization and correction. 
    Read [this section](https://timeor.brown.edu/app_direct/timeor/timeor_app_tutorial.html#run-timeor-using-simulated-data-starting-from-read-count-matrix).

Then simply follow the prompts. Fill out the **grey** boxes to begin
interacting with each stage and tab. 

#### Input file types:

  \* **SraRunTable from GEO** follow instructions in TIMEOR first tab
    (“Process Raw Data”)

   \*\* **metadata file** requires at least these columns.
    -   *ID, condition, time, batch*
        -   *ID*: a unique identifier (ID) for the user
            (e.g. case1min\_rep1)
        -   *condition*: one word description (e.g. case, control)
        -   *time*: numerical values e.g. (0, 20, 40)
        -   *batch*: string description of batch (e.g. b1, b2, b3)

  \*\*\* **count matrix** : rows should be unique gene identifiers
    (e.g. Flybase, Ensembl or Entrez IDs) and columns should be the IDs
    from metadata file.

## Input Data: Metadata File
- Make sure there are an equal number of replicates for each sample.
- Label control as “control”.
- Have time point data in order where control or time point 1 is at the top and the last time point is at the bottom.
- have unique IDs follow these formats:
  - NAMETIME_REPLICATE_TIME 
- Sometimes when using editors such as Excel, odd delimiters specific to the user's machine are added at the end of lines. We advise users to check that these are not present.

## Input Data: Read Count Matrix
- Pre-filter out any rows you are not interested to process (such as low count genes across all samples).
- Gene ID column should be named “ID” and populated with Ensembl IDs.
- Make sure columns are in the same order as the rows of metadata file.

## Inputs Data, SraRunTable
- Make sure that the resulting metadata file meets the requirements for the input metadata file.

## Input Data, .fastq Files
- compressed fastq files?

## Suggestions for How to Answer Six "Determine Adaptive Default Methods" Questions
"Please select at least the organism, sequencing, and experiment type, then load metadata or SraRunTable.txt.",
"Not applicable" = "matrix"
Keep “multiple method comparisons” true to see comparison plots and results to show TIMEOR’s full functionality

## Method and question choice assistance
- Types of time-series
- Comparing multiple methods for differential expression
  - Note that there is not a fold change cut-off, only an adjusted p-value cutoff. This allows the user to view significant differences in expression while the fold change might be smaller. This is useful to observe changes for genes including non-coding genes and genes involved in dosage compensation. 

- Normalization
  - Upper Quartile 
  - Spearman
  - Pearson
- Cluster choice
- Talk about types of clusters and outputs

## Going Through Each Tab
Enrichment needs to be performed for each cluster to go to the next tab.

