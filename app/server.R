# server.R
# Ashley Mae Conard
# Last Mod. 6 June 2020

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(autoplotly)
library(heatmaply)
library(d3heatmap)
library(gplots)
library(data.table)
library(dplyr)
library(ggfortify)
library(stats)
library(factoextra)
library(reticulate)
library(shinyLP)
library(shinyalert)
library(fpc)
library(stringr)
library(cluster)
library(vegan)
library(promises)
library(future)
plan(multiprocess)

# Set app dir
app_dir <- getwd()

# Set maximum file size upload to app
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)

source(paste(app_dir, "/scripts/clustermap_function.r", sep = ""))
source(paste(app_dir, "/scripts/DESeq2.r", sep = ""))
source(paste(app_dir, "/global.R", sep = ""))
#source(paste(app_dir, "/scripts/stringdb_indiv_cluster.r"))

# Runs this server code with each new session
function(input, output, session) {
  # Setting all reactive values
  # Data input
  local_results_folder  = reactiveVal(tempdir()) # tempdir()
  sim_demo_data         = reactiveVal(FALSE) # use simulated (demo) data
  real_demo_data        = reactiveVal(FALSE) # use real (demo) data
  user_input_data       = reactiveVal(FALSE) # user inputs data
  metadata_filepath     = reactiveVal() # uploaded metadata or created metadata file path from uploaded SraRunTable_file
  metadata_df           = reactiveVal() # metadata dataframe
  countMatrix_filepath  = reactiveVal() # count matrix file path
  countMatrix_df        = reactiveVal() # count matrix dataframe
  norm_corr_countMatrix_filepath = reactiveVal() # normalized and corrected count matrix filepath
  sra_input             = reactiveVal(FALSE) # if sra_input = T, user input raw data
  metadata_input        = reactiveVal(FALSE) # if metadata_input = T, user input count matrix

  # Adaptive defaults
  organism_ad           = reactiveVal()
  sequencing_ad         = reactiveVal()
  timeSeries_ad         = reactiveVal()
  experiment_ad         = reactiveVal()
  multiMethods_ad       = reactiveVal()
  trans_time_ad         = reactiveVal()
  
  # Pre-processing
  run_process_raw       = reactive({
    list(input$run_adaptive_defaults, sra_input())
  }) # if run pressed SraRunTable is input
  retrieveDone          = reactiveVal(FALSE)
  qcDone                = reactiveVal(FALSE)
  qcTableDone           = reactiveVal(FALSE)
  fastqOrgDone          = reactiveVal(FALSE)
  alignHdone            = reactiveVal(FALSE)
  plotAlignH            = reactiveVal(FALSE)
  run_Bowtie            = reactive({list(qcDone(), alignHdone())})
  alignBdone            = reactiveVal(FALSE)
  plotAlignB            = reactiveVal(FALSE)
  whichAligner          = reactiveVal('hisat2')
  #genCountMat           = reactiveVal(FALSE)
  countMatDone          = reactiveVal(FALSE)
  
  # Primary Analysis tab folders
  analysis_folder_name  = reactiveVal()
  adj_p_analysis_folder = reactiveVal()
  
  # Primary Analysis differential expression (DE) method results files
  impulseDE2File        = reactiveVal()
  nextMaSigProFile      = reactiveVal()
  deSeqFile             = reactiveVal()
  impulseDE2_out_file   = reactiveVal()
  deSeq_out_file        = reactiveVal()
  diffExpDone           = reactiveVal(0)
  numClusters           = reactiveVal()

  # Secondary Analysis 
  avg_prof_name         = reactiveVal()
  
  shinyalert("Welcome, it's about time!
              Quick start:",
              "TIMEOR accepts 2 input types: 
              (1) raw .fastq files
              (2) read count matrix. 
              
              For (1) in 'Example Data' (side-bar) under

              'Load raw data' click 'SraRunTable & .fastq files' button.
              
              This will guide you through the 'Process Raw Data' tab demo.",
            type = "info"
            )

  #################### If User Clicks on Demo (Simulated or Real) Data ###################

  # Use simulated data
  observeEvent(input$runSim, {
    print("running simulation")
    shinyalert(
      "Read Count Matrix Data and Metadata File Uploaded",
      "Metadata, adaptive defaults methods, and 
      read count matrix are loaded. 

      1) Go to 'Process Count Matrix' tab to see raw data metrics.
      2) Then go to 'Normalize and Correct Data' tab and 
      fill in the grey box.",
      type = "success"
    )
    
    # Use simulated metadata and count matrix (from demo folder)
    sim_demo_data(TRUE)
    real_demo_data(FALSE)
    metadata_input(TRUE)
    
    # Update switch to show metadata is being used
    updateSwitchInput(
      session = session,
      inputId = "sra_or_meta",
      value = TRUE,
      disabled = TRUE
    )
    
    # Set to simulated results directory
    local_results_folder(paste(app_dir, "/../demos/simulated_data/", sep = ""))
    
    # Enable buttons disabled from Load raw data button
    enable("runCorrection")
    enable("runDE")
    enable("renderVen")
    enable("prevStudy")
    enable("inputNumClust")
    enable("bigWigGo1")
    enable("bigWigGo2")
    enable("bigWigGo3")

    # Disable metadata file upload (using simulated data)
    disable("metadataFile")
    
    # Set adaptive default parameters and disable change
    updateTextInput(session, "organism", value = paste("dme"))
    disable("organism")
    updateTextInput(session, "sequencing", value = paste("se"))
    disable("sequencing")
    updateTextInput(session, "experiment", value = paste("cc"))
    disable("experiment")
    updateTextInput(session, "timeSeries", value = paste("ctstc"))
    disable("timeSeries")
    updateTextInput(session, "multipleMethods", value = paste("Yes"))
    disable("multipleMethods")
    updateTextInput(session, "trans_time", value = paste("2"))
    disable("trans_time")

    # Disable run (quality control button, choosing alignment method)
    toggleState("qcDownload", condition = F)
    disable("runGenCountMat")
    
    # Disable run (primary analysis: preprocessing) and read count matrix (primary analysis: load count matrix)
    disable("run_adaptive_defaults")
    disable("count_mat")
    
    # Set default analysis results folder name
    analysis_folder_name("simulated")
    
    # Insert results folder name in primary analysis and disable change
    updateTextInput(session, "analysisFolder", value = paste(analysis_folder_name()))
    disable("analysisFolder")
    
    # Set default adjusted p-value threshold
    adj_p_analysis_folder("0.05")
    
    # Insert adjusted p-value in primary analysis and disable change
    updateTextInput(session, "pValue", value = paste(adj_p_analysis_folder()))
    disable("pValue")
    
    # Disable toggle between differential expression methods
    disable("whichMethodInput")
    
    # Disable toggle between number of clusters
    disable("numClusts")

    # Set average profile for 1 transcription factor (pho) as example
    avg_prof_name("pho")
    updateTextInput(session, "tf1", value = paste("example: pho"))
    disable("tf1")
    disable("bigWig1")
    disable("bigWigGo1")
    
  }, ignoreInit = TRUE)
  
  # Render metadata table from SIMULATED (DEMO) data
  observeEvent(sim_demo_data(), {
    if (sim_demo_data() == TRUE) {
      # Create metadata dataframe
      df <- create_metadata_df()
      
      # Output metadata dataframe
      output_metadata_df(df)
      
      # Load .bigWig file
      #output_metadata_df
    }
  })
  
  # Use real data
  observeEvent(input$runReal, {
    shinyalert(
      "Raw Data and SraRunTable Uploaded",
      "SraRunTable will be converted to metadata, adaptive default methods will be set. 
      
      Press 'Run' to process raw RNA-seq time series data. 
      
      Follow pop-ups and fill in grey boxes.",
      type = "success"
    )
    
    # Use real data (so SraRunTable will be automatically converted to metadata) and count matrix (from demo folder)
    real_demo_data(TRUE)
    sim_demo_data(FALSE)
    sra_input(TRUE)
    
    # Update switch to show metadata is being used
    updateSwitchInput(
      session = session,
      inputId = "sra_or_meta",
      value = FALSE,
      disabled = TRUE
    )
    
    # Set to simulated results directory
    local_results_folder(paste(app_dir, "/../demos/real_data_subset/", sep = ""))
    
    # Disable metadata file upload (using simulated data)
    disable("metadataFile")
    
    # Set adaptive default parameters and disable change
    updateTextInput(session, "organism", value = paste("dme"))
    organism_ad("dme")
    disable("organism")
    updateTextInput(session, "sequencing", value = paste("se"))
    sequencing_ad("se")
    disable("sequencing")
    updateTextInput(session, "experiment", value = paste("cc"))
    experiment_ad("cc")
    disable("experiment")
    updateTextInput(session, "timeSeries", value = paste("ct"))
    timeSeries_ad("ct")
    disable("timeSeries")
    updateTextInput(session, "multipleMethods", value = paste("Yes"))
    multiMethods_ad("Yes")
    disable("multipleMethods")
    updateTextInput(session, "trans_time", value = paste("2"))
    trans_time_ad("2")
    disable("trans_time")
    
    # Disable load count matrix and correction
    disable("count_mat")
    disable("runCorrection")
    
    # Disable differential expression
    disable("analysisFolder")
    disable("runDE")
    disable("pValue")
    disable("renderVen")
    disable("prevStudy")
    
    # Disable toggle between differential expression methods
    disable("whichMethodInput")
    
    # Disable toggle between number of clusters
    disable("numClusts")
    disable("inputNumClust")
    
    # Disable factor binding Gos
    disable("bigWigGo1")
    disable("bigWigGo2")
    disable("bigWigGo3")

  })
  
  # Render SraRunTable from REAL (DEMO) data
  observeEvent(real_demo_data(), {
    if (real_demo_data() == TRUE) {

      # Create metadata dataframe
      df <- create_metadata_df()
      
      # Output metadata dataframe
      output_metadata_df(df)
      
      # Load .bigWig file
      #output_metadata_df
    }
  })
  
  #################### Pre-Process Stage ##################
  ##################### Process Raw Data ##################
  
  ## "Search and Retrieving" data box (top left had corner)
  
  ## Create Metadata ########
  
  # Read in SraRunTable or Metadata input file from .csv
  read_csv_sra_or_metadata_file <- reactive({
    print("read_csv_sra")
    
    # Get the metadata file from simulated_data folder
    if (input$runSim) {
      metadata_filepath(paste(
        local_results_folder(),
        "/timeor/data/metadata.csv",
        sep = ""
      ))
      
      # Get the SraRunTable file from real_data folder
    } else if (input$runReal) {
      metadata_filepath(paste(
        local_results_folder(),
        "/timeor/data/SraRunTable.csv",
        sep = ""
      ))
      
      # Set input folder to be user defined
    } else{
      metadata_filepath(input$metadataFile$datapath)
    }
    
    if (is.null(metadata_filepath())) {
      return()
    }
    read.csv(file = metadata_filepath(),
             sep = ",",
             header = TRUE)
  })
  
  # Check metadata file is in correct format and then save as data frame
  create_metadata_df <- function() {
    # Read .csv (SraRunTable or metadata file) and cast as dataframe
    sra_metadata_df <-
      as.data.frame(read_csv_sra_or_metadata_file())
    metadata_df(sra_metadata_df)
    
    # If metadata used
    if (sim_demo_data() == TRUE ||
        input$sra_or_meta == TRUE) {
      # if sim data used or user chooses metadata
      print("Metadata toggle")
      metadata_input(TRUE)
      
      # Check format of metadata file
      if (!(any(grepl(
        "ID", colnames(sra_metadata_df), ignore.case = TRUE
      )) &
      any(grepl(
        "condition", colnames(sra_metadata_df), ignore.case = TRUE
      )) &
      any(grepl(
        "time", colnames(sra_metadata_df), ignore.case = TRUE
      )) &
      any(grepl(
        "batch", colnames(sra_metadata_df), ignore.case = TRUE
      )))) {
        shinyalert(
          "Please check Metadata.",
          "1) Toggle 'Input file' to 'Metadata' in question 4 above.
                 2) Make sure metadata file includes these columns:
                   'id', 'condition', 'time', batch'.",
          type = "error"
        )
        # Create results folder structure for metadata + count matrix inputs and return metadata dataframe reactive
      } else{
        dir.create(file.path(local_results_folder(), "/timeor/"),
                   showWarnings = FALSE)
        dir.create(file.path(
          paste0(local_results_folder(), "/timeor/"),
          "results/"
        ), showWarnings = FALSE)
        dir.create(file.path(
          paste0(local_results_folder(), "/timeor/results/"),
          "analysis/"
        ), showWarnings = FALSE)
        return(metadata_df())
      }
      
      # If SraRunTable used
    } else{
      print("SraRunTable toggle")
      sra_input(TRUE)
      
      # Check format of SraRunTable file
      if (!(any(grepl(
        "treatment", colnames(sra_metadata_df), ignore.case = TRUE
      )) &
      any(grepl(
        "time", colnames(sra_metadata_df), ignore.case = TRUE
      )) &
      any(grepl(
        "replicate", colnames(sra_metadata_df), ignore.case = TRUE
      )) &
      any(grepl(
        "batch", colnames(sra_metadata_df), ignore.case = TRUE
      )) &
      any(grepl(
        "Run", colnames(sra_metadata_df), ignore.case = TRUE
      )))) {
        shinyalert(
          "Please check SraRunTable.",
          "1) Toggle 'Input file' to 'SraRunTable' in question 4 above.
                   2) Make sure SraRunTable file includes these columns:
                   'run', 'treatment', 'replicate', 'time', batch'.",
          type = "error"
        )
        
        # Parse SraRunTable into metadata, create timeor results directory structure, and return metadata dataframe reactive
      } else{
        df <- runParseMetaData(sra_metadata_df)
        metadata_df(df)
        return(metadata_df())
      }
    }
  }
  
  runParseMetaData <- function(updateProgress = NULL) {
    print("runParseMetadata")
    # Update progress bar
    if (is.function(updateProgress)) {
      text <- "Parsing Meta Data"
      updateProgress(detail = text)
    }
    
    # Call parse_sraRunTable script
    folder <-
      paste(local_results_folder(), "/timeor/data/", sep = "")
    SraRunTable_file <- metadata_filepath()
    script <-
      paste(app_dir, "/scripts/parse_sraRunTable.py", sep = "")
    command <-
      paste("python",
            script,
            local_results_folder(),
            SraRunTable_file,
            sep = " ")
    system(command, intern = TRUE)
    
    # Update filepath of converted SraRunTable to Metadata
    metadata_filepath(paste(folder, "metadata.csv", sep = "/"))
    print("metadata_filepath")
    print(metadata_filepath())
    
    # Update dataframe to show in UI
    sra_converted_to_metadata_df <-
      as.data.frame(read.csv(
        file = metadata_filepath(),
        sep = ",",
        header = TRUE
      ))
    
    # Create dataframe from converted SraRunTable into metadata file
    showNotification("Metadata has been parsed and folder structure has been created.")
    return(sra_converted_to_metadata_df)
  }
  
  # Output metadata table
  output_metadata_df <- function(metad_df) {
    output$metadataTable <- DT::renderDataTable({
      metad_df
    }, options = (list(
      pageLength = 5, scrollX = TRUE
    )))
  }
  
  # Render metadata table from USER INPUT data
  observeEvent(input$metadataFile, {
    if (!is.null(input$metadataFile)) {
      # Set user input data reactive to True
      user_input_data(TRUE)
      
      # Create metadata dataframe
      df <- create_metadata_df()
      
      # Output metadata dataframe
      output_metadata_df(df)
    }
  })
  
  ## "Determine Adaptive Default Parameters"
  
  ## Adaptive Defaults ########
  observeEvent(input$run_adaptive_defaults, {
    # Make sure all questions have been answered.
    if ((input$organism == "NA") |
        (input$sequencing == "NA") |
        (input$experiment == "NA")) {
      shinyalert(
        "Please answer questions.",
        "Please load data and then select the organism, sequencing and experiment type.",
        type = "error"
      )
    } else{
      # Output statement about data processing
      if (!is.null(metadata_filepath())) {
        output$dataProcessType <- renderText({
          if (metadata_input()) {
            paste0(
              "<p> <br>",
              "You uploaded a metadata file. No need to retrieve and process raw data. <br> Please proceed to 'Process Count Matrix' tab.",
              "</p>"
            )
          } else if (sra_input()) {
            paste0(
              "<p> <br>",
              "Raw data will be retrieved, quality checked, and aligned. 
              Choose the alignment method to then generate the read count matrix (below).",
              "</p>"
            )
          } else{
            shinyalert(
              "Please follow 'Search and Retrieve' steps 4 and 5).",
              "Please indicate (step #4) and upload (step #5) a Metadata or SraRunTable file.",
              type = "error"
            )
          }
        })
      } else{
        shinyalert(
          "Please follow 'Search and Retrieve' steps 4 and 5).",
          "Please indicate (step #4) and upload (step #5) a Metadata or SraRunTable file.",
          type = "error"
        )
      }
    }
  })
  
  # Get FastQ files
  getfastQFiles <- function(folder, updateProgress = NULL) {

    print("Running command to get .fastq files")
    
    script <-
      paste(app_dir, "/scripts/get_fastq_files.sh", sep = "")
    accessionList <-
      paste(local_results_folder(),
            "/timeor/data/SraAccList.csv",
            sep = "")
    command <- paste(script, folder, accessionList)
    system(command, intern = TRUE)
    showNotification("All .fastq files saved in personal analysis session folder.")
    
    retrieveDone(TRUE)
  }
  
  # Performing quality control
  qualityControl <- function(qc_results_folder, updateProgress = NULL) {
      req(retrieveDone())
      print("Running command to perform QC.")
      
      script <-  paste(app_dir, "/scripts/run_fastQC.sh", sep = "")
      fastQDir <-
        paste(local_results_folder(), "/timeor/data/fastq", sep = "")
      command <- paste(script, fastQDir, qc_results_folder)
      system(command, intern = TRUE)
      showNotification("Quality control results saved in personal analysis session folder.")
      
      qcDone(TRUE)
    }
  
  # Output MultiQC download html link
  output$qcDownload <- downloadHandler(
    filename = function() {
      paste("multiqc", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      qc_folder <-
        paste(local_results_folder(),
              "/timeor/results/preprocess/fastqc",
              sep = "")
      qcHTML_file <-
        Sys.glob(file.path(qc_folder, "/multiqc_report.html"))
      
      req(file.exists(qcHTML_file))
      writeLines(readLines(qcHTML_file), file)
    }
  )
  
  # Moving .fastq files into organized folders for alignment and further analyses
  fastq_files_into_folders <- function() {
    print("Running command to move .fastq files into folders")
    
    script <-
      paste("python ", app_dir, "/scripts/fastq_into_folders.py", sep = "")
    fastQDir <-
      paste(local_results_folder(), "/timeor/data/fastq/", sep = "")
    command <- paste(script, fastQDir, metadata_filepath())
    system(command, intern = TRUE)
    
    fastqOrgDone(TRUE)
  }
  
  # Run HISAT2
  alignHISAT2 <-
    function(results_dir, local_dir, seq) {
      #, updateProgress = NULL){
      print("Running command to align .fastq files using HISAT2")
      
      script <- paste(app_dir, "/scripts/run_HISAT2.sh", sep = "")
      fastQDir <- paste(local_dir, "/timeor/data/fastq/", sep = "")
      
      paired <- 1
      if (seq == "se") {
        paired <- 0
      }
      
      command <-
        paste(script,
              fastQDir,
              results_dir,
              paired,
              app_dir,
              input$organism)
      system(command, intern = TRUE)
      
      # Return true for 'future' command
      TRUE
    }
  
  # Plot HISAT2 alignment
  plotAlignHISAT2 <- function(bamFolder, seq) {
    print("Running command to plot HISAT2 alignment")
    
    paired_or_not <- 1
    if (seq == "se") {
      paired_or_not <- 0
    }
    #     compareMultiple <- 1
    #     if (input$multipleMethods == "No"){
    #       compareMultiple <- 0
    #     }
    
    script <-
      paste("python ", app_dir, "/scripts/plot_alignment.py", sep = "")
    command <-
      paste(script, 0, bamFolder, paired_or_not, "HISAT2", "NA")
    system(command, intern = TRUE)
    
    # Return true for 'future' command - to show plot immediately
    TRUE
  }
  
  # Run Bowtie2
  alignBowtie2 <-
    function(results_dir, local_dir, seq) {
      #updateProgress = NULL){
      print("Running command to align .fastq files using Bowtie2")
      
      script <- paste(app_dir, "/scripts/run_Bowtie2.sh", sep = "")
      fastQDir <-
        paste(local_results_folder(), "/timeor/data/fastq/", sep = "")
      
      paired <- 1
      if (seq == "se") {
        paired <- 0
      }
      
      command <-
        paste(script,
              fastQDir,
              results_dir,
              paired,
              app_dir,
              input$organism)
      system(command, intern = TRUE)
      
      # Return true for 'future' command - to show results immediately
      TRUE
    }
  
  # Plot Bowtie2 alignment
  plotAlignBowtie2 <- function(bamFolder, seq) {
    print("Running command to plot Bowtie2 alignment")
    
    paired_or_not <- 1
    if (seq == "se") {
      paired_or_not <- 0
    }
    
    script <-
      paste("python ", app_dir, "/scripts/plot_alignment.py", sep = "")
    command <-
      paste(script, 0, bamFolder, paired_or_not, "Bowtie2", "NA")
    system(command, intern = TRUE)
    
    # Return true for 'future' command - to show results immediately
    TRUE
  }
  
  # Run HTSeq
  genCountMatrix <- function(align_folder, htseq_folder) {
    print("aaRunning command to generate count matrix from chosen alignment method")
    print(align_folder)
    print(organism_ad())
    
    # Count matrices created for each sample
    script <- paste(app_dir, "scripts/run_HTSeq.sh", sep = "/")
    command <- paste(script, align_folder, organism_ad(), app_dir)
    
    system(command, intern = TRUE)
    
    
    # All count matrices merged for all samples
    script2 <-
      paste("python ", app_dir, "/scripts/htseq_merge.py", sep = "")
    command2 <- paste(script2, htseq_folder)
    system(command2, intern = TRUE)
    
    # Return true for 'future' command - to return results immediately
    TRUE
  }
  
  # If raw data are to be processed, get .fastq files
  observeEvent(run_process_raw(), {
    req(sra_input())
    req(input$run_adaptive_defaults)
    
    # Checking for .fastq files
    folder_fastQ <-
      paste(local_results_folder(), "/timeor/data/fastq/", sep = "")
    if (length(list.files(folder_fastQ, pattern = ".fastq")) == 0) {
      print("Calling fastq-dump")
      getfastQFiles(folder_fastQ)
    } else{
      print("Already performed fastq-dump")
      retrieveDone(TRUE) # .fastq files already retrieved
    }
  })
  
  # Show checkmark retrieved .fastq files finished
  output$textCheckMark_raw <- renderText({
    req(retrieveDone())
    session$sendCustomMessage(type = 'print',
                              message = list(selector = 'textCheckMark_raw', html = "✓"))
    return("✓")
  })
  
  # Checking for FastQC and MultiQC files
  observeEvent(retrieveDone(), {
    req(retrieveDone())
    results_folder_qc <-
      paste(local_results_folder(),
            "/timeor/results/preprocess/fastqc/",
            sep = "")
    if (length(list.dirs(results_folder_qc)) == 1) {
      print("Calling QC script")
      qualityControl(results_folder_qc)
    } else{
      print("Already performed QC")
      qcDone(TRUE) # QC results already produced
    }
    
    # Output QC table
    output$qcTable <- renderDataTable({
      req(!sim_demo_data())
      qc_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/preprocess/fastqc/multiqc_data/",
          sep = ""
        )
      qcTable_file <-
        Sys.glob(file.path(qc_folder, "/multiqc_general_stats.txt"))
      
      d <- read.table(qcTable_file, header = T)
      names(d) <-
        names(d) %>% str_replace("FastQC_mqc.generalstats.fastqc.", "")
      return(as.data.frame(d))
      
    }, options = (list(
      pageLength = 5,
      scrollX = TRUE,
      autoWidth = FALSE
    )))
    
    #Move .fastq files into folders for processing
    fastq_dir <-
      paste(local_results_folder(), "/timeor/data/fastq/", sep = "")
    if (length(list.dirs(fastq_dir)) == 1) {
      print("Calling script to organize .fastq files")
      fastq_files_into_folders()
    } else{
      print("Already organized .fastq files")
      fastqOrgDone(TRUE)
    }
    
    # Show checkmark for quality control finished
    output$textCheckMark_qc <- renderText({
      req(qcDone())
      session$sendCustomMessage(type = 'print',
                                message = list(selector = 'textCheckMark_qc', html = "✓"))
      return("✓")
    })
  })
  
  # Enable download MultiQC results
  observeEvent(qcDone(), {
    if (qcDone() == TRUE) {
      toggleState("qcDownload", condition = T)
    } else{
      toggleState("qcDownload", condition = F)
    }
  })
  
  # Checking to align using HISAT2
  observeEvent(qcDone(), {
    req(qcDone())
    results_folder_HISAT2 <-
      paste(local_results_folder(),
            "/timeor/results/preprocess/alignment/hisat2/",
            sep = "")
    local_dir <- local_results_folder()
    seq <- input$sequencing
    
    # Run HISAT2 if needed
    if (length(list.files(
      results_folder_HISAT2,
      recursive = "TRUE",
      pattern = "bam"
    )) == 0) {
      print("Calling align HISAT2")
      future({
        alignHISAT2(results_folder_HISAT2, local_dir, seq)
      }) %...>% alignHdone()
    } else{
      print("Already aligned using HISAT2")
      alignHdone(TRUE)
    }
    
    # Plot HISAT2 if needed
    req(alignHdone())
    align_file <-
      paste(results_folder_HISAT2,
            "/HISAT2_plot_alignment.svg",
            sep = "")
    if (!file.exists(align_file)) {
      print("Calling plot HISAT2 alignments")
      future({
        plotAlignHISAT2(results_folder_HISAT2, seq)
      }) %...>% plotAlignH()
    } else{
      print("Already plotted HISAT2 alignments")
      plotAlignH(TRUE)
    }
    
    output$alignment_graph_HISAT2 <- renderImage({
      req(plotAlignH())
      
      # Get plot
      align_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/preprocess/alignment/hisat2/",
          sep = ""
        )
      align_file <-
        paste(align_folder, "/HISAT2_plot_alignment.svg", sep = "")
      
      req(file.exists(align_file))
      list(
        src = align_file,
        style = "object-fit: contain; max-height: 100%; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
        contentType = 'image/svg+xml',
        alt = "No alignment to show"
      )
    }, deleteFile = FALSE)
    
  })
  
  # Checking to align using Bowtie2
  observeEvent(run_Bowtie(), {
    req(qcDone())
    req(alignHdone())
    
    results_folder_Bowtie2 <-
      paste(
        local_results_folder(),
        "/timeor/results/preprocess/alignment/bowtie2/",
        sep = ""
      )
    local_dir <- local_results_folder()
    seq <- input$sequencing
    
    # Run Bowtie2 if needed
    if (length(list.files(
      results_folder_Bowtie2,
      recursive = "TRUE",
      pattern = "bam"
    )) == 0) {
      print("Calling align Bowtie2")
      future({
        alignBowtie2(results_folder_Bowtie2, local_dir, seq)
      }) %...>% alignBdone()
    } else{
      print("Already aligned using Bowtie2")
      alignBdone(TRUE)
    }
    
    # Plot Bowtie2 if needed
    req(alignBdone())
    align_file <-
      paste(results_folder_Bowtie2,
            "/Bowtie2_plot_alignment.svg",
            sep = "")
    if (!file.exists(align_file)) {
      print("Calling plot Bowtie2 alignments")
      future({
        plotAlignBowtie2(results_folder_Bowtie2, seq)
      }) %...>% plotAlignB()
    } else{
      print("Already plotted Bowtie2 alignments")
      plotAlignB(TRUE)
    }
    
    output$alignment_graph_Bowtie2 <- renderImage({
      req(plotAlignB())
      
      # Get plot
      align_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/preprocess/alignment/bowtie2/",
          sep = ""
        )
      align_file <-
        paste(file.path(align_folder, "/Bowtie2_plot_alignment.svg"))
      
      req(file.exists(align_file))
      
      list(
        src = align_file,
        style = "object-fit: contain; max-height: 100%; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
        contentType = 'image/svg+xml',
        alt = "No alignment to show"
      )
    }, deleteFile = FALSE)
  })
  
  # Show checkmark alignment finished
  output$textCheckMark_align <- renderText({
    req(alignHdone())
    req(alignBdone())
    session$sendCustomMessage(
      type = 'print',
      message = list(selector = 'textCheckMark_align', html = "✓")
    )
    return("✓")
  })
  
  # Choosing an aligner (HISAT2 or Bowtie2)
  observeEvent(whichAligner(), {
    if (input$aligner == 'bowtie2') {
      whichAligner('bowtie2')
    }
  })
  
  # Checking to create count matrix using HTSeq
  observeEvent(alignBdone(), {
    req(!sim_demo_data())
    if (alignBdone() == FALSE) {
      toggleState("runGenCountMat", condition = F)
    } else{
      # Enable button to generate count matrix
      toggleState("runGenCountMat", condition = T)
    }
  })
  
  # Checking to run HTSeq only once 'Generate Count Matrix' button is enabled
  observeEvent(input$runGenCountMat, {
    req(input$runGenCountMat)
    shinyalert("Quick start:",
                  "You completed the 'Process Raw Data' tab demo.
                  
                  TIMEOR accepts 2 input types: 
                  (1) raw .fastq files
                  (2) read count matrix. 
              
                  For (2) in 'Example Data' (side-bar) under

                  'Load count matrix' click 'Metadata & read count file' button. 
                  
                  This will guide you through the rest of the full method demo.",
                  type = "info"
                )

    # 'Generate count matrix' button must be pressed
    results_folder_alignment <-
      paste(
        local_results_folder(),
        "/timeor/results/preprocess/alignment/",
        whichAligner(),
        sep = "/"
      )
    results_folder_htseq <-
      paste(
        local_results_folder(),
        "/timeor/results/preprocess/count_matrix/htseq/",
        sep = ""
      )
    
    if (length(list.files(
      results_folder_htseq,
      recursive = "TRUE",
      pattern = "htseq"
    )) == 0) {
      print("Calling count matrix script")
      future({
        genCountMatrix(results_folder_alignment, results_folder_htseq)
      }) %...>% countMatDone()
    } else{
      print("Already created a count matrix")
      countMatDone(TRUE)
    }
    
    # Show checkmark count matrix finished
    output$textCheckMark_count <- renderText({
      req(countMatDone())
      session$sendCustomMessage(
        type = 'print',
        message = list(selector = 'textCheckMark_count', html = "✓")
      )
      return("✓")
    })
  })
  
  #################### Pre-Process Stage ##################
  #################### Load Count Matrix ##################
  
  # Reading gene x sample input matrix from .csv
  data_count_matrix <- reactive ({
    print("countMatrix_filepath()")
    print(countMatrix_filepath())
    if (input$runSim) {
      countMatrix_filepath(paste(
        local_results_folder(),
        "/timeor/data/countMatrix.csv",
        sep = ""
      ))
      print("countMatrix:")
      print(countMatrix_filepath())
      countMatrix_df(
        read.table(
          file = countMatrix_filepath(),
          sep = input$sep,
          header = input$header,
          dec = ",",
          row.names = "ID"
        )
      )
    } else{
      print("else - not using simulation")
      req(input$count_mat)
      print("DATAPATH")
      print(input$count_mat$datapath)
      print(input$count_mat)
      
      countMatrix_filepath(input$count_mat$datapath)
      count_df <-
        read.table(
          file = countMatrix_filepath(),
          sep = input$sep,
          header = input$header,
          dec = ","
        )
      if (!(any(grepl("ID", colnames(count_df))))) {
        shinyalert(
          "Please check Count Matrix.",
          "1) Make sure 1st column is unique identifiers (IDs), and says 'ID'.
           2) Make sure metadata sample names match other column names",
          type = "error"
        )
        
        # Make ID the rownames and remove the column ID. Save df to reactive countMatrix_df
      } else{
        rownames(count_df) <- count_df$ID
        countMatrix_df(count_df[, !(names(count_df) %in% "ID")])
      }
    }
    if (is.null(countMatrix_df())) {
      return()
    }
    print("countMatrix_df")
    print(head(countMatrix_df()))
    # return count matrix dataframe reactive
    countMatrix_df()
  })
  
  # Saving .csv as a data frame
  df_count_matrix <- function() {
    return(as.data.frame(data_count_matrix()))
  }
  
  # Render countMatrix file as a DataTable
  output$countMatrix <- DT::renderDataTable({
    df_count_matrix()
  }, options = (list(
    pageLength = 5, scrollX = TRUE
  )))
  
  # Rendering heatmap plot
  output$heatmap <- renderPlotly({
    req(data_count_matrix())
    p <- heatmaply(
      df_count_matrix(),
      xlab = "",
      ylab = "",
      main = "",
      scale = "column",
      margins = c(60, 100, 40, 20),
      grid_color = "white",
      grid_width = 0.00001,
      hide_colorbar = FALSE,
      branches_lwd = 0.1,
      fontsize_row = 5,
      fontsize_col = 5,
      labCol = colnames(df_count_matrix()),
      labRow = rownames(df_count_matrix()),
      heatmap_layers = theme(axis.line = element_blank())
    )
    p
  })
  
  # PCA Scatter Before
  output$pcaScatBefore <- renderPlotly({
    req(data_count_matrix())
    dataSubset <- data_count_matrix() %>%
      dplyr::select(-starts_with("FlyBaseID"))
    p <-
      autoplotly(prcomp(dataSubset), data = data_count_matrix(), frame = FALSE)
    p
  })
  
  # PCA loadings bar plot before normalization
  output$pcaBarBefore <- renderPlot({
    req(data_count_matrix())
    fviz_eig(prcomp(data_count_matrix()),
             addlabels = TRUE,
             ylim = c(0, 100))
  })
  
  # Rendering correlation plot
  output$correlationBefore <- renderPlotly({
    req(data_count_matrix())
    req(input$correMethodBefore)
    val <- input$correMethodBefore
    if (input$correMethodBefore == "Spearman") {
      val <- "spearman"
    } else {
      val <- "pearson"
    }
    p <- heatmaply(
      cor(df_count_matrix(), method = c(val)),
      xlab = "Experiments",
      ylab = "Experiments",
      main = paste(
        input$correMethodBefore,
        "Correlation Between Experiments"
      ),
      margins = c(40, 40),
      limits = c(-1, 1),
      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red")
    )
  })
  
  #################### Pre-Process Stage ##################
  ############### Normalize and Correct Data ##############
  
  # Returning the normalized data as a global variable
  outputNormalizedData <- reactive({
    req(!is.null(countMatrix_df()))
    req(input$runCorrection)
    print("outputNoralizedData")
    
    shinyalert(
      "Completed Pre-processing",
      "Proceed to Primary Analysis (side-bar). 
      Fill in grey box and follow pop-ups.",
      type = "info"
    )

    # Computing normalization
    if (input$normMethods == "Trimmed Mean of M-Values") {
      print("TMM")
      normData <<-
        as.data.frame(heatmaply::normalize(countMatrix_df(), method = "tmm", trim = 0.3)) # Making normData a global variable for further access
      print("normData")
    } else if (input$normMethods == "Upper Quartile") {
      normData <<- as.data.frame(upperNormalization(countMatrix_df()))
    }
  })
  
  # Returning the corrected data as a global variable
  outputCorrectedData <- reactive({
    req(outputNormalizedData())
    req(input$runCorrection)
    
    # Saving corrected read count matrix
    correctedData <<-
      as.data.frame(reconstructData(harmanCorrection(normData, metadata_df()), this = "corrected"))
    print("correctedData")
    
    correctedData <-
      cbind(ID = rownames(correctedData), correctedData)
    rownames(correctedData) <- 1:nrow(correctedData)
    norm_corr_countMatrix_filepath(
      paste(
        local_results_folder(),
        "/timeor/results/analysis/countMatrix_norm_corrected.csv",
        sep = '/'
      )
    )
    write.csv(
      correctedData,
      norm_corr_countMatrix_filepath(),
      row.names = FALSE,
      quote = FALSE
    )
    print("norm_corr_countMatrix_filepath()")
    print(norm_corr_countMatrix_filepath())
  })
  
  # PCA Scatter after norm/correction
  output$pcaScatAfter <- renderPlotly({
    req(input$runCorrection)
    req(outputCorrectedData())
    dataSubset <- correctedData %>%
      dplyr::select(-starts_with("ID"))
    print("p")
    p <-
      autoplotly(prcomp(dataSubset), data = correctedData, frame = FALSE)
    p
  })
  
  # PCA Bar after normalization and correction
  output$pcaBarAfter <- renderPlot({
    req(input$runCorrection)
    req(outputCorrectedData())
    fviz_eig(
      prcomp(correctedData),
      title = "",
      addlabels = TRUE,
      ylim = c(0, 100)
    )
    
  })
  
  # Output correlation plot after normalization and correction
  output$correlationAfter <- renderPlotly({
    req(input$runCorrection)
    req(input$correMethodAfter)
    req(outputCorrectedData())
    val <- input$correMethodAfter
    if (input$correMethodBefore == "Spearman") {
      val <- "spearman"
    } else {
      val <- "pearson"
    }
    p <- heatmaply(
      cor(correctedData, method = c(val)),
      xlab = "Experiments",
      ylab = "Experiments",
      main = paste(input$correMethodAfter, "Correlation Between Experiments"),
      margins = c(40, 40),
      limits = c(-1, 1),
      scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red")
    )
  })
  
  #################### Primary Analysis ##################
  
  output$deParameters <- renderText({
    if (input$multipleMethods == "Yes") {
      paste(
        "You selected to compare multiple methods. TIMEOR will run ImpulseDE2, Next maSigPro and DESeq2."
      )
    } else if (input$distant) {
      paste(
        "You selected to not compare methods and have distant timepoints. TIMEOR will run DESeq2."
      )
    } else if (input$ct || input$ctstc) {
      paste(
        "You selected to not compare methods and have close timepoints. TIMEOR will run ImpulseDE2."
      )
    }
  })
  
  # Run ImpulseDE2
  runImpulseDE2 <- function() {
    impulsede2_script <-
      paste("Rscript ", app_dir, "/scripts/ImpulseDE2.r", sep = "")
    path_to_output <-
      paste(local_results_folder(),
            "/timeor/results/analysis/",
            sep = "")
    
    # Run ImpulseDE2.r script
    print("Running ImpulseDE2 command")
    command <-
      paste(
        impulsede2_script,
        metadata_filepath(),
        countMatrix_filepath(),
        path_to_output,
        analysis_folder_name(),
        adj_p_analysis_folder()
      )
    print(command)
    system(command, intern = TRUE)
  }
  
  # Run NextMaSigPro 
  runNextMaSigPro <- function() {
    nextMaSigPro_script <-
      paste("Rscript ", app_dir, "/scripts/next_maSigPro.r", sep = "")
    path_to_output <-
      paste(local_results_folder(),
            "/timeor/results/analysis/",
            sep = "")
    
    print("Normalized and corrected countMatrix filepath")
    print(norm_corr_countMatrix_filepath())
    
    # Run next_maSigPro.r script
    print("Running NextMaSigPro command")
    command <-
      paste(
        nextMaSigPro_script,
        metadata_filepath(),
        norm_corr_countMatrix_filepath(),
        path_to_output,
        analysis_folder_name(),
        adj_p_analysis_folder()
      )
    print(command)
    system(command, intern = TRUE)
  }
  
  # Run DESeq2
  runDESeq2 <- function() {
    path_to_output <-
      paste(local_results_folder(),
            "/timeor/results/analysis/",
            sep = "")
    batch_effect <- 1 # yes to consider batch effects
    time <- 1 # yes to consider time course (LRT test)
    
    print("Running DESeq2 command") # imported from "source" above
    print(
      paste(
        metadata_filepath(),
        countMatrix_filepath(),
        path_to_output,
        analysis_folder_name(),
        batch_effect,
        time,
        adj_p_analysis_folder(),
        sep = ","
      )
    )
    
    # Command
    run_DESeq2(
      metadata_filepath(),
      countMatrix_filepath(),
      path_to_output,
      analysis_folder_name(),
      batch_effect,
      time,
      adj_p_analysis_folder(),
      'control', 
      input$organism
    )
  }
  
  # Run all differential expression methods
  run_all_DE <- function() {
    # Requirements
    req(!file.exists(impulseDE2File()))
    req(!is.null(norm_corr_countMatrix_filepath()))
    req(!file.exists(nextMaSigProFile()))
    req(!file.exists(deSeqFile()))
    req(input$whichMethodInput)
    req(analysis_folder_name())
    req(input$organism != "NA")
    req(input$runDE) # Run button
    
    print("Run ImpulseDE2")
    runImpulseDE2()
    print("Run NextMaSigPro")
    runNextMaSigPro()
    print("Run DESeq2")
    runDESeq2()
    print("Completed running of ImpulseDE2, NextMaSigPro, and DESeq2.")
  }
  
  # Run differential expression (DE) method(s) depending on user input
  determine_DE_methods <- observe({
    # Requirements
    req(input$runDE) # Run button
    req(!is.null(countMatrix_filepath()))
    req(!is.null(metadata_filepath()))
    
    # Set analysis_folder_name if needed (i.e. not preset by clicking on demo dataset buttons
    if (is.null(analysis_folder_name()) & !sim_demo_data()) {
      analysis_folder_name(input$analysisFolder)
    }
    
    # If p-value has not been set from simulation or real data (demos)
    if (is.null(adj_p_analysis_folder()) & !sim_demo_data()) {
      adj_p_analysis_folder(as.numeric(input$pValue))
    }
    
    # Define all DE method filepaths
    # ImpulseDE2 outputs
    impulseDE2File(
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/impulsede2/impulsede2_clustermapInput_padj",
        adj_p_analysis_folder(),
        ".csv",
        sep = ""
      )
    )
    impulseDE2_out_file(
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/impulsede2/impulsede2_output_padj",
        adj_p_analysis_folder(),
        ".csv",
        sep = ""
      )
    )
    
    # NextMaSigPro output
    nextMaSigProFile(
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/nextmasigpro/nextMaSigPro_clustermapInput_padj",
        adj_p_analysis_folder(),
        ".csv",
        sep = ""
      )
    )
    
    # DESeq2 outputs
    deSeqFile(
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/deseq2/deseq2_noShrinkage_clustermapInput_padj",
        adj_p_analysis_folder(),
        ".csv",
        sep = ""
      )
    )
    deSeq_out_file(
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/deseq2/deseq2_output_noShrinkage_padj",
        adj_p_analysis_folder(),
        ".csv",
        sep = ""
      )
    )
    
    # Determine which DE method(s) to run
    if (input$multipleMethods == "Yes") {
      if (!file.exists(impulseDE2File()) ||
          !file.exists(nextMaSigProFile()) || !file.exists(deSeqFile())) {
        # Run all DE methods
        run_all_DE()
      }
    } else if (input$distant) {
      if (!file.exists(deSeqFile())) {
        # Run DESeq2 (categorical DE method)
        runDESeq2()
      }
    } else if (input$ct || input$ctstc) {
      if (!file.exists(impulseDE2File())) {
        # Run ImpulseDE2
        runImpulseDE2()
      }
    }
    
    # Notify user to use Venn diagram to compare differential expression results.
    showNotification("Click 'Render Venn Diagram' to compare differential expression methods' results.")
    diffExpDone(1)
  })
  
  # Create venn diagram
  run_venn <- function() {
    # Get previous study if provided
    path_to_prev_study <- ""
    if (!is.null(input$prevStudy)) {
      path_to_prev_study <- input$prevStudy$datapath
    } else {
      path_to_prev_study = "NA"
    }
    cat("path_to_prev_study", path_to_prev_study)
    past_study_name <- "Previous_Study"
    
    # Call Venn diagram script using Intervene
    output_dir <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/www/",
        sep = ""
      )
    script <-
      paste("Rscript ", app_dir, "/scripts/venns_intervene.r", sep = "")
    command <-
      paste(
        script,
        impulseDE2File(),
        deSeqFile(),
        nextMaSigProFile(),
        path_to_prev_study,
        past_study_name,
        output_dir
      )
    system(command, intern = TRUE)
  }
  
  # Output venn diagram
  output$vennMethods <- renderImage({
    # Requirements
    req(input$renderVen)
    req(!is.null(analysis_folder_name()))
    
    # Call to create Venn diagram
    run_venn()
    
    # Render png
    list(
      src = paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/www/Intervene_venn.png",
        sep = ""
      ),
      height = "100%",
      width = "auto",
      style = "display: block; margin-left: auto; margin-right: auto;",
      # style centers the image
      contentType = 'image/png/pdf',
      alt = "No data found"
    )
  }, deleteFile = FALSE)

  # Notify user to proceed to Secondary Analysis tab
  observeEvent(input$renderVen, {
    req(input$renderVen)
    shinyalert(
      "Completed Primary Analysis",
      "NOTE: demo chose 'ImpulseDE2' output and 
      'automatic' gene trajectory clustering. 
      On new data the user can choose these.
      
      Proceed to Secondary Analysis (side-bar). 
      
      Follow pop-ups and fill in grey boxes.",
      type = "info"
    )
  })
  
  # Output results of differential expression as a data fame
  output$deResults <- DT::renderDataTable({
    # Requirements
    print("analysis_folder_name()")
    print(analysis_folder_name())
    print(adj_p_analysis_folder())
    print(input$runDE)
    req(diffExpDone())
    req(!is.null(analysis_folder_name()))
    req(!is.null(adj_p_analysis_folder()))
    req(input$runDE)
    
    
    if (input$whichMethodInput == "ImpulseDE2") {
      deResult <<- read.csv(impulseDE2_out_file(), header = TRUE)
    } else if (input$whichMethodInput == "NextMaSigPro") {
      deResult <<- read.csv(nextMaSigProFile(), sep = ",", header = TRUE)
    } else{
      deResult <<- read.csv(deSeq_out_file(), sep = ",", header = TRUE)
    }
    deResult
  }, options = (list(
    pageLength = 5, scrollX = TRUE
  )))
  
  # Output clustermap
  output$deClustering <-
    output$deClustering1 <-
    output$deClustering2 <- output$deClustering3 <- renderPlotly({
      # Requirements
      req(analysis_folder_name())
      req(adj_p_analysis_folder())
      req(input$runDE)
      req(input$numClusts)
      
      # Setting output results folder
      output_dir <-
        paste(local_results_folder(),
              "/timeor/results/analysis/",
              sep = "")
      desired_output_name <- input$resultsFolder
      
      # Determine which method (from bottom left of Primary Analysis) to display (on bottom right of Primary Analysis)
      if (input$whichMethodInput == "ImpulseDE2") {
        clustermap <<- impulseDE2File()
      } else if (input$whichMethodInput == "NextMaSigPro") {
        clustermap <<- nextMaSigProFile()
      } else if (input$whichMethodInput == "DESeq2") {
        clustermap <<- deSeqFile()
      }
      
      # Setting parameters to call TIMEOR's clustermap script
      timepoint <- input$timeSeries
      if (timepoint == "ct" || timepoint == "ctstc") {
        timepoint <- 1
      } else {
        timepoint <- 0
      }
      if (input$numClusts == "automatic") {
        numClusters("0")
      } else{
        numClusters(input$numClusts)
      }
      distance <- "euclidean"
      clusterMethod <- "ward.D2"
      
      #updateTextInput(session, "numClusts", value = paste("dme"))
      
      # Calling clustermap script
      produceClusterMap(
        output_dir,
        clustermap,
        analysis_folder_name(),
        timepoint,
        numClusters(),
        distance,
        clusterMethod
      )
      
    })
  
  #################### Secondary Analysis ##################
  ######################## Enrichment ######################
  
  # Analysis button pressed (toggle between on and off)
  observeEvent(input$runEnrichment, {
    req(!real_demo_data())
    if (input$runEnrichment == TRUE) {
      toggleState("inputNumClust", condition = F)
      toggleState("memeHTML", condition = T)
    } else{
      toggleState("inputNumClust", condition = T)
      toggleState("memeHTML", condition = F)
    }
  })
  
  # Print gene list for user input cluster
  output$textWithHTML <- renderText({
    print(input$inputNumClust)
    req(!is.null(analysis_folder_name()))
    req(input$inputNumClust != "NA")
    
    currentClust <- input$inputNumClust
    parent_dir <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/clusters/",
        currentClust,
        "/",
        sep = ""
      )
    geneList_file <-
      Sys.glob(file.path(parent_dir, "*geneList*")) # get the geneList file
    print("geneList_file")
    print(geneList_file)
    # print(read.csv(geneList_file, skip=1))
    geneList <-
      as.character(read.table(geneList_file, skip = 1, sep = ",")$V1) # skip the header
    geneList_comma_sep <-
      toString(paste0("'", geneList, "'")) # add commas and quotes to make it easy to copy and paste the list elsewhere
    clus_num_to_colr <-
      read.table(
        paste(
          dirname(parent_dir),
          "cluster_color_to_number_map.txt",
          sep = "/"
        ),
        sep = "\r",
        comment.char = "-"
      )
    print(as.character(clus_num_to_colr[input$inputNumClust, 1]))
    #paste("<span style='color:red\'>This is red text</span>")
    print(
      paste(
        "<span style=\'color:",
        "\\'",
        as.character(clus_num_to_colr[input$inputNumClust, 1]),
        "\'>\'",
        geneList_comma_sep,
        "\'</span>",
        sep = ""
      )
    )
    return(
      paste(
        "<span style='color:",
        as.character(clus_num_to_colr[input$inputNumClust, 1]),
        "\'>\'",
        geneList_comma_sep,
        "\'</span>",
        sep = ""
      )
    )
  })
  
  # Run secondary analysis enrichment tab scripts
  check_to_run_secondary_enrichment <- observe({
    req(input$organism != "NA")
    req(input$inputNumClust != "NA")
    req(input$runEnrichment)
    
    # Identify folders to indicate that the cluster has been run (motif file will always be made even if nothing else)
    motif_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/clusters/",
        input$inputNumClust,
        "/MEME/",
        sep = ""
      )
    print("motif_folder")
    print(motif_folder)
    motif_file <-
      Sys.glob(file.path(motif_folder, "/*DNAseqs.fasta*")) # get the geneList file
    print("length ss-motif_file")
    print(motif_file)
    print(length(motif_file))
    
    # Run Enrichment if new individual cluster
    if (length(motif_file) == 0) {
      run_secondary_enrichment()
    }
  })
  
  run_secondary_enrichment <- function() {
    # Shared enrichment parameters
    condition <- paste(analysis_folder_name(), "_results/", sep = "")
    currentClust_dir <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        condition,
        "/clusters/",
        input$inputNumClust,
        "/",
        sep = ""
      )
    animal <- input$organism
    
    # GO and pathway script calls analysis_folder_name(),"_results/
    clusterProf_pathview_script <-
      paste("Rscript ",
            app_dir,
            "/scripts/clusterProfiler_pathview_indiv_cluster.r",
            sep = "")
    separate_timepoints <- input$timeSeries
    if (separate_timepoints == "ct" ||
        separate_timepoints == "ctstc") {
      separate_timepoints <- 0
    } else {
      separate_timepoints <- 1
    }
    p_value <- as.numeric(input$pValue)
    command_GO_path <-
      paste(
        clusterProf_pathview_script,
        currentClust_dir,
        analysis_folder_name(),
        separate_timepoints,
        animal,
        p_value
      )
    print(command_GO_path)
    system(command_GO_path, intern = TRUE)
    
    # Network script calls
    stringdb_script <-
      paste("Rscript ",
            app_dir,
            "/scripts/stringdb_indiv_cluster.r",
            sep = "")
    if (animal == "dme") {
      ncbi_id = 7227
    } else if (animal == "hse") {
      ncbi_id = 9606
    } else{
      ncbi_id = 10090
    }
    command_stringdb <-
      paste(stringdb_script, currentClust_dir, ncbi_id, app_dir)
    print(command_stringdb)
    system(command_stringdb, intern = TRUE)
    
    # De novo motif search script calls
    # preparing MEME data
    meme_prep_script <-
      paste(app_dir, "/scripts/meme_prep_indiv_cluster.py", sep = "")
    reformatted_gtf <-
      paste(app_dir,
            "/../genomes_info/dme/reformatted_genes_gtf.csv",
            sep = "")
    dm6_fa <- paste(app_dir, "/../genomes_info/dme/dm6.fa", sep = "")
    command_meme_prep <-
      paste(
        "python",
        meme_prep_script,
        reformatted_gtf,
        currentClust_dir,
        dm6_fa,
        0,
        animal,
        sep = " "
      )
    system(command_meme_prep, intern = TRUE)
    # running MEME
    meme_script <-
      paste(app_dir, "/scripts/run_meme_indiv_cluster.sh", sep = "")
    command_meme <- paste(meme_script, currentClust_dir)
    system(command_meme, intern = TRUE)
    showNotification("MEME finished.")
    
    print("Enrichment complete.")
  }
  
  # Molecular function box
  output$molecularFunc <- renderImage({
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    print("molecular function")
    GO_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        sep = ""
      )
    MF_folder <- paste(GO_folder, "/MF/", sep = "")
    MF_dotplot <-
      Sys.glob(file.path(MF_folder, "/www/*_dotplot_MF*")) # get the geneList file
    print("MF_dotplot")
    req(file.exists(MF_dotplot))
    print(MF_dotplot)
    list(
      src = paste(MF_dotplot),
      height = "100%",
      width = "100%",
      contentType = 'image/svg+xml',
      alt = "No data found"
    )
  }, deleteFile = FALSE)
  
  # Biological process box
  output$bioProc <- renderImage({
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    print("biological process")
    BP_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/BP/",
        sep = ""
      )
    BP_dotplot <-
      Sys.glob(file.path(BP_folder, "/www/*_dotplot_BP*")) # get the geneList file
    req(file.exists(BP_dotplot))
    list(
      src = paste(BP_dotplot),
      height = "100%",
      width = "100%",
      contentType = 'image/svg+xml',
      alt = "No data found"
    )
  }, deleteFile = FALSE)
  
  # Cellular component box
  output$cellComp <- renderImage({
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "na")
    print("cellular component")
    CC_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/CC/",
        sep = ""
      )
    CC_dotplot <-
      Sys.glob(file.path(CC_folder, "/www/*_dotplot_CC*")) # get the geneList file
    req(file.exists(CC_dotplot))
    list(
      src = paste(CC_dotplot),
      height = "100%",
      width = "100%",
      contentType = 'image/svg+xml',
      alt = "No data found"
    )
  }, deleteFile = FALSE)
  
  # Pathway plot
  output$pathway <- renderImage({
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    print("Pathway")
    pv_multi_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/www/",
        sep = ""
      )
    pv_multi_png <-
      Sys.glob(file.path(pv_multi_folder, "/*pathview.multi*png")) # get Pathview multi output png
    print("pv_multi_png")
    print(pv_multi_png)
    req(file.exists(pv_multi_png))
    
    list(
      src = pv_multi_png,
      style = "object-fit: contain; max-height: 100%; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Pathway does not exist for this cluster"
    )
  }, deleteFile = FALSE)
  
  # Network plot # aspect ratio and disable scrolling and remove white space
  output$network <- renderImage({
    # Set requirements
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    
    # Get png
    network_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/www/",
        sep = ""
      )
    network_png <-
      Sys.glob(file.path(network_folder, "/*stringdb_network*")) # get StringDB network path
    req(file.exists(network_png))
    
    # Format and get png
    list(
      src = network_png,
      style = "object-fit: contain; max-height: 100%; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
      #height = "auto", width = "auto",
      contentType = 'image/png',
      alt = "Pathway does not exist for this cluster"
    )
  }, deleteFile = FALSE)
  
  output$pval_network_enrichment <- renderText({
    # Set requirements
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    
    # Get p-value
    network_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/",
        sep = ""
      )
    network_pval_file <-
      file.path(network_folder, "/stringdb_info_table.tsv") # get StringDB network information path
    req(file.exists(network_pval_file))
    return(read.table(network_pval_file, header = T)$'pval_PPI_enrich'[1])
  })
  
  # Motif html link
  output$memeHTML <- downloadHandler(
    filename = function() {
      paste("meme", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      motif_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/analysis/",
          analysis_folder_name(),
          "_results/",
          "clusters/",
          input$inputNumClust,
          "/MEME/",
          sep = ""
        )
      motif_html <- Sys.glob(file.path(motif_folder, "/meme.html"))
      writeLines(readLines(motif_html), file)
    }
  )
  
  # Motif plot logo 1
  output$motif <- renderImage({
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    motif_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/MEME/",
        sep = ""
      )
    motif_file <- Sys.glob(file.path(motif_folder, "/logo1.png"))
    req(file.exists(motif_file))
    
    list(
      src = motif_file,
      style = "object-fit: contain; max-height: 95%; max-width: 95%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Enriched motif does not exist for this cluster."
    )
  }, deleteFile = FALSE)
  
  # Motif plot logo 2
  output$motif2 <- renderImage({
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    motif_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/MEME/",
        sep = ""
      )
    motif_file <- Sys.glob(file.path(motif_folder, "/logo2.png"))
    req(file.exists(motif_file))
    
    list(
      src = motif_file,
      style = "object-fit: contain; max-height: 95%; max-width: 95%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Enriched motif does not exist for this cluster."
    )
  }, deleteFile = FALSE)
  
  # Motif plot logo 3
  output$motif3 <- renderImage({
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    motif_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "clusters/",
        input$inputNumClust,
        "/MEME/",
        sep = ""
      )
    motif_file <- Sys.glob(file.path(motif_folder, "/logo3.png"))
    req(file.exists(motif_file))
    
    list(
      src = motif_file,
      style = "object-fit: contain; max-height: 95%; max-width: 95%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Enriched motif does not exist for this cluster."
    )
  }, deleteFile = FALSE)
  
  #################### Secondary Analysis ##################
  ###################### Factor Binding ####################
  
  # Return transcription factor  table
  output$tfTable <- output$tfTable2 <- renderDataTable({
    req(input$runDE)
    req(analysis_folder_name())
    req(local_results_folder())
    tf_table_file <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/factor_binding/observed_putative_tfs_n_encode.csv",
        sep = ""
      )
    
    if (!file.exists(tf_table_file)) {
      getTopTfs()
    }
    req(file.exists(tf_table_file))
    
    d <-
      read.table(
        tf_table_file,
        header = T,
        row.names = count.fields(tf_table_file, sep = ",")[1],
        sep = ","
      )
    data <- d[order(row.names(d)), ] # order rows by cluster number
    return(as.data.frame(data))
  }, options = (list(
    pageLength = 5,
    scrollX = TRUE,
    autoWidth = FALSE
  )))
  
  # Get the top transcription factors
  getTopTfs <- function() {
    avg_prof_script <-
      paste(app_dir, "/scripts/get_top_tfs.r", sep = "")
    res_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        sep = ""
      )
    command_avg_prof <-
      paste("Rscript",
            avg_prof_script,
            res_folder,
            input$organism,
            3,
            4,
            40,
            app_dir,
            sep = " ")
    system(command_avg_prof, intern = TRUE)
  }
  
  # Rcistarget interactive results download
  output$interactiveRcisResults <- downloadHandler(
    filename = function() {
      paste("all_method_motifs_interactive_rcistarget_results",
            Sys.Date(),
            ".html",
            sep = "")
    },
    content = function(file) {
      interRcis_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/analysis/",
          analysis_folder_name(),
          "_results/",
          "clusters/",
          input$rcisClustNum,
          "/top_tfs/",
          sep = ""
        )
      interRcis_html <-
        Sys.glob(file.path(
          interRcis_folder,
          "/interactive_rcistarget_results.html"
        ))
      writeLines(readLines(interRcis_html), file)
    }
  )
  
  # Return consensus RcisTarget results
  output$rcistargetConsensus <- DT::renderDataTable({
    req(diffExpDone())
    rcis_consens_file <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/clusters/",
        input$rcisClustNum,
        "/top_tfs/top_TFs_per_method_n_encode_nums.csv",
        sep = ""
      )
    rcis_consensus <<- read.csv(rcis_consens_file, header = TRUE)
    rcis_consensus
  }, options = (list(
    pageLength = 5, scrollX = TRUE
  )))
  
  # Create average profiles from .bigWig inputs
  runBigWig <- function(bw, name) {
    avg_prof_script <-
      paste(app_dir, "/scripts/create_avg_prof.sh", sep = "")
    res_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        sep = ""
      )
    command_avg_prof <-
      paste(
        avg_prof_script,
        bw,
        name,
        res_folder,
        paste(res_folder, "/factor_binding/", sep = ""),
        paste(
          res_folder,
          "/clusters/cluster_color_to_number_map.txt",
          sep = ""
        ),
        sep = " "
      )
    system(command_avg_prof, intern = TRUE)
  }
  
  # 1st .bigWig file
  output$bigWigResults1a <- renderImage({
    if(sim_demo_data() == TRUE){
      req(input$renderVen)
      avg_prof_png_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/analysis/",
          analysis_folder_name(),
          "_results/factor_binding/",
          sep = ""
        )
      pngName <- paste("avg_profile", avg_prof_name(), "png", sep = ".")
      avg_prof_png_file <- paste(avg_prof_png_folder, pngName, sep = "/")
      
      print("avgprof")
      print(avg_prof_png_file)
      
      list(
        src = avg_prof_png_file,
        width = "100%",
        style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
        contentType = 'image/png',
        alt = "Image for the first bigwig does not exist"
      )

    } else{
      req(input$bigWigGo1)
      req(!is.null(input$tf1))
      req(!is.null(input$bigWig1))
      
      avg_prof_png_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/analysis/",
          analysis_folder_name(),
          "_results/factor_binding/",
          sep = ""
        )
      pngName <- paste("avg_profile", input$tf1, "png", sep = ".")
      avg_prof_png_file <- paste(avg_prof_png_folder, pngName, sep = "/")
      
      print("avgprof")
      print(avg_prof_png_file)
      
      # Create average profiles if file does not exist
      if (!file.exists(avg_prof_png_file)) {
        print("creating profile")
        
        runBigWig(input$bigWig1$datapath, input$tf1)
      }
      req(file.exists(avg_prof_png_file))
      list(
        src = avg_prof_png_file,
        width = "100%",
        style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
        contentType = 'image/png',
        alt = "Image for the first bigwig does not exist"
      )
    }
  }, deleteFile = FALSE)
  
  output$bigWigResults1b <- renderImage({
    if(sim_demo_data() == TRUE){
      req(input$renderVen)
      
      heatmap_png_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/analysis/",
          analysis_folder_name(),
          "_results/factor_binding/",
          sep = ""
        )
      pngName <-
        paste("heatmap_genes.clusters", avg_prof_name(), "png", sep = ".")
      heatmap_png_file <- paste(heatmap_png_folder, pngName, sep = "/")
      
      print("heatmap")
      print(heatmap_png_file)
      req(file.exists(heatmap_png_file))
      list(
        src = heatmap_png_file,
        width = "100%",
        style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
        contentType = 'image/png',
        alt = "Image for the first bigwig does not exist"
      )
    } else{
      req(input$bigWigGo1)
      req(!is.null(input$tf1))
      req(!is.null(input$bigWig1))
      
      heatmap_png_folder <-
        paste(
          local_results_folder(),
          "/timeor/results/analysis/",
          analysis_folder_name(),
          "_results/factor_binding/",
          sep = ""
        )
      pngName <-
        paste("heatmap_genes.clusters", input$tf1, "png", sep = ".")
      heatmap_png_file <- paste(heatmap_png_folder, pngName, sep = "/")
      
      print("heatmap")
      print(heatmap_png_file)
      req(file.exists(heatmap_png_file))
      list(
        src = heatmap_png_file,
        width = "100%",
        style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
        contentType = 'image/png',
        alt = "Image for the first bigwig does not exist"
      )
    }
  }, deleteFile = FALSE)
  
  # 2nd .bigWig
  output$bigWigResults2a <- renderImage({
    req(input$bigWigGo2)
    req(!is.null(input$tf2))
    req(!is.null(input$bigWig2))
    
    avg_prof_png_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/factor_binding/",
        sep = ""
      )
    pngName <- paste("avg_profile", input$tf2, "png", sep = ".")
    avg_prof_png_file <- paste(avg_prof_png_folder, pngName, sep = "/")
    
    print("avgprof")
    print(avg_prof_png_file)
    
    # Create average profiles if file does not exist
    if (!file.exists(avg_prof_png_file)) {
      print("creating profile")
      
      runBigWig(input$bigWig2$datapath, input$tf2)
    }
    req(file.exists(avg_prof_png_file))
    list(
      src = avg_prof_png_file,
      width = "100%",
      style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Image for the first bigwig does not exist"
    )
  }, deleteFile = FALSE)
  
  output$bigWigResults2b <- renderImage({
    req(input$bigWigGo2)
    req(!is.null(input$tf2))
    req(!is.null(input$bigWig2))
    
    heatmap_png_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/factor_binding/",
        sep = ""
      )
    pngName <-
      paste("heatmap_genes.clusters", input$tf2, "png", sep = ".")
    heatmap_png_file <- paste(heatmap_png_folder, pngName, sep = "/")
    
    print("heatmap")
    print(heatmap_png_file)
    req(file.exists(heatmap_png_file))
    list(
      src = heatmap_png_file,
      width = "100%",
      style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Image for the first bigwig does not exist"
    )
  }, deleteFile = FALSE)
  
  # 3rd .bigWig
  output$bigWigResults3a <- renderImage({
    req(input$bigWigGo3)
    req(!is.null(input$tf3))
    req(!is.null(input$bigWig3))
    
    avg_prof_png_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/factor_binding/",
        sep = ""
      )
    pngName <- paste("avg_profile", input$tf3, "png", sep = ".")
    avg_prof_png_file <- paste(avg_prof_png_folder, pngName, sep = "/")
    
    print("avgprof")
    print(avg_prof_png_file)
    
    # Create average profiles if file does not exist
    if (!file.exists(avg_prof_png_file)) {
      print("creating profile")
      
      runBigWig(input$bigWig3$datapath, input$tf3)
    }
    req(file.exists(avg_prof_png_file))
    list(
      src = avg_prof_png_file,
      width = "100%",
      style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Image for the first bigwig does not exist"
    )
  }, deleteFile = FALSE)
  
  output$bigWigResults3b <- renderImage({
    req(input$bigWigGo3)
    req(!is.null(input$tf3))
    req(!is.null(input$bigWig3))
    
    heatmap_png_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/factor_binding/",
        sep = ""
      )
    pngName <-
      paste("heatmap_genes.clusters", input$tf3, "png", sep = ".")
    heatmap_png_file <- paste(heatmap_png_folder, pngName, sep = "/")
    
    print("heatmap")
    print(heatmap_png_file)
    req(file.exists(heatmap_png_file))
    list(
      src = heatmap_png_file,
      width = "100%",
      style = "object-fit: contain; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
      contentType = 'image/png',
      alt = "Image for the first bigwig does not exist"
    )
  }, deleteFile = FALSE)
  
  #################### Secondary Analysis ##################
  #################### Temporal Relations ##################
  
  # Read in gene sample input file from CSV
  tf_data <- reactive ({
    req(input$tfTable)
    file1 <- input$tfTable
    if (is.null(file1)) {
      return()
    }
    read.table(file = file1$datapath, dec = ",")
  })
  
  output$tfTempTable <- DT::renderDataTable({
    req(input$renderVen)
    req(analysis_folder_name())
    req(local_results_folder())
    req(diffExpDone())
    
    # TF_temp_rel.csv file
    tf_temp_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/temporal_relations/",
        sep = ""
      )
    tf_temp_file <-
      Sys.glob(file.path(tf_temp_folder, "/TF_temp_rel.csv"))
    
    # Run get_tf_relations.r if needed
    if (!file.exists(tf_temp_file)) {
      getTFrelations()
    }
    
    # Visualizing temporal relations table
    req(file.exists(tf_temp_file))
    df <- as.data.frame(read.csv(tf_temp_file, header = T))
    options(DT.options = list(pageLength = 6, scrollX = TRUE))
    datatable(df) %>% formatStyle('regulation_type', backgroundColor = styleEqual(
      c(
        'obs_obs_known_int',
        'obs_obs_pred_int',
        'pred_obs_known_int',
        'pred_obs_pred_int'
      ),
      c('#e65d8e', '#FEC02F', '#288AE2', '#508279')
    ))
  })
  
  getTFrelations <- function() {
    tf_rel_script <-
      paste("Rscript ", app_dir, "/scripts/get_tf_relations.r", sep = "")
    res_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        sep = ""
      )
    
    # Network (STRINGdb) script calls
    if (input$organism == "dme") {
      ncbi_id = 7227
    } else if (input$organism == "hse") {
      ncbi_id = 9606
    } else{
      ncbi_id = 10090
    }
    
    command_tf_rel <-
      paste(tf_rel_script, res_folder, ncbi_id, app_dir, sep = " ")
    system(command_tf_rel, intern = TRUE)
  }
  
  output$TF_pval_network_enrichment <- renderText({
    # Set requirements
    req(input$runEnrichment)
    req(analysis_folder_name())
    req(input$inputNumClust != "NA")
    
    # Get p-value
    TF_network_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "temporal_relations/",
        sep = ""
      )
    TF_network_pval_file <-
      file.path(TF_network_folder, "/stringdb_info_table.tsv") # get StringDB network information path
    req(file.exists(TF_network_pval_file))
    return(read.table(TF_network_pval_file, header = T)$'pval_PPI_enrich'[1])
  })
  
  # Network plot # aspect ratio and disable scrolling and remove white space
  output$TFnetwork <- renderImage({
    # Set requirements
    req(input$renderVen)
    req(analysis_folder_name())
    req(local_results_folder())
    req(diffExpDone())
    
    # Get png
    tf_network_folder <-
      paste(
        local_results_folder(),
        "/timeor/results/analysis/",
        analysis_folder_name(),
        "_results/",
        "temporal_relations/www/",
        sep = ""
      )
    tf_network_png <-
      Sys.glob(file.path(tf_network_folder, "/*stringdb_network*")) # get StringDB network path
    req(file.exists(tf_network_png))
    
    # Format and get png
    list(
      src = tf_network_png,
      style = "object-fit: contain; max-height: 100%; max-width: 100%; display: block; margin-left: auto; margin-right: auto;",
      #height = "auto", width = "auto",
      contentType = 'image/png',
      alt = "Network does not exist."
    )
  }, deleteFile = FALSE)
  
  # STRINGdb website for user to explore temporally related transcription factors and other genes (if desired)
  output$stringDB_web <- renderUI({
    string_db <-
      tags$iframe(src = "https://string-db.org/cgi/input",
                  height = "1300",
                  width = "100%")
    string_db
  })
  
  
  ######################## Tutorial ########################
  output$web <- renderUI({
    #Shiny takes .md so convert .Rmd to .md with: library("rmarkdown"), render("timeor_app_tutorial.Rmd", output_format = "md_document")
    withMathJax(includeMarkdown(
      paste(app_dir, "/tutorial/timeor_app_tutorial.md", sep = "")
    ))
  })
  
  output$commandLine <- renderUI({
    withMathJax(includeMarkdown(
      paste(
        app_dir,
        "/tutorial/timeor_command_line_tutorial.md",
        sep = ""
      )
    ))
  })
  
  # Download zipped results folder
  output$downloadResultsFolder <- downloadHandler(
    
    # Download simulated or real data used
    if(sim_demo_data() || real_demo_data()){
      filename = function() {
        paste("timeor", "tar.gz", sep = ".")
      },
      content = function(folderName) {
        print(paste0(
          "Local results folder: ",
          local_results_folder(),
          "/timeor/"
        ))
        file.copy(paste(local_results_folder(), "/timeor.tar.gz", sep = ""), folderName)},
        #tar(folderName, paste(local_results_folder(), "/timeor/", sep = ""))}
        contentType = "application/zip"
    }else{

      # Download new data 
      filename = function() {
        paste("timeor", "tar", "gz", sep = ".")
      },
      content = function(folderName) {
        command <- paste("tar -czvf ", local_results_folder(),"/timeor.tar.gz ",local_results_folder(),"/timeor", sep = "")
        cat(command)
        system(paste("tar --usage", "> /tmp/file_output.txt",  sep=" "))# , intern = TRUE)
        system(command)# , intern=TRUE)
        file.copy(paste(local_results_folder(), "/timeor.tar.gz", sep = ""), folderName)},
        contentType = "application/octet-stream"
    }
  )
}