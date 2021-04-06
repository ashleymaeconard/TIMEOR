# ui.R
# Ashley Mae Conard
# Last Mod. 6 June 2020

library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(heatmaply)
library(d3heatmap)
library(shinyjs)
library(shinyalert)
library(shinyWidgets)
library(shinycssloaders)
library(knitr)
options(shiny.sanitize.errors = FALSE)

# Displaying TIMEOR logo
title <-
  tags$a(tags$img(
    src = "timeor_logo.png",
    height = "auto",
    width = "100%"
  ))

# Runs this UI code with each new session
function(request) {
  
  # TIMEOR UI page
  fluidPage(
    tags$head(HTML("\n<!-- Global site tag (gtag.js) - Google Analytics -->
                  <script async src=\"https://www.googletagmanager.com/gtag/js?id=G-D4Z8R0914Q\"></script>
                  <script>
                    window.dataLayer = window.dataLayer || [];
                    function gtag(){dataLayer.push(arguments);}
                    gtag('js', new Date());

                    gtag('config', 'G-D4Z8R0914Q');
                  </script>")),
    tags$head(tags$script(HTML("Shiny.addCustomMessageHandler ('print',function (message) {
              $('#'+message.selector).html(message.html);
              console.log(message);});"))),
    useShinyjs(), # perform common useful JavaScript operations in Shiny apps
    useShinyalert(), # alert the user to important information
    
    # Setting dashboard style
    dashboardPage(
      title="TIMEOR", # browser title
      skin = "black",
      
      # Setting leftmost sidebar (i.e. dashboard)
      dashboardHeader(title = title),
      dashboardSidebar(
        sidebarMenu(
          HTML("",sep="<br/>"), # new line
          menuItem(
            "Pre-process",
            tabName = "pre-process",
            icon = icon("chart-bar")
          ),
          menuItem(
            "Primary Analysis",
            tabName = "primAnalysis",
            icon = icon("sitemap")
          ),
          menuItem(
            "Secondary Analysis",
            tabName = "secondAnalysis",
            icon = icon("project-diagram")
          ),
          menuItem("Tutorials", tabName = "tutorials_ft_bk", icon = icon("book")),
          HTML("",sep="<br/>"), # new line
          menuItem(
            "Example Data",
            tabName = "exampleData",
            icon = icon("table"),
            menuItem(
              "Load raw data",
              actionButton("runReal", "SraRunTable & .fastq files", style =
                             "margin: 6px 5px 6px 0px; width:90%")
            ),
            menuItem(
              "Load count matrix",
              actionButton("runSim", "Metadata & read count file", style =
                             "margin: 6px 5px 6px 0px; width:90%")
            )
          ),
          menuItem(
            "Download Results Folder",
            tabName = "resultsFolder",
            icon = icon("download"),
            downloadButton("downloadResultsFolder", "Download results", style =
                             "width:90%;color:#444")
          )
        )
      ),
      
      # Style settings for body of app (right of sidebar (i.e. dashboard))
      dashboardBody(
        tags$head(
          tags$style(type = "text/css",
                     "#image img {max-width: 100%; width: 100%; height: auto}")
        ),
        
        # Body of app (right of sidebar (i.e. dashboard))
        tabItems(
          
#################### Pre-Process Stage ##################
##################### Process Raw Data ##################
          
          tabItem(tabName = "pre-process",
                  tabsetPanel(
                    tabPanel(
                      "Process Raw Data",
                      fluidRow(
                        box(
                          
                          # Search and retrieve text within upper leftmost app body within Process Raw Data Tab
                          titlePanel(h3("Search and Retrieve")),
                          width = 6,
                          h4(
                            "Processing your own data? Jump to step 4 to load your metadata file, answer 6 questions on the right panel, then proceed to the 'Process Count Matrix' tab."
                          ),
                          tags$hr(), # horizontal line 
                          h4(
                            "1. Go to ",
                            tags$a(href = "https://www.ncbi.nlm.nih.gov/sra", "GEO", target = "_blank"),
                            "to search for a time series RNA-seq dataset."
                          ),
                          h4("2. Click on \"SRA Run Selector\" on bottom right."),
                          h4(
                            "3. Download data (all or selected replicates) by clicking on \"Metadata\" under \"Select\"."
                          ),
                          h4(
                            "4. Set input file type to either the resulting \"SraRunTable.txt\" or a pre-made metadata file."
                          ),
                          switchInput(
                            inputId = "sra_or_meta",
                            label = "Input file",
                            onLabel = "Metadata",
                            offLabel = "SraRunTable"
                          ),
                          h4("5. Upload input file. Note: SraRunTables are automatically converted to metadata files."),
                          fileInput(
                            "metadataFile",
                            label = NULL,
                            accept = c(
                              "text/csv",
                              "text/comma-separated-values,text/plain",
                              ".txt"
                            )
                          ),
                          DT::dataTableOutput("metadataTable"),
                          style="background: #E8E8E8"),
                        box(style="background: #E8E8E8",
                          
                          # Determine adaptive default parameters text within upper rightmost app body within Process Raw Data Tab
                          titlePanel(h3("Determine Adaptive Default Methods")),
                          h4("1. What type of organism?"),
                          selectInput(
                            "organism",
                            label = NULL,
                            choices = c(
                              "Select an organism" = "NA",
                              "Drosophila melanogaster" = "dme",
                              "Homo sapiens" = "hse",
                              "Mus musculus" = "mus"
                            )
                          ),
                          h4("2. What type of sequencing?"),
                          selectInput(
                            "sequencing",
                            label = NULL,
                            selected = "Select a sequencing type",
                            choices = c(
                              "Select a sequencing type" = "NA",
                              "Paired-end" = "pe",
                              "Single-end" = "se"
                            )
                          ),
                          h4("3. What type of experiment?"),
                          selectInput(
                            "experiment",
                            label = NULL,
                            selected = "Select an experiment type",
                            choices = c(
                              "Select an experiment type" = "NA",
                              "Case versus control" = "cc",
                              "Just case or control" = "jc"
                            )
                          ),
                          h4("4. What type of time series?"),
                          selectInput(
                            "timeSeries",
                            label = NULL,
                            selected = "Close time point and long time series",
                            choices = c(
                              "Close time point and long time series" = "ct",
                              "Close time point and short time series" = "ctstc",
                              "Distant time point" = "distant"
                            )
                          ),
                          h4(
                            "5. Compare multiple methods (alignment and differential expression)?"
                          ),
                          selectInput(
                            "multipleMethods",
                            label = NULL,
                            selected = "Yes",
                            choices = c("Yes", "No")
                          ),
                          h4(
                            "6. What is the maximum number of time steps over which one gene can influence the transcription of another gene?"
                          ),
                          selectInput(
                            "trans_time",
                            label = NULL,
                            selected = "1",
                            choices = c("1", "2", "3", "4", "5", "6", "7")
                          ),
                          actionButton("run_adaptive_defaults", label = "Run"),
                          htmlOutput("dataProcessType"), tags$head(tags$style("#dataProcessType{color: #377BB5;
                                 font-size: 18px;text-align: center; font-weight: bold}")
                          )
                        )
                      ),
                      tags$hr(style = "border-color: black;"), # horizontal line 
                      
                      # Process, quality control, and alignment quality text within lower leftmost text in app body within Process Raw Data Tab
                      fluidRow(
                        box(
                          titlePanel(h3("Process")),
                          height = 500,
                          width = 3,
                          tags$table(
                            tags$tr(tags$td(h4("1. Retrieve Data"), style = "padding: 15px"),
                                    tags$td(h2(
                                      span(textOutput('textCheckMark_raw'), style = "color:green")
                                    ))),
                            tags$tr(tags$td(h4(
                              "2. Quality Control"
                            ), style = "padding: 15px"),
                            tags$td(h2(
                              span(textOutput('textCheckMark_qc'), style = "color:green")
                            ))),
                            tags$tr(tags$td(h4("3. Align Data"), style =
                                              "padding: 15px"),
                                    tags$td(h2(
                                      span(textOutput('textCheckMark_align'), style = "color:green")
                                    ))),
                            tags$tr(tags$td(h4(
                              "4. Generate Count Matrix"
                            ), style = "padding: 15px"),
                            tags$td(h2(
                              span(textOutput('textCheckMark_count'), style = "color:green")
                            ))),
                            style = "width:100%"
                          ),
                          uiOutput("withprogress")
                        ),
                        box(
                          height = 500,
                          width = 9,
                          titlePanel(h3("Quality Control")),
                          downloadButton("qcDownload", "Download interactive results summary."),
                          fluidRow(tags$br(), tags$br(), 
                                   box(width=12,DT::dataTableOutput('qcTable'))
                          )
                        )),
                      fluidRow(box(style="background: #E8E8E8",
                        width=12, height=700, 
                        titlePanel(h3("Alignment Quality")),
			h4("Which alignment method?"),
                          fluidRow(column(4,column(8,selectInput( width='70%',
                            "aligner",
                            label = NULL,
                            choices = c(
                              "HISAT2" = "hisat2",
                              "Bowtie2" = "bowtie2"
                            ))),
                          column(4,actionButton("runGenCountMat", label = "Generate count matrix")))),
                        box(imageOutput("alignment_graph_HISAT2")), 
                        box(imageOutput("alignment_graph_Bowtie2"))))
                    ),
                    
#################### Pre-Process Stage ##################
#################### Process Count Matrix ##################
                    
                    tabPanel(
                      "Process Count Matrix",
                      fluidRow(box(width=4, # sidebar panel for read count matrix inputs
                          
                          # Input a read count matrix
                          fileInput(
                            "count_mat",
                            titlePanel(h3("Upload Count Matrix")),
                            accept = c(
                              "text/csv",
                              "text/comma-separated-values,text/plain",
                              ".csv"
                            )
                          ),
                          tags$hr(), # horizontal line 
                          checkboxInput("header", "Header", TRUE),# input: checkbox if file has header
                          radioButtons(# input select separator
                            "sep",
                            "Separator",
                            choices = c(
                              Comma = ",",
                              Semicolon = ";",
                              Tab = "\t"
                            ),
                            selected = ","
                          ),
                          radioButtons( # input: select quotes
                            "quote",
                            "Quote",
                            choices = c(
                              None = "",
                              "Double Quote" = '"',
                              "Single Quote" = "'"
                            ),
                            selected = '"'
                          ),
                          tags$hr(), # horizontal line 
                          radioButtons( # input select number of rows to display
                            "disp",
                            "Display",
                            choices = c(Head = "head",
                                        All = "all"),
                            selected = "head"
                          ), style="background: #E8E8E8"
                        ),
                        
                        # Display read count matrix
                        box(titlePanel(h3("Count Matrix")), width=8,
                            withSpinner(DT::dataTableOutput("countMatrix")))
                      ),
                      
                      # Show raw read count matrix principal component analysis (PCA)
                      fluidRow(
                        box(titlePanel(h3("Principal Component Analysis")), height=600,
                            withSpinner(plotlyOutput("pcaScatBefore"))),
                        
                        # Show raw read count matrix replicate (i.e. sample) correlations
                        box(titlePanel(h3("Replicate Correlations")), height=600,
                          selectInput(
                            "correMethodBefore",
                            label = "Correlation Methods",
                            choices = c("Spearman", "Pearson")
                          ),
                          withSpinner(plotlyOutput('correlationBefore'))
                        )
                      ),
                      
                      # Show accompanying loadings plot for read count matrix PCA
                      fluidRow(box(titlePanel(h3("Loadings")), height=600,
                                   withSpinner(
                                     plotOutput("pcaBarBefore")
                                   )))
                    ),
                    
#################### Pre-Process Stage ##################
############### Normalize and Correct Data ##############
                    
                    tabPanel(
                      "Normalize and Correct Data",
                      fluidRow(
                        
                        # Normalization and correction method user options
                        box(titlePanel(h3("Choose a Normalization and Correction Method")),
                          selectInput(
                            "normMethods",
                            label = "Normalization Methods",
                            choices = c("Upper Quartile", "Trimmed Mean of M-Values")
                          ),
                          selectInput(
                            "corrMethods",
                            label = "Correction Methods",
                            choices = c("Harman")
                          ),
                          actionButton("runCorrection", label = "Run"),
                          style="background: #E8E8E8"),
                        
                        # Show normalized and corrected read count matrix PCA
                        box(titlePanel(h3("Principal Component Analysis")), 
                            withSpinner(plotlyOutput("pcaScatAfter")))
                      ),
                      fluidRow(
                        
                        # Show normalized and corrected accompanying loadings plot for read count matrix PCA
                        box(titlePanel(h3("Loadings")), height=600,
                            withSpinner(plotOutput("pcaBarAfter"))),
                        
                        # Show normalized and corrected replicate (i.e. sample) correlations
                        box(titlePanel(h3("Replicate Correlations")), height=600,
                          selectInput(
                            "correMethodAfter",
                            label = "Correlation Methods",
                            choices = c("Spearman", "Pearson")
                          ),
                          withSpinner(plotlyOutput('correlationAfter'))
                        )
                      )
                    )
                  )),
          
#################### Primary Analysis ##################
          
          tabItem(tabName = "primAnalysis",
                  fluidRow(
                    
                    # Set differential expression options
                    box(titlePanel(h3("Run Differential Expression")),
                      textOutput("deParameters"),
                      height = "200",
                      textInput("analysisFolder", label = "What do you want to name your results folder?", value = ""),
                      textInput("pValue", label = "What is your adjusted p-value threshold?", value = ""),
                      tags$hr(), # horizontal line 
                      actionButton("runDE", label = "Run"),
                      style="background: #E8E8E8"),
                    
                    # Display Venn diagram overlap plot between methods
                    box( titlePanel(h3("Venn Diagram Between Methods")),
                      fileInput(
                        "prevStudy",
                        label = "Compare differential expression results (and with previous study).",
                        accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".txt"
                        )
                      ),
                      actionButton("renderVen", label = "Render Venn Diagram"),
                      imageOutput("vennMethods", width =
                                    "auto")
                    )
                  ),
                  fluidRow(
                    box(titlePanel(h3("Display Desired Differential Expression Method Results")),
                      selectInput(
                        "whichMethodInput",
                        label = "Which method?",
                        choices = c("ImpulseDE2", "NextMaSigPro", "DESeq2")
                      ),
                      withSpinner(DT::dataTableOutput("deResults"))
                    ),
                    box(titlePanel(h3("Cluster Gene Expression Trajectories")),
                      selectInput(
                        "numClusts",
                        "Select the number of clusters.",
                        choices = c(
                          "automatic", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"
                        )
                      ),
                      withSpinner(plotlyOutput("deClustering"))
                      )
                  )),
          
#################### Secondary Analysis ##################
######################## Enrichment ######################
          
          tabItem(tabName = "secondAnalysis",
                  tabsetPanel(
                    tabPanel(
                      "Enrichment",
                      fluidRow(
                        box(
                          titlePanel(h3("Gene Expression Trajectory Clusters")),
                          selectInput(
                            "inputNumClust",
                            label = paste("Clusters are labeled in ascending order from 1 for top-most cluster. Select a gene trajectory cluster to analyze:", sep = ""),
                            selected = "NA",
                            choices = c(
                              "NA" = "NA", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"
                            )
                          ),
                          withSpinner(plotlyOutput("deClustering1")),
                          style="background: #E8E8E8"),
                        box(
                          titlePanel(h3("Chosen Cluster Gene Set")),
                          span(uiOutput('textWithHTML')),
                          headerPanel(""),
                          switchInput(inputId = "runEnrichment", label = "Analyze"),
                          style="background: #E8E8E8")#,
                          #label=paste("Should the process finish and images not show, simply turn the toggle 'OFF' and then 'ON' to view results.")
                      ),
                      fluidRow(
                        box(
                          titlePanel(h3("Molecular Function")),
                          width = 4,
                          height = 500,
                          withSpinner(imageOutput("molecularFunc", width =
                                                    "auto"))
                        ),
                        box(
                          titlePanel(h3("Biological Process")),
                          width = 4,
                          height = 500,
                          withSpinner(imageOutput("bioProc"))
                        ),
                        box(
                          titlePanel(h3("Cellular Component")),
                          width = 4,
                          height = 500,
                          withSpinner(imageOutput("cellComp"))
                        )
                      ),
                      fluidRow(
                        box(
                          titlePanel(h3("Pathway")),
                          width = 6,
                          height = "500",
                          withSpinner(imageOutput("pathway"))
                        ),
                        box(
                          titlePanel(h3("Network")),
                          width = 5,
                          height = "500",
                          fluidRow(column(
                            6,
                            strong("Protein-protein interaction enrichment p-value:"),
                            uiOutput('pval_network_enrichment')
                          )),
                          #height = 500,
                          withSpinner(imageOutput("network"))
                        ),
                        box(
                          titlePanel(h3("Motif Analysis ")),
                          height = "auto",
                          width = 6,
                          downloadButton("memeHTML", "Download MEME interactive results."),
                          fluidRow(column(
                            6, tags$br(), strong("Top ", em("de novo"), " motif:")
                          ),
                          withSpinner(imageOutput("motif"))),
                          fluidRow(column(6, strong(
                            "2nd top ", em("de novo"), " motif:"
                          )),
                          withSpinner(imageOutput("motif2"))),
                          fluidRow(column(6, strong(
                            "3rd top ", em("de novo"), " motif:"
                          )),
                          withSpinner(imageOutput("motif3")))
                        )
                      )
                    ),
                    
#################### Secondary Analysis ##################
###################### Factor Binding ####################
                    
                    tabPanel("Factor Binding",
                             fluidRow(
                               ## Data table for transcription factor info table
                               box(titlePanel(h3("Gene Expression Trajectory Clusters")),
                                   withSpinner(plotlyOutput("deClustering2"))),
                               
                               box(titlePanel(h3("Observed and Top Predicted Transcription Factors",h4("At least 40% of methods must agree on their top predicted transcription factors, otherwise cell is left blank. Row names are gene expression trajectory clusters."))),
                                   withSpinner(DT::dataTableOutput("tfTable"))),
                               boxPlus(collapsed = TRUE, collapsible = TRUE, closable = FALSE, title="See each method's predicted transcription factors:", titlePanel(
                                   h3("Methods' Transcription Factor Prediction Rankings",
                                      h4("Blanks indicate an enriched motif is not assigned to a transcription factor region (to see motif click 'Download interactive cluster motif result'). Row names are top 1 - 4 transcription factors."))),
                                       selectInput(
                                        "rcisClustNum",
                                        label = "Which cluster?",
                                        choices = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"), 
                                      ),
                                       downloadButton("interactiveRcisResults", "Download interactive cluster motif results."),
                                      HTML("<br/>"),
                                      withSpinner(DT::dataTableOutput("rcistargetConsensus"))
                                )
                             ),
                             fluidRow(
                               box(
                                 
                                 # Create average profiles of user chosen transcription factors
                                 titlePanel(h3("Average Profiles Across Each Gene Expression Trajectory Cluster")),
                                 width = 12,
                                 box(
                                   h4(strong("Upload .bigWig Files")),
                                   height = "auto",
                                   width = 2,
                                   h5(strong("(private or public ChIP-seq or CUT&RUN)")),
                                   tags$br(),
                                   h4("If using public data:"),
                                   tags$br(),
                                   h4(
                                     "1. Go to ",
                                     tags$a(href = "https://www.encodeproject.org", target = "_blank", "ENCODE"),"."
                                   ),
                                   h4("2. Choose ENCODE data from provided ENCODE ID or TF name above."),
                                   h4("3. Download .bigWig file."),
                                   h4("4. Upload to TIMEOR to visualize."),
                                   tags$br(),
                                   h4(
                                     strong("Note"),
                                     ": each distribution is the same color as in clustermap."
                                   ),
                                   h5("ID: identifier, TF: transcription factor"),
                                   style="background: #E8E8E8"),
                                 
                                 box(
                                   textInput("tf1", label = "TF name:"),
                                   height = "auto",
                                   width = 3,
                                   fileInput("bigWig1", label = "Upload .bigWig"),
                                   actionButton("bigWigGo1", label = "Go"),
                                   tags$div(
                                     imageOutput("bigWigResults1a", width = "100%", height = "auto"),
                                     imageOutput("bigWigResults1b", width =
                                                   "100%", height = "auto")
                                   )
                                 ),
                                 
                                 box(
                                   textInput("tf2", label = "TF name:"),
                                   height = "auto",
                                   width = 3,
                                   fileInput("bigWig2", label = "Upload .bigWig"),
                                   actionButton("bigWigGo2", label = "Go"),
                                   tags$div(
                                     imageOutput("bigWigResults2a", width = "100%", height = "auto"),
                                     imageOutput("bigWigResults2b", width =
                                                   "100%", height = "auto")
                                   )
                                 ),
                                 
                                 box(
                                   textInput("tf3", label = "TF name:"),
                                   height = "auto",
                                   width = 3,
                                   fileInput("bigWig3", label = "Upload .bigWig"),
                                   actionButton("bigWigGo3", label = "Go"),
                                   tags$div(
                                     imageOutput("bigWigResults3a", width = "100%", height = "auto"),
                                     imageOutput("bigWigResults3b", width =
                                                   "100%", height = "auto")
                                   )
                                 )
                               )
                             )),
                    
#################### Secondary Analysis ##################
#################### Temporal Relations ################## 
                      
                    tabPanel("Temporal Relations",
                      fluidRow(
                        box(titlePanel(h3("Gene Expression Trajectory Clusters")),
                            withSpinner(plotlyOutput("deClustering3"))),
                        box(titlePanel(h3("Observed and Top Predicted Transcription Factors")),
                            DT::dataTableOutput("tfTable2"))
                      ),       
                      fluidRow(
                            column(
                          6,
                          box(width = "auto", height = "500", titlePanel(h3("Transcription Factor Network")),
                           column(9,
                                   height = "500",
                                   fluidRow(column(12,
                                     strong("Protein-protein interaction enrichment p-value:"),
                                     uiOutput('TF_pval_network_enrichment')
                                   )),
                                   withSpinner(
                                     imageOutput("TFnetwork", width = "100%", height = "400")
                                   )),
                            column(3,
                                   titlePanel(tags$a(
                                     tags$img(
                                       src = "TIMEOR_STRINGdbnetwork_web_vert.svg",
                                       width = "100%",
                                       height = "400"
                                     )
                                   ))
                            )
                          )
                          ),
                            column(
                          6,
                          box(width = "auto", height = "500", titlePanel(h3("Temporal Relations Between Observed and Predicted Transcription Factors"), "For edge_type 'a' is activation, 'r' is repression."),
                          fluidRow(
                           column(9,
                                    height = "500",
                                    DT::dataTableOutput("tfTempTable")
                                  ),
                            column(3,
                                   titlePanel(tags$a(
                                     tags$img(
                                       src = "TIMEOR_web_vert.svg",
                                       width = "100%",
                                       height = "350"
                                     )
                                   )))
                           )
                          )
                      )
                      ),    
                      fluidRow(height=50, 
                               box(width=12, h3("Network Customization: move and add desired genes to describe temporal relation"), "Click 'Search' then 'Multiple proteins' to begin adding gene names from above. User can also visit website directly (https://string-db.org/cgi/input).",style="background: #E8E8E8"),
                      fluidRow(height = 1500,
                              htmlOutput("stringDB_web")))
                    )
                  )),
          
######################## Tutorial ########################
          
          tabItem(tabName = "tutorials_ft_bk",
                  tabsetPanel(
                    tabPanel( "Overview and Contact",
                    fluidRow(box(width=12,
                    titlePanel(
                    tags$iframe(src="./overview_timeor.html", width='100%', height='1000px',frameborder=0,scrolling='auto') 
                    )),
                    #box(titlePanel(h3("Command Line TIMEOR", htmlOutput("commandLine"))))
                  )),tabPanel( "Web Server",
                    fluidRow(box(width=12,
                    titlePanel( 
                    tags$iframe(src="./timeor_app_tutorial.html", width='100%', height='1000px',frameborder=0,scrolling='auto')
                    ))
                  )),tabPanel( "Command Line",
                    fluidRow(box(width=12,
                    titlePanel( 
                    tags$iframe(src="./timeor_command_line_tutorial.html", width='100%', height='1000px',frameborder=0,scrolling='auto')
                    )))
                  ))
          )
        )
      )
    ),
    #bookmarkButton(label = "Save your place.")
  )
}