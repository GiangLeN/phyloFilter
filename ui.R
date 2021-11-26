library("shiny")
library("tidyverse")
library("phyloseq")
library("BiocManager")
library("DT")
library("methods")
library("ape")
library("Biostrings")

options(repos = BiocManager::repositories())

taxLevel <- c("Kingdom","Phylum", "Class", "Order")

fluidPage(
  
  titlePanel("Taxa filtering"),
  
  sidebarPanel(
    
    textInput("name", "Project's name", value = 'phyloseq'),
    fileInput("file1", "Choose Phyloseq File", accept = ".rds"),
    radioButtons("taxaLevel", "Select taxonomic level", selected = "Phylum", taxLevel),
    numericInput("prevaTaxa", "Prevalance, % of samples with ASVs", value = 5, min = 0, max = 100),
    numericInput("abunTaxa", "Abundance, % of ASV overall", value = 0.01, min = 0, max = 100),
    sliderInput("pointSize", "Sample's size", value = 1.5, min = 1, max = 15),
    sliderInput("pointTrans", "Transparent", value = 0.5, min = 0, max = 1),
    h4("How to Exclude Taxa:"),
    h5("Advance option: Do not use unless you know what you are doing."),
    h5("Type in the taxa names seperate by | without any space"),
    h5("eg: Fusobacteriota|Cyanobacteria "),
    textInput("excludeTaxa", "Exclude taxa"),
    downloadButton("downloadTable", "ASV table"),
    downloadButton("downloadData", "Final phyloseq"),
    downloadButton("downloadTree", "Tree"),
    h4("Rhea inputs"),
    downloadButton("downloadRheaASV", "Rhea ASVs"),
    downloadButton("downloadRheaMapping", "Rhea Mapping"),
    downloadButton("downloadRheaSeq", "Rhea Sequences"),
    h4("MicrobiomeAnalyst inputs"),
    downloadButton("downloadMAnalystASV", "ASV table"),
    downloadButton("downloadMAnalystTaxa", "Taxa table"),
    downloadButton("downloadMAnalystMapping", "Metafile"),
    h5("Note: The metafile should contain no columns with NA or blank spaces.")

    
  ),
  
  mainPanel(
    verbatimTextOutput("inPs"),
    textOutput("summary"),
    textOutput("phyloseq"),
    textOutput("filterSummary"),
    plotOutput("taxaFilter", height=1000),
    DT::dataTableOutput("taxaTable"),
    textOutput("phyloSave"),
    verbatimTextOutput("saveFile"),
#    verbatimTextOutput("Rhea"),
    textOutput("unique")
    
  )
  
)

##rsconnect::deployApp('Work/learnedToday/App/', appName = "phyloFilter")