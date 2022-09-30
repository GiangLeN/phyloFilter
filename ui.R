library("shiny")
library("tidyverse")
library("phyloseq")
library("BiocManager")
library("DT")
library("methods")
library("ape")
library("Biostrings")

options(repos = BiocManager::repositories())

taxLevel <- c("Kingdom","Phylum", "Class")

fluidPage(
  
  titlePanel("Phyloseq taxa filter"),
  
  sidebarPanel(
    
    fileInput("file1", "Choose Phyloseq File", accept = ".rds"),
    textInput("name", "Project's name", value = 'phyloseq'),
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
    downloadButton("downloadRheaASV", "Rhea_ASV"),
    downloadButton("downloadRheaMapping", "Rhea_Mapping"),
    downloadButton("downloadRheaSeq", "Rhea_Sequences"),
    h4("MicrobiomeAnalyst inputs"),
    downloadButton("downloadMAnalystASV", "MA_ASVs"),
    downloadButton("downloadMAnalystTaxa", "Taxa_table"),
    downloadButton("downloadMAnalystMapping", "Metafile"),
    h5("Note: The metafile should contain no columns with NA or blank spaces."),
    h5("Version: 0.1"),
    tags$footer(
      "For detail instruction please visit ",
      tags$a(
        "https://github.com/GiangLeN/phyloFilter",
        target = "_blank",
        href = "https://github.com/GiangLeN/phyloFilter"
      ),
      style = "position: absolute; width: 100%; color: black; text-align: left;"
    ),
    


  ),
  
  mainPanel(
    strong(textOutput("info")),
    verbatimTextOutput("inPs"),
    textOutput("samdat"),
    textOutput("refseq"),
    textOutput("phytree"),
    textOutput("phyloseq"),
    br(),
    strong(textOutput("method")),
    textOutput("filterSummary"),
    br(),
    plotOutput("taxaFilter", height=1000),
    DT::dataTableOutput("taxaTable"),
    textOutput("phyloSave"),
    verbatimTextOutput("saveFile"),
    verbatimTextOutput("Rhea"),
    textOutput("manual"),
    
    
  )
  
)

##rsconnect::deployApp('Work/learnedToday/App/', appName = "phyloFilter", account = "giangle")