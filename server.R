library("shiny")
library("tidyverse")
library("phyloseq")
library("BiocManager")
library("DT")
library("methods")
library("ape")
library("Biostrings")

# fix package issue
options(repos = BiocManager::repositories())

writeFastaFile<-function(data, filename){
  
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste0(">", data[rowNum,"name"])))
    fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
  }
  
  ## Function to create, open and close connections (files, URLs etc)
  fileConnect<-file(filename)
  ## Write to file
  writeLines(fastaLines, fileConnect)
  ## Close once done
  close(fileConnect)
}


## Process rds file
function(input, output, session) {
  
  ## Load in phyloseq
  ps <- reactive({
    
    file <- input$file1
    ## Check file extension
    ext <- tools::file_ext(file$datapath)
    
    ## Require input file to run
    req(file)
    ## Validate phloseq file
    validate(need(ext == c("rds","RDS"), "Please upload a rds/RDS file"))
    
    ## Load the input file
    readRDS(file$datapath)
  })

  ## Phyloseq information
  output$inPs <- renderPrint({
    
    ## Check rds is phyloseq
    if (methods::is(ps(), "phyloseq")){
      
      ## Remove empty ASVs
      ps_clean <- reactive(prune_taxa(taxa_sums(ps()) > 0, ps()))
      
      ## Slots info from phyloseq 
      phyloSlots <- reactive(getslots.phyloseq(ps_clean()))
      
      if ("otu_table" %in% phyloSlots() & "tax_table" %in% phyloSlots()){
        ## Compute number of sample for each ASV, store as data.frame
        abunfltdf <- reactive({
          
          apply(X = otu_table(ps_clean()),
                MARGIN = ifelse(taxa_are_rows(ps_clean()), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
        })
        
        ## Add taxonomy and total read counts to this data.frame
        abunfltdata <- reactive({
          
          data.frame(Prevalence = abunfltdf(),
                     TotalAbundance = taxa_sums(ps_clean()),
                     tax_table(ps_clean())[, input$taxaLevel])
        })
        
        ## Calculate total reads and ignore na entries
        totalReads <- reactive({ sum(abunfltdata()$TotalAbundance, na.rm = TRUE) })
        
        ## Plot figures based on UI input
        output$taxaFilter <- renderPlot({
          
          ggplot(abunfltdata(), aes(TotalAbundance, Prevalence / nsamples(ps_clean()),color=get(input$taxaLevel))) +
            geom_hline(yintercept = input$prevaTaxa/100, alpha = 1, linetype = 2) +
            geom_vline(xintercept = input$abunTaxa * totalReads() / 100, alpha = 1, linetype = 3) +
            geom_point(size = input$pointSize, alpha = input$pointTrans) +
            scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
            facet_wrap(~get(input$taxaLevel), ncol = 2) + theme(legend.position="none")
          ## input$taxaLevel is a constant characters string, so need to have 'get'.
        })
        
        ## Create new table with the filtering settings
        prevaContamFreeTable <- reactive(
          
          abunfltdata() %>% rownames_to_column("ASVs") %>% 
            mutate(Prevalence_percentage = Prevalence/nsamples(ps_clean()) * 100, Abundance_percentage = TotalAbundance/totalReads()*100) %>%
            select(ASVs,input$taxaLevel,Prevalence,Prevalence_percentage,TotalAbundance,Abundance_percentage) %>%
            filter(Prevalence_percentage >= input$prevaTaxa & Abundance_percentage >= input$abunTaxa)
        )
        
        output$summary <- renderText(
          "Phyloseq information:"
        )
        
        ## ASV general information
        output$filterSummary <- renderText(
          
          paste0("An online shiny app (https://giangle.shinyapps.io/phyloFilter/) was used for ASV filtering. Any with number of reads below ",
                 ceiling(input$abunTaxa * totalReads() / 100), " (", input$abunTaxa, " %) and appeared in less than " ,
                 ceiling(input$prevaTaxa * nsamples(ps_clean()) / 100), " samples (", input$prevaTaxa, " %) were removed. A total of ",
                 nrow(prevaContamFreeTable ()), " ASVs were maintained for downstream analysis.",
                 "There were ", length(get_taxa_unique(finalPs(), taxonomic.rank = "Genus")), " unique genus and ",
                 length(get_taxa_unique(finalPs(), taxonomic.rank = "Species")), " unique species in the final phyloseq. About ",
                 round((ntaxa(ps_clean()) - nrow(prevaContamFreeTable ())) / ntaxa(ps_clean()) *100, digits = 2), " % of ASV was removed.")
        )
        
        output$method <- renderText({"Method:"})
        
        ## Render table
        output$taxaTable <- DT::renderDataTable({
          
          DT::datatable(prevaContamFreeTable(), options = list(orderClasses = TRUE)) %>%
            formatRound(columns='Prevalence_percentage', 2) %>%
            formatRound(column ='Abundance_percentage', 4)
        })
        
        ## Final phyloseq 
        finalPs <- reactive ({
          
          ## if no exclusive taxa
          if (input$excludeTaxa==""){
            
            ## Taxa from newly filtered table
            keeptaxa <- prevaContamFreeTable()$ASVs
            prune_taxa(keeptaxa, ps_clean())
            
          } else {
            
            removeTaxa <- reactive({
              
              ## Remove ASV based on search term
              prevaContamFreeTable() %>%
                filter(!grepl(input$excludeTaxa, get(input$taxaLevel))) %>%
                select(ASVs)
            })
            
            keepTaxa <- removeTaxa()$ASVs
            prune_taxa(keepTaxa, ps_clean())
          }
        })
        
        output$phyloSave <- renderText(
          
          paste0("The information of the final phyloseq is shown below:")
        )
        
        ## final phyloseq information
        output$saveFile <- renderPrint(
          
          finalPs()
        )
        
        # Download final phyloseq file
        output$downloadData <- downloadHandler(
          
          filename = function() {
            
            paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_filtered.rds")
          },
          content = function(file) {
            
            saveRDS(finalPs(), file)
          }
        )
        
        # Download csv abundance table
        output$downloadTable <- downloadHandler(
          
          filename = function() {
            
            paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_table.csv")
          },
          
          content = function(file) {
            # Can also re-extract from finalPs
            write.csv(prevaContamFreeTable(), file)
          }
        )
        
        ## Rhea Create taxonomy column
        ## Format data if Species exists
        if (any(colnames(tax_table(finalPs())) %in% "Species")) {
          
          taxaTable <- tax_table(finalPs())[,1:7] %>% as.data.frame() %>%
            mutate_all(tibble::lst(~str_replace(., ".__", ""))) %>%
            select(contains("str_replace")) %>% rename_all(gsub, pattern = '_.*', replacement = '') %>%
            unite(taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";",remove=TRUE)
          
        } else {
          
          ## Normal Rhea table
          taxaTable <- tax_table(finalPs())[,1:6] %>% as.data.frame() %>%
            mutate_all(tibble::lst(~str_replace(., ".__", ""))) %>%
            select(contains("str_replace")) %>% rename_all(gsub, pattern = '_.*', replacement = '') %>%
            unite(taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";",remove=TRUE)
          
        }
        
        ## Replace NA with blank  
        rheaASV <- sapply(taxaTable, gsub, pattern = "NA", replacement = "", fixed = TRUE)
        
        ## Combine ASVs to asv table
        rheaASV1 <- cbind(t(otu_table(finalPs())),rheaASV)
        rheaASV2 <- rheaASV1 %>% as.data.frame %>% rownames_to_column("#OTUId")
        
        output$downloadRheaASV <- downloadHandler(
          
          filename = function() {
            
            paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_ASV_table.tab")
          },
          
          content = function(file) {
            
            write.table(rheaASV2, file, row.names=FALSE, sep= "\t")
          }
        )

        ## MicrobiomeAmalyst parsing
        mAnalystASV <- t(otu_table(finalPs())) %>%
          as.data.frame %>% rownames_to_column("#NAME")
        row.names(mAnalystASV) <- NULL
        
        output$downloadMAnalystASV <- downloadHandler(
          
          filename = function() {
            
            paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_MA_ASV_table.txt")
          },
          
          content = function(file) {
            
            write.table(mAnalystASV,file,row.names=FALSE, sep= "\t")
          }
        )
        
        if (any(colnames(tax_table(finalPs())) %in% "Species")) {
          
          mAnalystTax <- tax_table(finalPs())[,1:7] %>% as.data.frame() %>%
            mutate_all(tibble::lst(~str_replace(., ".__", ""))) %>%
            select(contains("str_replace")) %>% rename_all(gsub, pattern = '_.*', replacement = '') %>%
            rownames_to_column("#TAXONOMY")
        } else {
          
          mAnalystTax <- tax_table(finalPs())[,1:6] %>% as.data.frame() %>%
            mutate_all(tibble::lst(~str_replace(., ".__", ""))) %>%
            select(contains("str_replace")) %>% rename_all(gsub, pattern = '_.*', replacement = '') %>%
            rownames_to_column("#TAXONOMY")
        }
        
        ## MA taxonomy table
        output$downloadMAnalystTaxa <- downloadHandler(
          
          filename = function() {
            
            paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_MA_taxonomy_table.txt")
          },
          
          content = function(file) {
            
            write.table(mAnalystTax ,file,row.names=FALSE, sep= "\t")
          }
        )
        
        output$samdat <- renderText(
          
          if ("sam_data" %in% phyloSlots()){
            
            ## Mapping file   
            rheaMapping <- sample_data(finalPs())
            rheaMapping1 <- cbind(rownames(rheaMapping), data.frame(rheaMapping, row.names=NULL))
            names(rheaMapping1)[1]<-"#SampleId"
            
            output$downloadRheaMapping <- downloadHandler(
              
              filename = function() {
                
                paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_mapping.tab")
              },
              
              content = function(file) {
                
                write.table(rheaMapping1,file,row.names=FALSE, sep= "\t")
              }
            )
            
            
            ## Turns blank to NA
            mAnalystMapping <- rheaMapping1
            names(mAnalystMapping)[1]<-"#NAME"
            
            mAnalystMapping1 <- mAnalystMapping %>% select_if(~ !any(is.na(.)))
            
            output$downloadMAnalystMapping <- downloadHandler(
              
              filename = function() {
                
                paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_MA_Metadata.txt")
              },
              
              content = function(file) {
                
                write.table(mAnalystMapping,file,row.names=FALSE, sep= "\t")
              }
            )
            
            "Can export metafile"
          } else {
            
            "Note: Missing meta data. Unable to use Rhea or MicrobiomeAnalyst !"
          }
        )
        
        output$refseq <- renderText(
          
          if ("refseq" %in% phyloSlots()){
            
            ## Fasta sequences
            seqFasta <- data.frame(refseq(finalPs()))
            seqFasta <- seqFasta %>% rownames_to_column()
            colnames(seqFasta) <- c('name','seq')
            
            #output$Rhea <- renderPrint(seqFasta)
            
            output$downloadRheaSeq <- downloadHandler(
              
              filename = function() {
                
                paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_seq.fasta")
              },
              
              content = function(file) {
                
                writeFastaFile(seqFasta, file)
              }
            )
            
            "Can export fasta sequences"
          } else {
            
            "Note: Missing fasta information (ref_seq). Unable to export this data !"
          }
        )

        output$phytree <- renderText(
          
          if ("phy_tree" %in% phyloSlots()){
            
            output$downloadTree <- downloadHandler(
              
              filename = function() {
                
                paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_tree.tre")
              },
              
              content = function(file) {
                
                ape::write.tree(phy_tree(finalPs()), file)
              }
            )
            
            "Can export phylogenetic tree"
          } else {
            
            "Note: Missing phylogenetic tree (phy_tree)). Unable to export this data !"
          }
        )
        
      } else {
        
        "Missing otu_table or tax_table. Unable to use this app !"
      }
    } else {
      
      "Upload phyloseq file to use this app !"
    }
  })
}