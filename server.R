library("shiny")
library("tidyverse")
library("phyloseq")
library("BiocManager")
library("DT")
library("methods")
library("ape")
library("Biostrings")

options(repos = BiocManager::repositories())

make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
  is.na(x) <- x=="NA"; x} else {
    x}

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


function(input, output, session) {
  
  ## Load in phyloseq
  ps <- reactive({
    
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == c("rds","RDS"), "Please upload a rds file"))
    
    readRDS(file$datapath)
    
  })
  
  ## Phyloseq information
  output$inPs <- renderPrint({
    
    ## phyloseq file loaded
    if (!is.null(ps())){
      
      output$summary <- renderText(
        
        ## Check if file is phyloseq
        if (methods::is(ps(), "phyloseq")){
          
          ## Check phyloseq slot names
          phyloSlots <- reactive(getslots.phyloseq(ps()))
          
          ## Slot require for full analysis
          mustHave <- c("otu_table", "tax_table", "sam_data", "phy_tree", "refseq")
          
          output$phyloseq <- renderText(
            
            ## Check phyloseq has required slots
            if (any(mustHave %in% phyloSlots())){
              
              ## Compute number of sample for each ASV, store as data.frame
              abunfltdf <- reactive({
                
                apply(X = otu_table(ps()),
                      MARGIN = ifelse(taxa_are_rows(ps()), yes = 1, no = 2),
                      FUN = function(x){sum(x > 0)})
              }) 
              
              ## Add taxonomy and total read counts to this data.frame
              abunfltdata <- reactive({
                
                data.frame(Prevalence = abunfltdf(),
                           TotalAbundance = taxa_sums(ps()),
                           tax_table(ps())[, input$taxaLevel])
              })
              
              totalReads <- reactive({ sum(abunfltdata()$TotalAbundance) })

              output$taxaFilter <- renderPlot({
                
                ggplot(abunfltdata(), aes(TotalAbundance, Prevalence / nsamples(ps()),color=get(input$taxaLevel))) +
                  geom_hline(yintercept = input$prevaTaxa/100, alpha = 1, linetype = 2) +
                  geom_vline(xintercept = input$abunTaxa * totalReads() / 100, alpha = 1, linetype = 3) +
                  geom_point(size = input$pointSize, alpha = input$pointTrans) +
                  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
                  facet_wrap(~get(input$taxaLevel), ncol = 2) + theme(legend.position="none")
                ## input$taxaLevel is a constant characters string, so need to have 'get'.
              })
              
              prevaContamFreeTable <- reactive(
                
                abunfltdata() %>% rownames_to_column("ASVs") %>% 
                  mutate(Prevalence_percentage = Prevalence/nsamples(ps()) * 100, Abundance_percentage = TotalAbundance/totalReads()*100) %>%
                  select(ASVs,input$taxaLevel,Prevalence,Prevalence_percentage,TotalAbundance,Abundance_percentage) %>%
                  filter(Prevalence_percentage >= input$prevaTaxa & Abundance_percentage >= input$abunTaxa)
              )
              
              output$filterSummary <- renderText(
                
                paste0("ASVs with number of reads below ", ceiling(input$abunTaxa * totalReads() / 100), " (", input$abunTaxa, " %) and appeared in less than " , ceiling(input$prevaTaxa * nsamples(ps()) / 100), " samples (", input$prevaTaxa, " %)  were filtered out. There are ", nrow(prevaContamFreeTable ()), " ASVs left (", round((ntaxa(ps()) - nrow(prevaContamFreeTable ())) / ntaxa(ps()) *100, digits = 2), " % was removed)." )
              )
              
              output$taxaTable <- DT::renderDataTable({
                
                DT::datatable(prevaContamFreeTable(), options = list(orderClasses = TRUE)) %>%
                  formatRound(columns='Prevalence_percentage', 2) %>%
                  formatRound(column ='Abundance_percentage', 4)
                
              })
              
              saveData <- reactive ({
                
                if (input$excludeTaxa==""){
                  
                  keeptaxa <- prevaContamFreeTable()$ASVs
                  prune_taxa(keeptaxa, ps())
                  
                } else {
                  
                  removeTaxa <- reactive({
                    
                    prevaContamFreeTable() %>%
                      filter(!grepl(input$excludeTaxa, get(input$taxaLevel))) %>%
                      select(ASVs)
                  })
                    
                  keepTaxa <- removeTaxa()$ASVs
                  prune_taxa(keepTaxa, ps())
                }
              })
              
              output$phyloSave <- renderText(
                
                paste0("The information of the final phyloseq is shown below:")
              )
              
              output$saveFile <- renderPrint(
                
                saveData()
              )

              # Download csv abundance table
              output$downloadTable <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_table.csv")
                },
                
                content = function(file) {
                  
                  write.csv(prevaContamFreeTable(), file)
                }
              )
              
              # Download filtered files
              output$downloadData <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_filtered.rds")
                },
                content = function(file) {
                  
                  saveRDS(saveData(), file)
                }
              )
              
              output$downloadTree <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_tree.tre")
                },
                
                content = function(file) {
                  
                  ape::write.tree(phy_tree(saveData()), file)
                }
              )
              
              ## Rhea Create taxonomy column
              taxaTable <- tax_table(saveData())[,1:6] %>% as.data.frame() %>%
                mutate_all(tibble::lst(~str_replace(., ".__", ""))) %>%
                select(contains("str_replace")) %>% rename_all(gsub, pattern = '_.*', replacement = '') %>%
                unite(taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus"), sep=";",remove=TRUE)

              ## Replace NA to blank  
              rheaASV <- sapply(taxaTable, gsub, pattern = "NA", replacement = "", fixed = TRUE)
              ## Combine ASVs to asv table
              rheaASV1 <- cbind(t(otu_table(saveData())),rheaASV)
              rheaASV2 <- rheaASV1 %>% as.data.frame %>% rownames_to_column("#OTUId")

#              output$Rhea <- renderPrint(rheaASV2)
              
              ## Mapping file
              rheaMapping <- sample_data(saveData())
              rheaMapping1 <- cbind(rownames(rheaMapping), data.frame(rheaMapping, row.names=NULL))
              names(rheaMapping1)[1]<-"#SampleId"
              
              ## Fasta sequences
              seqFasta <- data.frame(refseq(saveData()))
              seqFasta <- seqFasta %>% rownames_to_column()
              colnames(seqFasta) <- c('name','seq')

              output$downloadRheaASV <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_ASV_table.tab")
                },
                
                content = function(file) {
                  
                  write.table(rheaASV2, file, row.names=FALSE, sep= "\t")
                }
              )
              
              output$downloadRheaMapping <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_mapping.tab")
                },
                
                content = function(file) {
                  
                  write.table(rheaMapping1,file,row.names=FALSE, sep= "\t")
                }
              )
              
              output$downloadRheaSeq <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_Rhea_seq.fasta")
                },
                
                content = function(file) {
                  
                  writeFastaFile(seqFasta, file)
                }
              )
              
              ## MicrobiomeAmalyst parsing
              mAnalystASV <- t(otu_table(saveData())) %>%
                as.data.frame %>% rownames_to_column("#NAME")
              row.names(mAnalystASV) <- NULL
              
              mAnalystTax <- tax_table(saveData())[,1:6] %>% as.data.frame() %>%
                mutate_all(tibble::lst(~str_replace(., ".__", ""))) %>%
                select(contains("str_replace")) %>% rename_all(gsub, pattern = '_.*', replacement = '') %>%
                rownames_to_column("#TAXONOMY")

              
              ## Turns blank to NA
              mAnalystMapping <- rheaMapping1
              names(mAnalystMapping)[1]<-"#NAME"
              
              mAnalystMapping1 <- mAnalystMapping %>% select_if(~ !any(is.na(.)))

              output$downloadMAnalystASV <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_MA_ASV_table.txt")
                },
                
                content = function(file) {
                  
                  write.table(mAnalystASV,file,row.names=FALSE, sep= "\t")
                }
              )
              
              output$downloadMAnalystTaxa <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_MA_taxa_table.txt")
                },
                
                content = function(file) {
                  
                  write.table(mAnalystTax ,file,row.names=FALSE, sep= "\t")
                }
              )
              
              output$downloadMAnalystMapping <- downloadHandler(
                
                filename = function() {
                  
                  paste0(input$name, "_p", input$prevaTaxa, "_a", input$abunTaxa, "_MA_mapping.txt")
                },
                
                content = function(file) {
                  
                  write.table(mAnalystMapping,file,row.names=FALSE, sep= "\t")
                }
              )
              
              ###output$unique <- renderPrint(
              output$unique <- renderText(
                
                paste0("An online shiny app (https://giangle.shinyapps.io/phyloFilter/) was used to filter ASVs, which presented in less than ",
                       ceiling(input$abunTaxa * totalReads() / 100), " (", input$abunTaxa, " %) reads and appeared in less than " ,
                       ceiling(input$prevaTaxa * nsamples(ps()) / 100), " samples (", input$prevaTaxa, " %). A total of ",
                       nrow(prevaContamFreeTable ()), " ASVs were maintained for downstream analysis.",
                       "There were ", length(get_taxa_unique(saveData(), taxonomic.rank = "Genus")), " unique genus and ",
                       length(get_taxa_unique(saveData(), taxonomic.rank = "Species")), " unique species in the final phyloseq." )
              )
              
              ## Text need to be at the end so it does not disrupt the flow
              paste0("The input file has ", nsamples(ps()), " samples, with ", totalReads(), " reads, across ", ntaxa(ps()), " ASVs")
              
            } else {
            
              "Please recheck to make sure that your phyloseq contains the following: otu_table, tax_table, sam_data, phy_tree and/or refseq"
            })
          
          "Phyloseq file: "
          
        } else {
          
          "Please upload a phyloseq file !"
        } 
      )
    }
  })
}