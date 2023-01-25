options(stringsAsFactors=FALSE)

gits <- "/home/rstudio/scripts/"
gitp <- "/home/rstudio/data/"
datp <- "/home/rstudio/results/"
resp <- "/home/rstudio/results/"

source(file.path(gits,"deseq2_analysis_hydra.R"))

## create output path for results
outp <- file.path(datp,"deseq2")
dir.create(outp, showWarnings=FALSE)


#### DESEQ2 ####

### RUN PARAMETERS
## threshold for pre-filtering
threshold <- 0
## mapping ID: use raw or deduplicated?
MID <- "kallistoCountsDedup" # "kallistoCounts" #


## get common gene IDs
target_ids <- file.path(gitp,"result_tables","EcoKD1",MID,
                        "I25727-L1","abundance.tsv")
id_df <- read.table(target_ids, header=TRUE)
tx2gene <- data.frame(TXNAME=id_df$target_id,GENEID=id_df$target_id)

## get experiment sample IDs
experiments <- read.csv(file.path(gitp,"sample_descriptions.csv"))

## get control conditions for each experiment
controls <- read.csv(file.path(gitp, "sample_controls.csv"))

### PERFORM DeSeq2 FOR ALL EXPERIMENTS
for ( exp in unique(experiments$EXP) ) {
  
  cat(paste("[*] Handling experiment", exp, "\n"))
  
  ## results directory
  dir <- file.path(gitp,exp,MID)
  
  ## sample IDs
  samples <- experiments[experiments$EXP==exp,-1]
  
  if ( exp%in%c("HydraAHL","HydraRecolonization") ) 
    samples$ID <- sub("-S1","-L1",samples$ID)
  if ( exp%in%"HydraTemperature" )
    samples$ID <- paste0(samples$ID,"-L1")
  
  ## read raw read-counts
  dds <- deseq2_pipeline(dir = dir,samples = samples, tx2gene = tx2gene,
                         ref = FALSE, threshold = threshold)
  
  ## PCA Plot
  
  rld <- rlog(dds, blind=TRUE)
  cat(paste("[*] PRODUCING PCA",file.path(outp,paste0(exp,"_pca.png")),"\n"))
  png(file=file.path(outp,paste0(exp,"_pca.png")), width=800, height=550)
  pca <- plotPCA(rld, intgroup="condition")
  print(pca)
  dev.off()
  
  cat("[*] DONE\n")
  ## perform differential expression analysis for each control group
  
  ## control groups
  cntrl <- controls[controls$EXP==exp,-1]
  
  ## obtaining results
  for ( i in 1:nrow(cntrl) ) {
    res <- results(dds, contrast=c("condition",
                                   cntrl$treatment.group[i],
                                   cntrl$control.group[i]))
    resOrdered <- res[order(res$pvalue),]
    
    ## writing output file
    write.csv(as.data.frame(resOrdered), 
              file=file.path(outp,paste0(exp,"_",
                                         cntrl$treatment.group[i],"_vs_",
                                         cntrl$control.group[i],
                                         "_results.csv")))
  

    ## volcano plot
    
    cat(paste("[*] PRODUCING VOLCANO PLOT",file.path(outp,paste0(exp,"_",
                                                                    cntrl$treatment.group[i],"_vs_",
                                                                    cntrl$control.group[i],
                                                                    "_volcano_plot.png")),"\n"))
    resOrdered[["id"]] <- row.names(resOrdered)
    show <- as.data.frame(resOrdered[1:10, c("log2FoldChange", "padj","id")])#"id"
    png(file=file.path(outp,paste0(exp,"_",
                                   cntrl$treatment.group[i],"_vs_",
                                   cntrl$control.group[i],
                                   "_volcano_plot.png")),width=800, height=550)
    print(degVolcano(resOrdered[,c("log2FoldChange", "padj")],plot_text = show))
    cat("[*] DONE\n")
    dev.off()
  }
}


