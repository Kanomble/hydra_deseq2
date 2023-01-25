library(tximport)
library(rhdf5)
library(DESeq2)
library(DEGreport)
#library(reshape2)
library(ggplot2)
#library(viridis)

#renames kallisto directories to samples defined in the samples.csv file
#this is important for correct file abundance file input : files <- file.path(dir, samples$ID, "abundance.h5")
filter_and_rename_kallisto_dir_names <- function(dir,samples){
  files <- list.dirs(dir)
  #reanming procedure
  for (filePath in files[2:length(files)]) {
    splitted <- strsplit(filePath,'_')
    if(file.exists(splitted[[1]][1]) == FALSE){
      print("[*] Changing Kallisto Sample Abundance Directory Name To Sample ID")
      file.rename(filePath,splitted[[1]][1])
    } 
  }
}

#reads the target transcript identifiers that have to be analyzed
#this is important for tximport 
read_target_transcript_ids <- function(target_ids_filepath){
  print("[*] Reading Target IDs")
  id_df <- read.table(target_ids_filepath,header=TRUE)
  return(id_df)
}


#importing abundance files into deseq2 and performing dese2 analysis
deseq2_pipeline <- function(dir,samples,tx2gene,threshold=10,reference="control",ref=TRUE){
  cat("[*] Starting DESeq2 Pipeline\n")
  #tximport file input
  files <- file.path(dir, samples$ID, "abundance.h5")
  names(files) <- samples$ID
  #checks if files exist
  all(file.exists(files))
  cat("[*] Reading Files With tximport\n")
  #reading files with tximport, assigning TXNAME and GENEID to the same transcript identifiers
  txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  #replace ID to treatment for deseq2 analysis
  sampleTable <- data.frame(condition = factor(samples$treatment))
  rownames(sampleTable) <- colnames(txi.kallisto.tsv$counts)
  cat("[*] Reading tximport Dataframe with DESeqDataSetFromTximport\n")
  #data input to deseq2
  dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable, ~condition)
  #pre filtering
  keep <- rowSums(counts(dds)) >= threshold
  dds <- dds[keep,]
  cat("[*] Starting DESeq Analysis\n")
  #deseq2 analysis
  if(ref==TRUE){
    cat("[*] Releveling To Real Reference\n")
    dds$condition <- relevel(dds$condition, ref = reference)
  }
  dds <- DESeq(dds)
  cat("[+] DONE\n")
  return(dds)
}