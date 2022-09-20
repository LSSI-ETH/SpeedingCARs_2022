
###################################
# title: "SpeedingCARs: Demultiplexing of SpeedingCARs single cell CAR sequencing data"
# date: 09.09.22
# author: Rocio Castellanos Rueda
#######################################

library(Biostrings)
library(tidyverse)
library(dbplyr)
library(reshape2)

setwd(" --- ")

#######################################
### Load sample names

sample <- c("D1_2DomLib", 
            "D1_2DomLib_WT", 
            "D2_2DomLib_WT",
            "D3_2DomLib_1", 
            "D3_2DomLib_2", 
            "D3_2DomLib_WT")
BC_list <- c('D1_2DomLib_barcodes.tsv', 
             'D1_2DomLib_WT_barcodes.tsv', 
             'D2_2DomLib_WT_barcodes.tsv',
             'D3_2DomLib_1_barcodes.tsv', 
             'D3_2DomLib_2_barcodes.tsv', 
             'D3_2DomLib_WT_barcodes.tsv')
seqfile <- c("D1_2DomLib.fq", 
            "D1_2DomLib_WT.fq", 
            "D2_2DomLib_WT.fq",
            "D3_2DomLib_1.fq", 
            "D3_2DomLib_2.fq", 
            "D3_2DomLib_WT.fq")


for (j in 1:length(seqfile)) {
  
  #######################################
  ### 1. Load raw sequences
  
  # Read Hifi PacBio sequences
  
  rawseqs <- readDNAStringSet(seqfile[j] , format="fastq")
  
  #######################################
  ### 2 identify F and R reads based on common adoaptor sequence and reverse complement the R reads
  
  T.E_end_R <- DNAString("CTACACGACGCTCTTCCGAT")
  T.E_end_F <- DNAString("ATCGGAAGAGCGTCGTGTAG")
  
  # reverse complement the sequences that have a CD3z sequence in reverse
  rawseqs[(which(vcountPattern(T.E_end_R,rawseqs) == 1))] <- reverseComplement(rawseqs[(which(vcountPattern(T.E_end_R,rawseqs) == 1))])
  
  #######################################
  ### 3. filter reads that have a CAR gene
  
  CD28TM_F <- DNAString("ATTGAAGTTATGTAT")
  CAR_rawseq <- rawseqs[(which(vcountPattern(CD28TM_F,rawseqs) == 1))]
  
  #######################################
  ### 4. Identify 10X BC and CAR variant for each read
  
  # Identify 10X BC sequence
  after_adaptor <- DNAString("AGATCGGAAGAGCGT")  # sdquence downstream of 16 nt 10X BC
  CAR_rawseq <- CAR_rawseq[(which(vcountPattern(after_adaptor,CAR_rawseq) == 1))]
  hits <- vmatchPattern(after_adaptor,CAR_rawseq)
  
  #Add the 16 nt barcode to the selected sequence
  hits@width0 <- as.integer(hits@width0 + rep(c(16), length = length(CAR_rawseq)))  
  
  #extract barcode strings from stringset and reverse complement to match 10X barcodes
  barcodes <- reverseComplement(unlist(DNAStringSetList(extractAt(CAR_rawseq, hits))) %>% subseq(end= 16))
  
  # Read a list of possible CAR A/B_linker domains (13ntA/6ntlink/12ntB) in FASTA format
  linker_lib <- readDNAStringSet("2DomLinkerLib.fasta")
  
  # identify CAR variant
  Variant <-sapply(vwhichPDict(linker_lib,CAR_rawseq,max.mismatch=1), `[`, 1)
  # convert number to variant name
  Variant_names <- names(linker_lib)[match(Variant, (1:180), nomatch = NA)]
  
  #######################################
  ### 5. merge Variant and Barcode information into a dataframe
  
  seq_names <- names(barcodes) # or seq_names <- names(CAR_rawseq) , they should be the same
  BC_sequence <- paste(barcodes)
  annotation_df <- data.frame(seq_names, BC_DNA = as.character(BC_sequence), Variant_names)
  annotation_df <- filter(annotation_df ,!is.na(annotation_df$Variant_names))
  
  write.csv(annotation_df, paste("T.E.raw_", sample[j], ".csv" ,sep="")) 
  
  #######################################
  ### 6. CAR variant assignment to list of 10X barcodes 
  
  all_10x_bc <- read.table(file = BC_list[j], sep = '\t')
  all_10x_bc$barcodes <- substr(all_10x_bc$V1, 1, 16)
  
  Ture_BC_df <- filter(annotation_df , BC_DNA %in% all_10x_bc$barcodes)
  
  # Extract variant - barcode combinations and the read counts of each
  count_matrix <- with(Ture_BC_df, table(Variant_names, BC_DNA))  #colnames(count_matrix)
  temp <- reshape2::melt(count_matrix) %>%
    filter(value != 0)
  nrow(temp)
  
  # Due to sequencing mistakes often more than one variant is assigned to the same barcode --> Need a filtering step:
  temp$unique <- c(NA)
  
  nest <- dplyr::group_by(temp, BC_DNA)%>% 
    nest(.key = "BC_DNA")
  
  #filtering of sequencing mistakes: (depending on sequencing depth filtering criteria should be tuned)
  
  for (i in 1:lengths(nest[1])) {
    
    a <- nest[[2]][[i]]
    if (nrow(a) != 1) {
      a <- arrange(a, desc(value)) # order by occurence
      
      if (a$value[1]/a$value[2] > 9) {    # if there is a difference of 10 fold we accept 
        nest[[2]][[i]] <- a[1,]
      }
      else {
        a$unique<- "D"
        nest[[2]][[i]] <- a
        if (a$value[1] == 1) {
          a$unique<- "Error"      ## if only 1 read is detected this is identified as an error
          nest[[2]][[i]] <- a
        }
      }
    }
  }
  
  
  # Unnest data frame:
  annotation_df <- data.frame(unnest(nest, cols = colnames(nest[2])) )
  annotation_df <- filter(annotation_df , is.na(annotation_df $unique) | unique != "Error")
  mean(annotation_df$value)
  
  #######################################
  ### 8. create .csv to use in seurat
  
  barcodes <- select(all_10x_bc, -V1)
  colnames(barcodes) <- "barcodes_assigned"
  barcodes<- cbind(barcodes,CAR = NA)
  
  for (i in 1:nrow(annotation_df)){
    read<-DNAString(annotation_df[i,1])
    n <- which(vcountPattern(read,barcodes$barcodes_assigned,max.mismatch=0) == 1)
    barcodes$CAR[n] <- as.character(annotation_df[i,2])
    if (!is.na(annotation_df$unique[i])) {
      barcodes$CAR[n] <- "D"
    }
  }
  
  
  name <- paste(sample[j],"_barcodes_assigned.csv", sep="")
  write.csv(barcodes, name)
  
}







