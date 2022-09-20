
###################################
# title: "SpeedingCARs: Analysis of CAR library diversity"
# date: 09.09.22
# author: Rocio Castellanos Rueda
#######################################

library(Biostrings)
library(tidyverse)
library(dbplyr)

#######################################
### 1. Load and filter reads

# Read a list of CAR sequences in FASTA format 
seqfile <- "HER2_2DomLib.fq"
HER2_rawseq <- readDNAStringSet(seqfile, format="fastq")

#provide primer sequences
CD3_end_F <- DNAString("CTTCACATGCAGGCCCTGCCCCCTCGCTAA")
CD3_end_R <- DNAString("TTAGCGAGGGGGCAGGGCCTGCATGTGAAG")

# reverse complement the sequences that have a CD3z sequence in reverse
rawseqs[(which(vcountPattern(CD3_end_R,rawseqs) == 1))] <- reverseComplement(rawseqs[(which(vcountPattern(CD3_end_R,rawseqs) == 1))])

# extract the sequences that have a CD3z sequence
HER2_rawseq <- rawseqs[(which(vcountPattern(CD3_end_F,rawseqs) == 1))]

#######################################
### 2. Identify CAR Variant for each read

seqfile <- "HER2_2DomLib.fq"
HER2_rawseq <- readDNAStringSet(seqfile, format="fastq")

# Read a list of possible CAR A/B_linker domains (8ntA/6ntlink/8ntB) in FASTA format
linker_lib <- readDNAStringSet("2DomLinkerLib.fasta")

#create a dataframe 
lib_name <- names(linker_lib)
lib_sequence <- paste(linker_lib)
linker_df <- data.frame(lib_name, DNA = as.character(lib_sequence))
linker_df$count <- c(NA)

for (i in 1:nrow(linker_df)) {
  linker_df$count[i] <- sum(vcountPattern(DNAString(linker_df[i,2]),HER2_rawseq))
}


linker_df$lib_name <- gsub("4-1BB","41BB", linker_df$lib_name)
linker_df$lib_name <- gsub(" ","", linker_df$lib_name)
linker_df <- data.frame(linker_df, str_split_fixed(linker_df$lib_name, "-", 2))

m <- spread(select(linker_df,count,X1,X2), X2, count) 
rownames(m) <- m[,1]
m <- select(m, -X1)

#order matrix
col.order <-  c("CD3Z",	"LMP2",	"K1",	"FCGR2A",	"FCER1G",	"DAP12",	"CD79B",	"CD79A",	"CD3G",	"CD3E",	"CD3D",	"GP")
row.order <- c("CD28", "41BB", "CD4", "CD150", "CD226", "HVEM", "CD30", "CD84", "CD357", "CD244", "FCRL6","FCRL1",
               "CD223", "DR3","TIM1")

m <- m[row.order,col.order]
#write.csv(m, "PacBio_frequencies.csv")

library(circlize)
library("viridis") 

grid.col <- c((plasma(15)),(plasma(15)))
names( grid.col ) =  c("CD28", "41BB", "CD4", "CD150", "CD226", "HVEM", "CD30", "CD84", "CD357", "CD244", "FCRL6","FCRL1",
                       "CD223", "DR3","TIM1","CD3Z",	"LMP2",	"K1",	"FCGR2A",	"FCER1G",	"DAP12",	"CD79B",	"CD79A",	"CD3G",	"CD3E",	"CD3D",	"GP")
chordDiagram(t(m),grid.col = grid.col)
#chordDiagram(m)

pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE) # save 4x5

# frequencies
pheatmap((m/sum(m))*100, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = round((m/sum(m))*100,2)) # save 4x5













