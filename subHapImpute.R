#!/usr/bin/env Rscript

# Script for fastPHASE input generation
suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(tidyverse)))

args <- commandArgs(TRUE)
# transposed matrix (marker as row)
FILE <- args[1]


# make vector with major or minor allel 
f.allelCount <- function(x, major = T) {
  # major=F for minor allele
  String <- paste(na.omit(unlist(x)),collapse="")
  # counting
  countTable <- as.data.frame(matrix(c('A', 'T', 'C', 'G'), ncol = 1))
  countTable$count <- NA
  for (LET in c('A', 'T', 'C', 'G')){
    countTable[countTable$V1 == LET, ]$count <- str_count(String, pattern = LET)
  }
  # major or minor?
  rankedNt <- countTable[order(countTable$count, decreasing = T), 1]
  if(major == T){
    OUT <- rankedNt[1]
  }else{
    OUT <- rankedNt[2]
  }
  return(as.character(OUT))
}


## MATRIX
matrice <- fread(FILE, data.table = F, header = T) 
matrice <-matrice[, -c(2:11)]

matrice[matrice == 'NN'] <- NA
colnames(matrice)[1] <- 'Marker'

# duplicated names
colnames(matrice) <- make.unique(colnames(matrice))

# get major allele in a vector
MAJOR <- as.vector(unlist(apply(matrice[,-1], 1, f.allelCount)))


# get major allele in a vector
MINOR <- as.vector(unlist(apply(matrice[,-1], 1, f.allelCount, major = F)))


# OUT header
NRindividuals <- ncol(matrice)-1
NRmarker <- nrow(matrice)
cat(paste0(NRindividuals, "\n"))
cat(paste0(NRmarker, "\n"))



################################################


# format genotyping
for (ID in colnames(matrice[,-1])){
# run each genotype
  cat(paste0('# ', ID, "\n"))
  ROW <- matrice[,ID]
  
  upperLine <- c()
  bottomLine <- c()
  for (markerPOS in 1:length(MAJOR)){
    genotyping <- ROW[markerPOS]
    # split alleles
    alleles <- sort(unlist(strsplit(genotyping, split = '')))
    majorREF <- MAJOR[markerPOS]
    minorREF <- MINOR[markerPOS]
    # test if NA 
    if (is.na(genotyping)){
      upperLine <- c(upperLine, '?')
      bottomLine <- c(bottomLine, '?')
    } else if (alleles[1] %in% c(majorREF, minorREF) & alleles[1] %in% c(majorREF, minorREF)) {   
      # test hetero 
      if (length(unique(alleles)) == 2){
        upperLine <- c(upperLine, 0)
        bottomLine <- c(bottomLine, 1)
      } else if (alleles[1] == majorREF & alleles[2] == majorREF) {
        upperLine <- c(upperLine, 0)
        bottomLine <- c(bottomLine, 0)
        #cat(paste0('CUCU con: ', markerPOS))
      } else if (alleles[1] == minorREF & alleles[2] == minorREF) {
        upperLine <- c(upperLine, 1)
        bottomLine <- c(bottomLine, 1)
      }
      # if other alleles
    } else {
      upperLine <- c(upperLine, '?')
      bottomLine <- c(bottomLine, '?')
    }
  }  
  
  cat(paste0(cat(upperLine), "\n"))
  cat(paste0(cat(bottomLine), "\n"))
}

