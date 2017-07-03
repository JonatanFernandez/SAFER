#---------------------------------------------------------------
#  Extractor from MaxQuant Output "Protein_Groups" to SAFER input
#
#Coded by:    Jonatan Fernandez Garcia        Jonifg@usal.es
#Reviewed by: Rodrigo Garcia Valiente         Rodrigo.garcia.valiente@gmail.com
#
#Version: 1. 2017-06-26
#
#---------------------------------------------------------------

#------------------------
#         Data
#------------------------

ProteinGroups = read.csv(file.choose(), header = TRUE, sep = "\t", dec = ".", stringsAsFactors = FALSE)


#-----Current Example--------
setwd("~/R")
ProteinGroups = read.csv("proteinGroups.txt", header = TRUE, sep = "\t", dec = ".", stringsAsFactors = FALSE)
#----------------------------




#' MQ_ProteinGroups_2_SAFER
#'
#' @param df  <- main input. A Data Frame created by the previous function reading the proteinGroups.txt output of MaxQuant
#' @param Intensities. Choose between "LFQ", "iBAQ" or "Intensities".
#' @param CountsCutoff <- Discards ProteinIds if Unique or Razor peptide coutns are lower than the Cutoff
#' @param RCBS <- If true (by Default), applies the RCBSFilter: Discard Proteins if "Reverse", "Potential.contaminant" or "Only.identified.by.site" are positive.
#'
#' @return A dataframe of Proteins IDs and related intensities suited for latter analysis with SAFER.
#   

MQ_ProteinGroups_2_SAFER <- function(df, Intensities = "LFQ", CountsCutoff=2, RCBS = TRUE) {
  
  #-----------------------
  #   Helper Functions
  #-----------------------
  
  #RCBSFilter: Discarded if "Reverse", "Potential.contaminant" or "Only.identified.by.site" are positive.
  if (RCBS) {
    RCBSFilter <- function(df) {
      return(subset(df, ProteinGroups$Reverse != "+" 
                    & ProteinGroups$Potential.contaminant != "+" 
                    & ProteinGroups$Only.identified.by.site != "+"))
    }
  }
  
  
  #PeptideCoutnsFilter: Discards ProteinIds if Unique or Razor peptide coutns are lower than the Cutoff (n=2 by default, defined in Corefunction)
  
  PeptideCountsFilter <- function(df, CountsCutoff=CountsCutoff) {
    CountsCutoff = as.integer(CountsCutoff)
    if (CountsCutoff == 0) {
      dfFiltered = df
    } else if (CountsCutoff >= 1) {
      df$Peptide.counts..razor.unique.List = strsplit(df$Peptide.counts..razor.unique.,";")
      df$Protein.IDs.List = strsplit(df$Protein.IDs,";")
      if(all(lengths(df$Peptide.counts..razor.unique.List) == lengths(df$Protein.IDs.List))){
        BooleanList = lapply(df$Peptide.counts..razor.unique.List, function(x) x >= CountsCutoff)
        for (i in 1:length(BooleanList)) {
          df$Protein.IDs.List[[i]] = df$Protein.IDs.List[[i]][BooleanList[[i]]]
          df$Peptide.counts..razor.unique.List[[i]] = df$Peptide.counts..razor.unique.List[[i]][BooleanList[[i]]]
        }
        
        n<-1
        deleteindexes <- list()
        
        for (i in 1:length(df$Protein.IDs.List)){
          if (identical(df$Protein.IDs.List[[i]], character(0))){
            deleteindexes[n] <- i
            n <- n+1
          }
        }
        
        dfFiltered <- df[-unlist(deleteindexes),]
        dfFiltered$Protein.IDs <- sapply(dfFiltered$Protein.IDs.List ,function(x) paste(as.character(x), collapse=";"))
        dfFiltered$Peptide.counts..razor.unique. <- sapply(dfFiltered$Peptide.counts..razor.unique.List ,function(x) paste(as.character(x), collapse=";"))
        DropColumns <- c("Protein.IDs.List","Peptide.counts..razor.unique.List")
        dfFiltered <- dfFiltered[,!(names(dfFiltered) %in% DropColumns)]
      }else{
        return(print("Error, mismatch with peptide counts and Protein IDs, please check"))
      }
    } else {
      print("CountsCutOff Must be a positive integer")
      return(NULL)
    }
    return(dfFiltered)
  }
  
  #Intensities extractor functions:
  
  ExtractLFQIntensities <- function(x) {
    return(cbind(x$Protein.IDs, x[,startsWith(colnames(x), "LFQ")]/10^5))
  }
  
  ExtractiBAQIntensities <- function(x) {
    return(cbind(x$Protein.IDs, x[,startsWith(colnames(x), "iBAQ.")]/10^5))
  }
  
  ExtractIntensities <- function(x) {
    return(cbind(x$Protein.IDs, x[,startsWith(colnames(x), "Intensity.")]/10^5))
  }
  
  #------------------------
  #   Core execution
  #------------------------
  originalEntries <- length(unlist(df$Protein.IDs))
  if(RCBS){
    df <- RCBSFilter(df)
  }
  df <- PeptideCountsFilter(df, CountsCutoff = CountsCutoff)
  EntriesEliminated = originalEntries - length(unlist(df$Protein.IDs))
  print(paste(EntriesEliminated," entries were eliminated from the dataset."))
  switch(Intensities, 
         "LFQ" = {df = ExtractLFQIntensities(df)},
         "iBAQ" = {df = ExtractiBAQIntensities(df)},
         "Intensities" = {df = ExtractIntensities(df)}
  )
}

LFQ2 <- MQ_ProteinGroups_2_SAFER(ProteinGroups,Intensities = "LFQ", CountsCutoff = 2)
LFQ4 <- MQ_ProteinGroups_2_SAFER(ProteinGroups,Intensities = "LFQ", CountsCutoff = 4)
LFQ6 <- MQ_ProteinGroups_2_SAFER(ProteinGroups,Intensities = "LFQ", CountsCutoff = 6)
iBAQ2 <- MQ_ProteinGroups_2_SAFER(ProteinGroups,Intensities = "iBAQ", CountsCutoff = 2)
INT2 <- MQ_ProteinGroups_2_SAFER(ProteinGroups,Intensities = "Intensities", CountsCutoff = 2)
