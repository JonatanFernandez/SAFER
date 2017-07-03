#------------------------------------------------------
#Implementation of SAFER: mass Spectrometry data Analysis by Filtering of Experimental Replicates (Zhou et al. 2016)
#
#
#Coded by:    Jonatan Fernández García        Jonifg@usal.es
#Reviewed by: Rodrigo García Valiente         Rodrigo.garcia.valiente@gmail.com
#
#Version: 1.2 2017-06-30
#
#To Do: -Further error handling. Further documentation of the output.
#------------------------------------------------------



#----------------------------------------------------------------
#Input Data Frame as in the original article for example purposes
#----------------------------------------------------------------

control1 <- c(0,0,0,0,27.5,43.5,0,0,0,0)
control2 <- c(0,0,0,0,22.8,89.3,0,0,0,0)
control3 <- c(0,0,0,3.7,26.6,82.9,13.3,0,0,0)
control4 <- c(0,0,0,11.9,48.5,115.7,11.5,19.5,0,0)
exp1 <- c(0,0,0,13.3,22.3,212.6,0,7.9,0,5)
exp2 <- c(4.6,0,29.7,16.8,59.8,148.8,45.7,16.4,14.7,131.9)
exp3 <- c(0,0,0,15.8,64.5,140.2,57.8,14.4,8.3,113.7)
exp4 <- c(9.1,2.4,47.4,8.1,96.9,176.4,34.6,11.2,4.9,44.8)
ProteinData = data.frame(control1,control2,control3,control4,exp1,exp2,exp3,exp4,
                         row.names = c("DEB1","UBC20","RPL31","HAF4","RPL7","RACK1","LARP1","ACS5","ATG3","ATG7"))
rm(control1,control2,control3,control4,exp1,exp2,exp3,exp4)

Controls <- c("control1","control2","control3","control4")
Exps <- c("exp1","exp2","exp3","exp4")
#----------------------------------------------------------------

#' SAFER FUNCTION
#'
#' @param x  <- main input. A Data Frame containing the signal intensity of peptide precursor ions belonging to a given protein in all replicates.
#' @param Controls <- An array containing the colnames for the Control Replicates.
#' @param Exps <- An array containing the colnames for the Experimental biological bait replicates .
#' @param FirstColIsRownames <- TRUE: If the first column of the dataframe contains the protein IDs. FALSE: if the Protein IDs are in the rownames (such as the example).
#' @param ExportCSV <- TRUE: export all dataframes selected in the filters in CSVs for later analysis in a spreadsheet. Generated in R working directory.
#'
#' @return A dataframe of proteins selected as potential interactive protein candidates, and three other dataframes of protein discarded in the SAFER filters,
#'         all of them together with their fold change.
#     
SAFER <- function(x, Controls, Exps, FirstColAreProteinIds = TRUE, ExportTables = TRUE,OutputName=NULL) {
  
  #-------------------------------------------------
  #                 Helper Functions 
  #-------------------------------------------------
  
  #Fold Change Calculation: Ratio of the average values, without taking zeroes into account,
  #of the MSI of the biological bait replicates and of the control replicates.
  #Return: values rounded to second significative digit sorted by Fold Change.
  
  nozeromean <- function(x) {
    if (all(x==0)) 0 else mean(x[x!=0])
  }
  
  FoldChangeCalculator <- function(x, Exps, Controls) {
    x$FC <- round(apply(x[Exps],1,nozeromean)/apply(x[Controls],1,nozeromean), 2)
    return(x[order(-x$FC),])
  }
  
  #SAFER algorithm: First step.
  #Selects all proteins that are present in at least 3 biological bait replicates
  #Return: 2 dataframes with the selected and discarded proteins.
  
  SAFERFirstFilter <- function(x, Exps) {
    ProteinsSelectedFirstFilter <- subset(x, rowSums(x[c(Exps)] != 0) >= 3)
    ProteinsDiscardedFirstFilter <- subset(x, rowSums(x[c(Exps)] != 0) < 3)
    return(list(Selected = ProteinsSelectedFirstFilter, Discarded = ProteinsDiscardedFirstFilter))
  }
  
  #SAFER algorithm: Second step.
  #Keeps all of the proteins that are absent in most control replicates (at least one nonzero value among the control replicates)
  #Return: 2 dataframes with the mostly absent and the mostly present.
  
  SAFERSecondFilter <- function(x, Controls){
    if (length(Controls)==3) {
      FilterControlMostlyAbsent <- x[rowSums(x[Controls] == 0) == 3,]
      FilterControlMostlyPresent <- x[rowSums(x[Controls] == 0) < 3,]
      return(list(MostlyAbsent = FilterControlMostlyAbsent, MostlyPresent = FilterControlMostlyPresent))
    } else if (length(Controls)==4) {
      FilterControlMostlyAbsent <- x[rowSums(x[Controls] == 0) >= 3,]
      FilterControlMostlyPresent <- x[rowSums(x[Controls] != 0) >= 2,]
      return(list(MostlyAbsent = FilterControlMostlyAbsent, MostlyPresent = FilterControlMostlyPresent))
    } else {
      print("SAFER algorith only allows as input experiments with triplicates or quadruplicates. Please check data input.")
      return(NULL)
    }
  }
  
  #SAFER algorithm: Third step.
  #Fold Change of  Proteins. Those discarded proteins with Fold Change greater than or equal to 2 are retained.
  #Also, if the retained proteins in filter two have a Fold Change lower than 2, they are discarded.
  #Return: 2 dataframes with the selected and discarded proteins.
  
  SAFERFCFilter <- function(x) {
    FoldChangeFilterSelected <- subset(x, x$FC >= 2)
    FoldChangeFilterDiscarded <- subset(x, x$FC < 2)
    return(list(Selected = FoldChangeFilterSelected,Discarded = FoldChangeFilterDiscarded))
  }
  
  #For a triplicate set of data, a supplementary filter is added:
  #Discard If the maximum value of of the three negative control replicates is higher than the minimum value of the three bait replicates
  #Return: 2 dataframes with the selected and discarded proteins.
  
  SAFERTriplicateFilter <- function(x, Controls, Exps) {
    MaxValueControls <- apply(x[Controls],1,function(x) max(x))
    MinValueExps <- apply(x[Exps],1,function(x) min(x))
    TriplicateFilterSelected <- x[MaxValueControls <= MinValueExps,]
    TriplicateFilterDiscarded <- x[MaxValueControls > MinValueExps,]
    return(list(Selected = TriplicateFilterSelected, Discarded = TriplicateFilterDiscarded))
  }
  
  #Output
  
  Output <- function(df,filename,output,OutputName=NULL) {
    write.table(df, file = paste(OutputName,"_",filename,"_",deparse(substitute(output)),".txt",sep = ""), sep = "\t", dec=".")
  }
  
  #-------------------------------------------------
  #           SAFER Core Function Execution
  #-------------------------------------------------
  #Save Dataframe name
  DataframeName <- deparse(substitute(x))
  
  #Select columns from dataset
  
  if(FirstColAreProteinIds){
    rownames(x) <- x[,1]
    x <- x[c(Controls,Exps)]
  } else {
    x <- x[c(Controls,Exps)]
  }
  
  #Execution of filters
  df_FC <- FoldChangeCalculator(x, Exps = Exps, Controls = Controls)
  FirstFiltered <- SAFERFirstFilter(df_FC, Exps)
  SecondFiltered <- SAFERSecondFilter(FirstFiltered$Selected, Controls)
  FCFiltered <- SAFERFCFilter(SecondFiltered$MostlyPresent)
  SelectedPreys <- rbind(SecondFiltered$MostlyAbsent, FCFiltered$Selected)
  
  #Selecting execution for Triplicates or Quadruplicates
  
  if (length(Controls) == 3 | length(Exps) == 3) {
    TriplicateFilter <- SAFERTriplicateFilter(SelectedPreys,Controls,Exps)
    SelectedPreys <- TriplicateFilter$Selected
    SelectedPreys = setNames(data.frame(rownames(SelectedPreys),SelectedPreys$FC),c("Selected Preys","|Fold Change"))
    FirstFilterDiscarded = setNames(data.frame(rownames(FirstFiltered$Discarded), FirstFiltered$Discarded$FC),c("Discarded: # of Exps < 3", "|Fold Change"))
    FCFilterDiscarded = setNames(data.frame(rownames(FCFiltered$Discarded),FCFiltered$Discarded$FC),c("Discarded: Fold Change < 2","|Fold Change"))
    TriplicateFilterDiscarded = setNames(data.frame(rownames(TriplicateFilter$Discarded),TriplicateFilter$Discarded$FC),c("Discarded: Max(Control) > Min(Exps)","|Fold Change"))
    if(ExportTables){
      mainDir <- getwd()
      subDir <- paste("/",paste(OutputName,DataframeName,"output", sep = "_"),sep="")
      dir.create(file.path(mainDir, subDir))
      setwd(file.path(mainDir, subDir))
      Output(df = SelectedPreys, filename = DataframeName, output = SelectedPreys,OutputName = OutputName)
      Output(df = FirstFilterDiscarded, filename = DataframeName, output = FirstFilterDiscarded,OutputName = OutputName)
      Output(df = FCFilterDiscarded, filename = DataframeName, output = FoldChangeDiscarded,OutputName = OutputName)
      Output(df = TriplicateFilterDiscarded, filename = DataframeName, output = TriplicateFilterDiscarded,OutputName = OutputName)
      setwd(mainDir)
    }
    return(list(SelectedPreys = SelectedPreys, FirstFilterDiscarded = FirstFilterDiscarded,
                FCFilterDiscarded = FCFilterDiscarded,TriplicateFilterDiscarded = TriplicateFilterDiscarded
    ))
  } else if (length(Controls) < 3 | length(Exps) < 3) {
    print("SAFER algorith only allows as input experiments with triplicates or quadruplicates. Please check data input.")
  } else {
    SelectedPreys = setNames(data.frame(rownames(SelectedPreys),SelectedPreys$FC),c("Selected Preys","|Fold Change"))
    FirstFilterDiscarded = setNames(data.frame(rownames(FirstFiltered$Discarded), FirstFiltered$Discarded$FC),c("Discarded: # of Exps < 3", "|Fold Change"))
    FCFilterDiscarded = setNames(data.frame(rownames(FCFiltered$Discarded),FCFiltered$Discarded$FC),c("Discarded: Fold Change < 2","|Fold Change"))
    if(ExportTables){
      mainDir <- getwd()
      subDir <- paste(OutputName,DataframeName,"output", sep = "_")
      dir.create(file.path(mainDir, subDir))
      setwd(file.path(mainDir, subDir))
      Output(df = SelectedPreys, filename = DataframeName, output = SelectedPreys,OutputName = OutputName)
      Output(df = FirstFilterDiscarded, filename = DataframeName, output = FirstFilterDiscarded,OutputName = OutputName)
      Output(df = FCFilterDiscarded, filename = DataframeName, output = FoldChangeDiscarded,OutputName = OutputName)
      setwd(mainDir)
    }
    return(list(SelectedPreys = SelectedPreys, 
                FirstFilterDiscarded = FirstFilterDiscarded,
                FCFilterDiscarded = FCFilterDiscarded
    ))
  }
}
  
#UsageExample:
setwd("~/R")

ExampleResults = SAFER(x=ProteinData,Controls = Controls, Exps = Exps, FirstColAreProteinIds = FALSE, ExportCSVs = TRUE)
