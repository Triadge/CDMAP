## Context Mutation Analysis Software Package version 1.0.1

# Author: David Logan Patton
# Contributors: Way Sung (Principal Investigator), Thomas Cardenas (Web Visualization), 
# Perrin Mele (Theta-W Calculations, Web Visualization), Jon Navarro (Web Visualization)

# Purpose: This is the Main executable script for the CDMAP Package that runs the Single Organism Analysis Pipeline.
# In this script, CDMAP first checks for, and if necessary installs all required packages. Then CDMAP prompts the user
# for relevant input information such as the cleaned VCF file, GBK file, and reference FASTA file and then dynamically
# creates all the relevant child directories needed to store all output files. Then CDMAP proceeds to calculate the counts
# and rates for all 64 site-specific contexts, codon usage, both with respect to the leading and lagging strand.



### INSTALLATION AND ATTACHMENT OF REQUIRED PACKAGES ###
#========================================#
packages <- c("seqinr", "BiocManager", "pracma", "beepr", "lattice", "tidyverse", "vcfR", "stringr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
BiocManager::install("genbankr")
library("genbankr") # BiocManager package for genbank file parsing and manipulation

username <- Sys.info()[7]
User <- paste("/User/", username, sep ="")
MainDir <- paste("/Users/", username, "/Desktop/CDMAP", sep = "")
LibDir <- paste("/Users/", username, "/Desktop/CDMAP/CDMAP_Library", sep = "")

#========================================#



### USER INPUT OF REQUIRED FILES, MA INFORMATION AND POINTERS FOR OUTPUT DIRECTORIES ###
#========================================#
#name of the organism for file naming purposes
organism <- readline("What is the name of the organism (and chromosome number)?")
BaseCheck <- c("basecall", "base call", "Base call","Base Call")
VCFCheck <- c("vcf", "Vcf", "VCF")
VcfOrBaseCall <- readline("Do you have a VCF file or a Mutation Base Call file?")

#/Users/triadge/Desktop/CDMAP_Release_v1/CDMAP_Library

#Section 1: intaking the data
DirCheck <- readline("Do you use Default or Customized Directories for CDMAP?")
DefaultCheck <- c("Default", "default", "Def", "def")
CustomCheck <- c("Custom", "custom", "Customized", "customized")
MasterCheck <- c("Triadge", "triadge")

if(any(grepl(DirCheck, DefaultCheck, ignore.case = TRUE)))
{
  setwd(LibDir)
  source("DirectoryCheck.r")
}

if(any(grepl(DirCheck, MasterCheck, ignore.case = TRUE)))
{
  setwd(LibDir)
  source("DirectoryCheck.r")
}

if(any(grepl(DirCheck, DefaultCheck, ignore.case = TRUE)))
{
  Path_to_scripts <- LibDir
  Path_wd <- MainDir
  Path_output <- readline("Where do you want to output your data?")
  Path_RefFile <- readline("What is your reference sequence? (please provide the full Path)")
  Path_GBFile <- readline("What is your Genbank file? (please provide the full Path)")
  Path_correlate_repo <- readline("Where would you like to store multi-organism output?")
  
  if(any(grepl(VcfOrBaseCall, VCFCheck, ignore.case = TRUE)))
  {
    Path_InFile <- readline("What is your VCF File? (please provide the full Path)")
  }
  if(any(grepl(VcfOrBaseCall, BaseCheck, ignore.case = TRUE)))
  {
    Path_InFile <- readline("What is your Base Call File? (please provide the full Path)")
  }
  
  Path_output_organism <- paste(Path_output, "/", organism, "/", sep = "")
  Path_output_triplet <- paste(Path_output_organism, "Triplet", sep = "")
  Path_output_upstream <- paste(Path_output_organism, "Upstream", sep = "")
  Path_output_downstream <- paste(Path_output_organism, "Downstream", sep = "")
  
  Path_correlate_triplet  <- paste(Path_correlate_repo, "triplet", sep = "/")
  Path_correlate_repo_down  <- paste(Path_correlate_repo, "Downstream", sep = "/")
  Path_correlate_repo_up  <- paste(Path_correlate_repo, "Upstream", sep = "/")
  
}

if(any(grepl(DirCheck, CustomCheck, ignore.case = TRUE)))
{
  Path_to_scripts <- readline("Where is your CDMAP Library Located?")
  Path_wd <- readline("What is your Working Directory?")
  Path_output <- readline("Where do you want to output your data?")
  Path_RefFile <- readline("What is your reference sequence? (please provide the full Path)")
  Path_GBFile <- readline("What is your Genbank file? (please provide the full Path)")
  Path_correlate_repo <- readline("Where would you like to store multi-organism output?")
  
  if(any(grepl(VcfOrBaseCall, VCFCheck, ignore.case = TRUE)))
  {
    Path_InFile <- readline("What is your VCF File? (please provide the full Path)")
  }
  if(any(grepl(VcfOrBaseCall, BaseCheck, ignore.case = TRUE)))
  {
    Path_InFile <- readline("What is your Base Call File? (please provide the full Path)")
  }

  
  Path_output_organism <- paste(Path_output, "/", organism, "/", sep = "")
  Path_output_triplet <- paste(Path_output_organism, "Triplet", sep = "")
  Path_output_upstream <- paste(Path_output_organism, "Upstream", sep = "")
  Path_output_downstream <- paste(Path_output_organism, "Downstream", sep = "")
  
  Path_correlate_triplet  <- paste(Path_correlate_repo, "triplet", sep = "/")
  Path_correlate_repo_down  <- paste(Path_correlate_repo, "Downstream", sep = "/")
  Path_correlate_repo_up  <- paste(Path_correlate_repo, "Upstream", sep = "/")
  
 }



generations <- readline("How many generations did you carry out your experiment? ")
malines <- readline("How many Mutation Accumulation Lines were run during your experiment? ")
malines <- as.numeric(malines)
generations <- as.numeric(generations)
param_flag <- readline("How would you like to scale your output data for the mutation counts, triplets, and mutation rates? 
         (0 for default, 1 for scaling to the average mean, or 2, for custom parameters): ")

#setwd(Path_to_scripts)
source("createDirectories.r") #CREATE ALL PREREQUISITE DIRECTORIES FOR OUTPUT

#========================================#



### REFERENCE FASTA FILE, SINGLE ORGANISM BASECALL, AND GENBANK, DATA MUNGING AND WRANGLING ###
#========================================#
# Read in the Fasta File and 
# Extract the sequence data from the fasta file and Partition the Base call CSV data
# NOTE: THIS CURRENTLY ONLY WORKS WITH *.FASTA extension files
RefFile <- getSequence(read.fasta(Path_RefFile)) #Path to fasta file
RefSeq_arr <<- unlist(RefFile) # converting the reference into Character array
RefSeq_arr <- toupper(RefSeq_arr) #conversion of reference sequence to all uppercase
len_refseq <<- length(RefSeq_arr) #value for length of the reference fasta sequence

if(any(grepl(VcfOrBaseCall, BaseCheck, ignore.case = TRUE)))
{
MutBaseCalls <- read.csv(file=Path_InFile, header=FALSE, sep=",")  #read in the basecall file
colnames(MutBaseCalls) <- c("position", "original", "mutant")
MutBaseCalls$position <- as.numeric(MutBaseCalls$position) #conversion of position columns to numerics for future operations
MutBaseCalls$original <- as.character(MutBaseCalls$original) #conversion of mutant and original columns to chars for 
MutBaseCalls$mutant <- as.character(MutBaseCalls$mutant)     #future operations
MutBaseCalls <- MutBaseCalls[order(MutBaseCalls$position),] #order the basecall file by the position column
}

if(any(grepl(VcfOrBaseCall, VCFCheck, ignore.case = TRUE)))
{
  source("VCFMutationParser.R")
  MutBaseCalls$position <- as.numeric(MutBaseCalls$position) #conversion of position columns to numerics for future operations
  MutBaseCalls$original <- as.character(MutBaseCalls$original) #conversion of mutant and original columns to chars for 
  MutBaseCalls$mutant <- as.character(MutBaseCalls$mutant)     #future operations
  MutBaseCalls <- MutBaseCalls[order(MutBaseCalls$position),] #order the basecall file by the position column
}
#========================================#



#==================================================
# Section 2: Determining Origin and terminus of the Query Sequence, and splitting
# them into left and right replichores
#=================================================

#TESTCODE TO PARSE THE GENBANK FILE FOR REPLICATION PROTEIN A (REPA) OR FOR DNA REPLICATION PROTEIN A (DNAA)
#################################

#for i in length(OrganismGB$FEATURES)
#{
#testobject <- paste("OrganismGB$FEATURES$'", "i", "'$gene", sep = "") #grabs the ith feature of the genbank file
#assign("geneloc", eval(parse(text = testobject))) #parses the gene column for text
#if(geneloc == 'dnaA' || 'repA')
#{
  #ori_pos <- matobj <- data.matrix(OrganismGB[["FEATURES"]][[i]])
  #featstart <- as.numeric(matobj[1,2]) #start of replication origin
  #featend <- as.numeric(matobj[1,3]) #end of replication origin region
#}

ori_ref <- oriloc(gbk = Path_GBFile) #this generates the oriloc data, this is a test variable
ori_pos <- which.max(ori_ref$skew) # find the position of the origin of replication
ori_bp <- ori_pos*1000   #base pair position of the origin of replication
ori_value <- ori_ref$skew[ori_pos] #value of the origination position
term_pos <- which.min(ori_ref$skew) #find the position of the terminus of replication
term_value <- ori_ref$skew[term_pos] #value of the terminii
term_bp <- term_pos*1000 #base pair position of the terminus of replication

#generates a dummy indices vector based on the length of the reference sequence
RefSeq_inds <<- c()
for(i in 1:len_refseq) 
{
  RefSeq_inds[i] <- i
}

#Use for Linear Chromosomes only (agro_radiobbacter_Chr1)
# ori_pos <- 496.696
# ori_bp <- 496696
# term_pos <-2499.261
# term_bp <- 2499261

#Use for Linear Chromosomes only (agro_radiobbacter_Chr2)
# ori_pos <- .147
# ori_bp <- 147
# term_pos <-1325.604
# term_bp <- 1325604


#Use for Linear Chromosomes only (agro_vitis_Chr1)
# ori_pos <- 336.131
# ori_bp <- 336131
# term_pos <- 2199.319
# term_bp <- 2199319

#Use for Linear Chromosomes only (agro_vitis_Chr2)
# ori_pos <- .918
# ori_bp <- 918
# term_pos <- 642.512
# term_bp <- 642512

#Use for Linear Chromosomes only
 # ori_pos <- 0
 # ori_bp <- 0
 # term_pos <- len_refseq/2
 # term_bp <- len_refseq/2

ori_bp <- as.numeric(ori_bp) #coercion of the ORI and TERM to numeric values for future operations
term_bp <- as.numeric(term_bp)

setwd(Path_to_scripts)
source("generateContextCounts.r")

#========================================================================
# CONVERSION INTO LEFT AND RIGHT REPLICHORES
#========================================================

# to convert them into left and right replichores, we have to
# handle it in a multicase if statement

query_Rcore_full <<- matrix( nrow =0, ncol =9)
colnames(query_Rcore_full) <- context_col_names


query_Lcore_full <<- matrix( nrow =0, ncol =9)
colnames(query_Lcore_full) <- context_col_names

#CASES
#Scenario 1: Origin < Terminus
i <-1
if(ori_bp < term_bp)
{
  setwd(Path_to_scripts)
  source("Core_OriTerm.r")
}


#Scenario 2: Origin > Terminus
i <- 1
if(ori_bp > term_bp)
{
  setwd(Path_to_scripts)
  source("Core_TermOri.r")
} 

setwd(Path_Basecall_output)
write.csv(query_Lcore_full, "Left_Context_Core_Full.csv")
write.csv(query_Rcore_full, "Right_Context_Core_Full.csv")


setwd(Path_Basecall_output_up)
write.csv(query_Lcore_full, "Left_Context_Core_Full.csv")
write.csv(query_Rcore_full, "Right_Context_Core_Full.csv")


setwd(Path_Basecall_output_down)
write.csv(query_Lcore_full, "Left_Context_Core_Full.csv")
write.csv(query_Rcore_full, "Right_Context_Core_Full.csv")


## Calculating GWTC, GC3C and GC4C upstream and downstream analysis ##
#======================#

starttime <- Sys.time()
GC_counter <- 0
CodingRegionSize <<- 0
setwd(Path_to_scripts)
source("GWTC.r")
setwd(Path_to_scripts)
source("Rev_GWTC.r")
setwd(Path_to_scripts)
source("GC4C.r")

endtime <- Sys.time()
runtime <- endtime - starttime

print(paste("The total runtime for GC4C calculation: ", runtime, sep = ""))
#=======================#


## Recording the Context dependent mutations and preparing for mutation rate calculation ##
## =================== ##

setwd(Path_to_scripts)
source("4fold_Mutation_Record.r")
#source("Record_Mutations.R")

## =================== ##

# Calculation of mutation rates using a comvination of data refactoring scripts ##
## =================== ##
setwd(Path_to_scripts)
source("3mer_refactor.r")


#GC3C Refactor
Flag_3mer <- TRUE


#Sung style Context output matrices
setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)
source("3mer_store.R") #storing the chromosome wide triplets for later

# Genomic Coding Region 4mers Rates and Visualization (GC4C)

#============= A upstream and downstream =============== #
Flag_3mer <- FALSE
setwd(Path_to_scripts)    
source("GC4C_refactor_Aup.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)

source("GC4C_refactor_Adown.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)

## =================== ##

#============= C upstream and downstream =============== #
Flag_3mer <- FALSE
setwd(Path_to_scripts)    
source("GC4C_refactor_Cup.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)

source("GC4C_refactor_Cdown.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)

## =================== ##

#============= G upstream and downstream =============== #
Flag_3mer <- FALSE
setwd(Path_to_scripts)    
source("GC4C_refactor_Gup.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)

source("GC4C_refactor_Gdown.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)

## =================== ##

#============= T upstream and downstream =============== #
Flag_3mer <- FALSE
setwd(Path_to_scripts)    
source("GC4C_refactor_Tup.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)

source("GC4C_refactor_Tdown.r") #coerce the objects into staging variables for individual analysis and visualization.

setwd(Path_to_scripts)    
source("Sung_Matrix_Left.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Right.R")
setwd(Path_to_scripts) 
source("Sung_Matrix_Chromosome.R")

setwd(Path_to_scripts)
source("MutationRate.r")
setwd(Path_to_scripts)
source("RevCompliment_MutationRate.r")
setwd(Path_to_scripts)
source("Context_Data_Visualization.R")
setwd(Path_to_scripts)
## =================== ##
CodingRegionPercent <- CodingRegionSize/len_refseq
CodonGCcontent <- GC_counter/CodingRegionSize
CodonGCcontent <- CodonGCcontent * 100
CodonATcontent <- 100 - CodonGCcontent
source("Config_Dump.R")


## =================== ##

print(paste("Context Dependent Mutation Pipeline Analysis of ", organism, " is complete!!", sep = ""))
beep(3) #Play the final fantasy victory music to remind user who alt tabbed the process is running is complete.