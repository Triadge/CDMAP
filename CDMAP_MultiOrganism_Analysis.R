#Heatmap Correlations in R using dynamic directories

packages <- c("seqinr", "BiocManager", "pracma", "beepr", "lattice", "tidyverse", "vcfR", "stringr", "reshape2", "Scales")
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
library("genbankr")

source("GC_Analysis_CreateDir.r")

setwd(path_to_scripts)
path_analyze <- path_to_chrome
source("GCcontent.r")
#path_analyze_rates <- path_to_chrome_rate
setwd(path_to_scripts)
name_addon <- "Chromosome"
source("Correlation_Script.r")
#source("ChiSquare_Script.r")

setwd(path_to_scripts)
path_analyze <- path_to_RevComp_chrome
source("GCcontent.r")
#path_analyze_rates <- path_to_chrome_rate
setwd(path_to_scripts)
name_addon <- "Counterclockwise_Chromosome"
source("Correlation_Script.r")

#----- insert mutation simulation, and proportion graphing scripts here
#==========#

#==Left Replichore Analysis==#
setwd(path_to_scripts)
path_analyze <- path_to_Lcore
#path_analyze_rates <- path_to_Lcore_rate
name_addon <- "Left_Replichore"
source("Correlation_Script.r")
##source("ChiSquare_Script.r")

setwd(path_to_scripts)
path_analyze <- path_to_RevComp_Lcore
##path_analyze_rates <- path_to_Lcore_rate
name_addon <- "Left_Replichore_ReverseCompliment"
source("Correlation_Script.r")
##source("ChiSquare_Script.r")

setwd(path_to_scripts)
path_analyze <- path_to_Lcore_Rep
##path_analyze_rates <- path_to_Lcore_rate
name_addon <- "Left_Replichore_RepSpecific"
source("Correlation_Script.r")

setwd(path_to_scripts)
path_analyze <- path_to_RevComp_Lcore_Rep
##path_analyze_rates <- path_to_Lcore_rate
name_addon <- "Left_Replichore_RepSpecific_ReverseCompliment"
source("Correlation_Script.r")

#----- insert mutation simulation, and proportion graphing scripts here
#==========#

#==Right Replichore Analysis==#
setwd(path_to_scripts)
path_analyze <- path_to_Rcore
##path_analyze_rates <- path_to_Rcore_rate
name_addon <- "Right Replichore"
source("Correlation_Script.r")
##source("ChiSquare_Script.r")
#----- insert mutation simulation, and proportion graphing scripts here

setwd(path_to_scripts)
path_analyze <- path_to_RevComp_Rcore
#path_analyze_rates <- path_to_Rcore_rate
name_addon <- "Right_Replichore_ReverseCompliment"
source("Correlation_Script.r")
#source("ChiSquare_Script.r")

setwd(path_to_scripts)
path_analyze <- path_to_Rcore_Rep
#path_analyze_rates <- path_to_Rcore_rate
name_addon <- "Right_Replichore_RepSpecific"
source("Correlation_Script.r")

setwd(path_to_scripts)
path_analyze <- path_to_RevComp_Rcore_Rep
#path_analyze_rates <- path_to_Rcore_rate
name_addon <- "Right_Replichore_RepSpecific_ReverseCompliment"
source("Correlation_Script.r")
#==========#





##Scratch space======================##

# check if a string contains a variable
# grep('Cat', testvar) where testvar <- "This is a Cat."
# note: grep will return 1 (true) or integer(0) for false

#list files in a directory#
# list.files(testdir), where testdir is a path to a given directory

#access the ith object in list.files as a string
#list.files(testdir)[2], where testdir is a directory path, and 2 is the 2nd element of the array

#round function for the heatmaps (for the titles)
# round(varname)

