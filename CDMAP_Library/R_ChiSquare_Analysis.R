
packages <- c("seqinr", "BiocManager", "pracma", "beepr", "lattice", "tidyverse", "vcfR", "stringr", "reshape2", "scales")
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

path_output <- "/Users/triadge/Desktop/PhD_Thesis/Data Analysis/MultiOrganismAnalysis/ChiSquare_Fisher_Dumps"

dataTable <- read.csv("/Users/triadge/Desktop/PhD_Thesis/Data Analysis/MultiOrganismAnalysis/MultiOrganism_Rupload.csv")

#dataTable <- sort(dataTable[,1])
#OrganismChrome <- sort(unique(dataTable[,1]))
HeaderNames <- names(dataTable)

FishNames <- c("organism", "Fisher75", "Fisher50", "Fisher25", "Fisher10")
ChiNames <- c("organism", "Chi75", "Chi50", "Chi25", "Chi10")

ChiSquarePval <- matrix( nrow = 0, ncol = 5)
FisherExactPval <- matrix( nrow = 0, ncol = 5)
colnames(ChiSquarePval) <- ChiNames
colnames(FisherExactPval) <- FishNames

ChiSquareStat <- matrix( nrow = 0, ncol = 5)
FisherExactStat <- matrix( nrow = 0, ncol = 5)
colnames(ChiSquareStat) <- ChiNames
colnames(FisherExactStat) <- FishNames


i <- 1
for (i in 1:nrow(dataTable))
{
  
  organismSubset <- dataTable[i:(i+63),]
  organism <- organismSubset[i,1]
  observed <- organismSubset[,10]
  expected75 <- organismSubset[,3]*organismSubset[,6]
  expected50 <- organismSubset[,3]*organismSubset[,7]
  expected25 <- organismSubset[,3]*organismSubset[,8]
  expected10 <- organismSubset[,3]*organismSubset[,9]

  Chi75 <- chisq.test(observed, expected75)
  Chi50 <- chisq.test(observed, expected50)
  Chi25 <- chisq.test(observed, expected25)
  Chi10 <- chisq.test(observed, expected10)
  
  Fish75 <- fisher.test(observed, expected75, simulate.p.value = TRUE)
  Fish50 <- fisher.test(observed, expected50, simulate.p.value = TRUE)
  Fish25 <- fisher.test(observed, expected25, simulate.p.value = TRUE)
  Fish10 <- fisher.test(observed, expected10, simulate.p.value = TRUE)
  
  chiPval <- c(organism, Chi75$p.value, Chi50$p.value, Chi25$p.value, Chi10$p.value)
  fishPval <- c(organism, Fish75$p.value, Fish50$p.value, Fish25$p.value, Fish10$p.value)
  
  chiStat <- c(organism, Chi75$statistic, Chi50$statistic, Chi25$statistic, Chi10$statistic)
  fishStat <- c(organism, Fish75$statistic, Fish50$statistic, Fish25$statistic, Fish10$statistic)
  
  ChiSquarePval <- rbind(ChiSquarePval, chiPval)
  FisherExactPval <- rbind(FisherExactPval, fishPval)
  
  ChiSquareStat <- rbind(ChiSquareStat, chiStat)
  FisherExactStat <- rbind(FisherExactStat, fishStat)
  
  #context_output_downstream_matrix <- rbind(context_output_downstream_matrix, newrow_down) #appends new mutant to abbreviated text
  
  
  i <- i+64
}

cd(path_output)

write.csv(ChiSquarePval, "ChiSquare_AllOrganisms_Pval.csv")
write.csv(FisherExactPval, "FisherTest_AllOrganisms_Pval.csv")

write.csv(ChiSquareStat, "ChiSquare_AllOrganisms_Stat.csv")
write.csv(FisherExactStat, "FisherTest_AllOrganisms_Stat.csv")


