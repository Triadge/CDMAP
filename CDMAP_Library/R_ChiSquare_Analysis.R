
# packages <- c("seqinr", "BiocManager", "pracma", "beepr", "lattice", "tidyverse", "vcfR", "stringr", "reshape2", "scales")
# package.check <- lapply(
#   packages,
#   FUN = function(x) {
#     if (!require(x, character.only = TRUE)) {
#       install.packages(x, dependencies = TRUE)
#       library(x, character.only = TRUE)
#     }
#   }
# )
# #BiocManager::install("genbankr")
# library("genbankr")

path_output <- "/Users/triadge/Desktop/PhD_Thesis/Data Analysis/MultiOrganismAnalysis/ChiSquare_Fisher_Dumps"
path_output_SO <- paste(path_output, "/Single_Organism_Output", sep = "")
Path_to_scripts <- "/Users/triadge/Desktop/CDMAP/CDMAP_Library"

if(!(dir.exists(path_output_SO)))
{
  dir.create(path_output_SO)
}

dataTable <- read.csv("/Users/triadge/Desktop/PhD_Thesis/Data Analysis/MultiOrganismAnalysis/ChiSquare_Fisher_Dumps/input_files/MultiOrganism_Rupload.csv")
#dataTable <- read.csv("/Users/triadge/Desktop/PhD_Thesis/Data Analysis/MultiOrganismAnalysis/Agro_Chr1_Test.csv")

#generations <- readline("How many generations did you carry out your experiment? ")
#generations <- as.numeric(generations)
#malines <- readline("How many Mutation Accumulation Lines were run during your experiment? ")
#malines <- as.numeric(malines)

#dataTable <- sort(dataTable[,1])
OrganismChrome <- sort(unique(dataTable[,1]))
HeaderNames <- names(dataTable)


ChiNames <- c("Organism", "RawChiStat", "RawChiPval")



ChiSquarePval <- matrix( nrow = 0, ncol = 3)
ChiSquarePvalCodon <- matrix( nrow = 0, ncol = 3)
ChiSquarePvalGWTC <- matrix( nrow = 0, ncol = 3)
OrgPvalCodon <- matrix( nrow = 0, ncol = 3)
OrgPvalGWTC <- matrix( nrow = 0, ncol = 3)

colnames(ChiSquarePvalCodon) <- ChiNames
colnames(ChiSquarePvalGWTC) <- ChiNames
colnames(OrgPvalCodon) <- ChiNames
colnames(OrgPvalGWTC) <- ChiNames



#Amino Acid sites and data structures by Fold Designation
#============================

#6-fold Sites
Argindices <- c("C[G]T", "C[G]G", "C[G]C", "C[G]A", "A[G]A", "A[G]G")
Leuindices <- c("C[T]T", "C[T]G", "C[T]C", "C[T]A", "T[T]G", "T[T]A")
ArgCodons <- c()
LeuCodons <- c()

#4-fold sites
Serindices <- c("T[C]T", "T[C]G", "T[C]C", "T[C]A") 
Thrindices <- c("A[C]T", "A[C]G", "A[C]C", "A[C]A")
Valindices <- c("G[T]T", "G[T]G", "G[T]C", "G[T]A")
Proindices <- c("C[C]T", "C[C]G", "C[C]C", "C[C]A")
Alaindices <- c("G[C]T", "G[C]G", "G[C]C", "G[C]A")
Glyindices <- c("G[G]T", "G[G]G", "G[G]C", "G[G]A")
SerCodons <- c() 
ThrCodons <- c()
ValCodons <- c()
ProCodons <- c()
AlaCodons <- c()
GlyCodons <- c()

#2-fold sites
Pheindices <- c("T[T]T", "T[T]C")
Tyrindices <- c("T[A]C", "T[A]T")
Hisindices <- c("C[A]T", "C[A]C")
Glnindices <- c("C[A]A", "C[A]G")
Asnindices <- c("A[A]T", "A[A]C")
Lysindices <- c("A[A]A", "A[A]G")
Aspindices <- c("G[A]T", "G[A]C")
Gluindices <- c("G[A]A", "G[A]G")
Cysindices <- c("T[G]T", "T[G]C")
Serindices <- c("A[G]T", "A[G]C")
PheCodons <- c()
TyrCodons <- c()
HisCodons <- c()
GlnCodons <- c()
AsnCodons <- c()
LysCodons <- c()
AspCodons <- c()
GluCodons <- c()
CysCodons <- c()
SerCodons <- c()

#3-fold sites
Stopindices <- c("T[A]G", "T[A]A", "T[G]A")
Ileindices <- c( "A[T]T", "A[T]C", "A[T]A")
StopCodons <- c()
IleCodons <- c()

#1 fold
Trpindices <- c("T[G]G")
Startindices <- c("A[T]G")
TrpCodons <- c()
StartCodons <- c()

#======================================

Atriplet <- c()
Ttriplet <- c()
Ctriplet <- c()
Gtriplet <- c()


i <<- 1
j <<- 1

for (j in 1:length(OrganismChrome))
  {
    Atriplet <- c()
    Ttriplet <- c()
    Ctriplet <- c()
    Gtriplet <- c()
  
    organismSubset <- dataTable[i:(i+63),]
    organism <- OrganismChrome[j]
    Triplet <- organismSubset[,2]
    print(paste('running: ', organism, sep = ""))
  
    for(k in 1:length(Triplet))
      {
        if(grepl("[A]", Triplet[k], fixed = TRUE))
          {
            Atriplet <- rbind(Atriplet, organismSubset[k,])
          }
      
        if(grepl("[C]", Triplet[k], fixed = TRUE))
          {
            Ctriplet <- rbind(Ctriplet, organismSubset[k,])
          }
      
        if(grepl("[G]", Triplet[k], fixed = TRUE))
          {
          Gtriplet <- rbind(Gtriplet, organismSubset[k,])
          }
      
        if(grepl("[T]", Triplet[k], fixed = TRUE))
          {
          Ttriplet <- rbind(Ttriplet, organismSubset[k,])
          }
      }
  #Triplet
  observed <- c(Atriplet[,5], Ctriplet[,5], Gtriplet[,5], Ttriplet[,5])
  observed <- as.matrix(observed)
  
  #AsumCodon <- sum(Atriplet[,3])
  #CsumCodon <- sum(Ctriplet[,3])
  #GsumCodon <- sum(Gtriplet[,3])
  #TsumCodon <- sum(Ttriplet[,3])
  
  #CodonExpectedA <- Atriplet[,3]/AsumCodon
  #CodonExpectedC <- Ctriplet[,3]/CsumCodon
  #CodonExpectedG <- Gtriplet[,3]/GsumCodon
  #CodonExpectedT <- Ttriplet[,3]/TsumCodon
  #CodonExpected <- c(CodonExpectedA, CodonExpectedC, CodonExpectedG, CodonExpectedT)
  CodonExpected <- organismSubset[,3]/sum(organismSubset[,3])
  CodonExpected < as.matrix(CodonExpected)
  
  #AsumGWTC <- sum(Atriplet[,4])
  #CsumGWTC <- sum(Ctriplet[,4])
  #GsumGWTC <- sum(Gtriplet[,4])
  #TsumGWTC <- sum(Ttriplet[,4])
  
  #GWTCExpectedA <- Atriplet[,4]/AsumGWTC
  #GWTCExpectedC <- Ctriplet[,4]/CsumGWTC
  #GWTCExpectedG <- Gtriplet[,4]/GsumGWTC
  #GWTCExpectedT <- Ttriplet[,4]/TsumGWTC
  #GWTCExpected <- c(GWTCExpectedA, GWTCExpectedC, GWTCExpectedG, GWTCExpectedT)
  GWTCExpected <- organismSubset[,4]/sum(organismSubset[,4])
  GWTCExpected <- as.matrix(GWTCExpected)
  
  path_output_organism <- paste(path_output_SO, organism, sep = "/")
  
  if(!(dir.exists(path_output_organism)))
    {
      dir.create(path_output_organism)
    }

  
 # print(paste("Organism", organism,"iterator", j, sep = " "))
 #  print(paste("Triplet iterator", i, sep = " "))
  
   ChiRawCodon <- chisq.test(observed, p = CodonExpected)
   #print(ChiRawCodon$p.value)
   ChiRawGWTC <- chisq.test(observed, p = GWTCExpected)
   #print(ChiRawGWTC$p.value)

     
   chiPvalCodon <- c(organism, ChiRawCodon$statistic, ChiRawCodon$p.value)
   chiPvalGWTC <- c(organism, ChiRawGWTC$statistic, ChiRawGWTC$p.value)

  #  
   ChiSquarePvalCodon <- rbind(ChiSquarePvalCodon, chiPvalCodon)
   ChiSquarePvalGWTC <- rbind(ChiSquarePvalGWTC, chiPvalGWTC)

   setwd(path_output_organism)
   codonout <- cbind(Triplet, CodonExpected)
   gwtcout <- cbind(Triplet, GWTCExpected)
   write.csv(codonout, "ExpectedCodons.csv")
   write.csv(gwtcout, "ExpectedGWTC.csv")
   

   
   setwd(Path_to_scripts)
   print("running Fold Analysis")
   source("ChiSq_Foldsite.R")

  
  #context_output_downstream_matrix <- rbind(context_output_downstream_matrix, newrow_down) #appends new mutant to abbreviated text
  
  i = i+64

}

rownames(ChiSquarePvalCodon) <- NULL
rownames(ChiSquarePvalGWTC) <- NULL

setwd(path_output_SO)

write.csv(ChiSquarePvalCodon, "AllOrganisms_Pval_Codon.csv")
write.csv(ChiSquarePvalGWTC, "AllOrganisms_Pval_GWTC.csv")





