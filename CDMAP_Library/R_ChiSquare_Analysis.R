
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

GCChiSquarePval <- matrix( nrow = 0, ncol = 3)
GCChiSquarePvalCodon <- matrix( nrow = 0, ncol = 3)
colnames(GCChiSquarePvalCodon) <- ChiNames

MutChiSquarePval <- matrix( nrow = 0, ncol = 3)
MutChiSquarePvalCodon <- matrix( nrow = 0, ncol = 3)
colnames(GCChiSquarePvalCodon) <- ChiNames


#Amino Acid sites and data structures by Fold Designation
#============================

#6-fold Sites
Argindices <- c("C[G]T", "C[G]G", "C[G]C", "C[G]A", "A[G]A", "A[G]G")
Leuindices <- c("C[T]T", "C[T]G", "C[T]C", "C[T]A", "T[T]G", "T[T]A")
Serindices <- c("T[C]T", "T[C]G", "T[C]C", "T[C]A", "A[G]T", "A[G]C")
ArgCodons <- c()
LeuCodons <- c()

#4-fold sites
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
PheCodons <- c()
TyrCodons <- c()
HisCodons <- c()
GlnCodons <- c()
AsnCodons <- c()
LysCodons <- c()
AspCodons <- c()
GluCodons <- c()
CysCodons <- c()

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

#Generate GC content for each Organism#
#===============================#
setwd(Path_to_scripts)
source("GC_content_Multifold.R")
#================================#

i <<- 1

#Begin Single Organism analysis#
#=================================#
for (j in 1:length(OrganismChrome))
  {
  
  #Insert GC content calculation script here#
  #pass OrganismChrome to the script#
  GC_content <- as.numeric(GC_output_matrix[j,2])
  AT_content <- 1-GC_content
  
  #===================#
  
    Atriplet <- c()
    Ttriplet <- c()
    Ctriplet <- c()
    Gtriplet <- c()
  
    organismSubset <- dataTable[i:(i+63),]
    organism <- OrganismChrome[j]
    Triplet <- organismSubset[,2]
    print(paste('running: ', organism, sep = ""))
    
    k <<- 1
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

  ObservedMaster <- organismSubset
  
  #Observed Mutation Rates
   observedRate <- c(Atriplet[,6], Ctriplet[,6], Gtriplet[,6], Ttriplet[,6])
   #observedRate <- as.matrix(observed)
    
  #Observed Mutations
  observedMut <- c(Atriplet[,5], Ctriplet[,5], Gtriplet[,5], Ttriplet[,5])
  #observedMut <- as.matrix(observed)
 
  #observed Codons 
  observedCodon <- c(Atriplet[,3], Ctriplet[,3], Gtriplet[,3], Ttriplet[,3])
  #observedCodon <- as.matrix(observed)
  
  ObsFlag <- "Codon"
  
  if(ObsFlag == "Rates")
  {
    observed <- observedRate
  }
  if(ObsFlag == "Counts")
  {
    observed <- observedMut
  }
  if(ObsFlag == "Codon")
  {
    observed <- observedCodon
  }
  

  CodonExpected <- organismSubset[,3]/sum(organismSubset[,3])
  CodonExpected <- as.matrix(CodonExpected)
  
  CodonTotal <- sum(as.numeric(organismSubset[,3]))
  GWTCExpected <- organismSubset[,4]/sum(organismSubset[,4])
  GWTCExpected <- as.matrix(GWTCExpected)
  
  #Insert GC specific analysis Script here#
  
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
   print(chiPvalCodon)
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
   source("GC_Multifold_Analysis.R")
   setwd(Path_to_scripts)
   source("MutRate_ChiSqAnalysis.R")
   print("running Fold Analysis")
   setwd(Path_to_scripts)
   #source("ChiSq_Foldsite.R")

  
  #context_output_downstream_matrix <- rbind(context_output_downstream_matrix, newrow_down) #appends new mutant to abbreviated text
  
  i = i+64

}

rownames(ChiSquarePvalCodon) <- NULL
rownames(MutChiSquarePvalCodon) <- NULL
rownames(GCChiSquarePvalCodon) <- NULL
rownames(ChiSquarePvalGWTC) <- NULL

setwd(path_output_SO)

write.csv(ChiSquarePvalCodon, "AllOrganisms_Pval_Codon.csv")
write.csv(MutChiSquarePvalCodon, "AllOrganisms_Mut_Pval_Codon.csv")
write.csv(GCChiSquarePvalCodon, "AllOrganisms_GC_Pval_Codon.csv")
write.csv(ChiSquarePvalGWTC, "AllOrganisms_Pval_GWTC.csv")





