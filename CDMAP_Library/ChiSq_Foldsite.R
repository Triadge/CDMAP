#chisquare in depth analysis for each organism


#steps
ArgCodons <- c()
LeuCodons <- c()
SerCodons <- c() 
ThrCodons <- c()
ValCodons <- c()
ProCodons <- c()
AlaCodons <- c()
GlyCodons <- c()
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
StopCodons <- c()
IleCodons <- c()
TrpCodons <- c()
StartCodons <- c()

OrgPvalCodon <- c()
OrgPvalGWTC <- c()

for(k in 1:length(organismSubset[,2]))
{
codon <- organismSubset[k,2]

#6-fold Sites
#===============================
if(any(codon == Argindices))
{
  ArgCodons <- rbind(ArgCodons, organismSubset[k,])
}
if(any(codon == Leuindices))
{
  LeuCodons <- rbind(LeuCodons, organismSubset[k,])
}
SixFold <- rbind(ArgCodons, LeuCodons, SerCodons)
#==================================

#4-fold sites
if(any(codon == Thrindices))
{
  ThrCodons <- rbind(ThrCodons, organismSubset[k,])
}
if(any(codon == Valindices))
{
  ValCodons <- rbind(ValCodons, organismSubset[k,])
}
if(any(codon == Proindices))
{
  ProCodons <- rbind(ProCodons, organismSubset[k,])
}
if(any(codon == Alaindices))
{
  AlaCodons <- rbind(AlaCodons, organismSubset[k,])
}
if(any(codon == Glyindices))
{
  GlyCodons <- rbind(GlyCodons, organismSubset[k,])
}

FourFold <- rbind(GlyCodons, AlaCodons, ProCodons, ValCodons, ThrCodons)
#============================

#2-fold sites
if(any(codon == Pheindices))
{
  PheCodons <- rbind(PheCodons, organismSubset[k,])
}
if(any(codon == Tyrindices))
{
  TyrCodons <- rbind(TyrCodons, organismSubset[k,])
}
if(any(codon == Hisindices))
{
  HisCodons <- rbind(HisCodons, organismSubset[k,])
}
if(any(codon == Glnindices))
{
  GlnCodons <- rbind(GlnCodons, organismSubset[k,])
}
if(any(codon == Asnindices))
{
  AsnCodons <- rbind(AsnCodons, organismSubset[k,])
}
if(any(codon == Lysindices))
{
  LysCodons <- rbind(LysCodons, organismSubset[k,])
}
if(any(codon == Aspindices))
{
  AspCodons <- rbind(AspCodons, organismSubset[k,])
}
if(any(codon == Gluindices))
{
  GluCodons <- rbind(GluCodons, organismSubset[k,])
}
if(any(codon == Cysindices))
{
  CysCodons <- rbind(CysCodons, organismSubset[k,])
}
TwoFold <- rbind( CysCodons, GluCodons, AspCodons, LysCodons, AsnCodons, GlnCodons, HisCodons, TyrCodons, PheCodons)
#====================================

#3-fold sites
#====================================
if(any(codon == Stopindices))
{
  StopCodons <- rbind(StopCodons, organismSubset[k,])
}
if(any(codon == Ileindices))
{
  IleCodons <- rbind(IleCodons, organismSubset[k,])
}
ThreeFold <- rbind( IleCodons, StopCodons)
#====================================

#1 fold
#=========================================
if(any(codon == Trpindices))
{
  TrpCodons <- rbind(TrpCodons, organismSubset[k,])
}
if(any(codon == Startindices))
{
  StartCodons <- rbind(StartCodons, organismSubset[k,])
}
OneFold <- rbind(StartCodons, TrpCodons)
#============================================

}  

#Compute ChisSquare analysis for each organism
#==========================================

#conversion for analysis of each fold site
SixFold <- data.frame(SixFold)
FourFold <- data.frame(FourFold)
ThreeFold <- data.frame(ThreeFold)
TwoFold <- data.frame(TwoFold)
OneFold <- data.frame(OneFold)

#Developer Notes: when trying to analyze the each individual N-fold sites as a group, some (or all) have a zero mutation rate
# this leads to a divide by zero error in the denominator in the Chi Squared analysis. need to think about how to analyze this.

#============= 6 Fold Analysis=============================#
  Atriplet <- c()
  Ttriplet <- c()
  Ctriplet <- c()
  Gtriplet <- c()
  Triplet <- SixFold[,2]
  fold <- "Sixfold"

  for(k in 1:length(Triplet))
  {
    if(grepl("[A]", Triplet[k], fixed = TRUE))
    {
      Atriplet <- rbind(Atriplet, SixFold[k,])
    }

    if(grepl("[C]", Triplet[k], fixed = TRUE))
    {
      Ctriplet <- rbind(Ctriplet, SixFold[k,])
    }

    if(grepl("[G]", Triplet[k], fixed = TRUE))
    {
      Gtriplet <- rbind(Gtriplet, SixFold[k,])
    }

    if(grepl("[T]", Triplet[k], fixed = TRUE))
    {
      Ttriplet <- rbind(Ttriplet, SixFold[k,])
    }
  }
  
  #observed Mutations
  observed <- c(Atriplet[,5], Ctriplet[,5], Gtriplet[,5], Ttriplet[,5])
  observed <- as.matrix(observed)

  #observed Codons
  #observed <- c(Atriplet[,3], Ctriplet[,3], Gtriplet[,3], Ttriplet[,3])
  #observed <- as.matrix(observed)
  
  
  CodonExpected <- SixFold[,3]/sum(SixFold[,3])
  CodonExpected <- as.matrix(CodonExpected)
  
  GWTCExpected <- SixFold[,4]/sum(SixFold[,4])
  GWTCExpected <- as.matrix(CodonExpected)
  
  setwd(Path_to_scripts)
  source("GC_Multifold_Analysis_Nfold.R")

if(sum(observed) != 0)
{
  
SixGCChiRawCodon <- chisq.test(observed, as.numeric(CodonProbMatrix[,4]))
SixChiRawCodon <- chisq.test(observed, p = CodonExpected) #, simulate.p.value = TRUE, B= 5000)
SixChiRawGWTC <- chisq.test(observed, p = GWTCExpected) #, simulate.p.value = TRUE, B = 5000)
SixchiPvalCodon <- c("Six Fold Sites", SixChiRawCodon$statistic, SixChiRawCodon$p.value, SixGCChiRawCodon$statistic, SixGCChiRawCodon$p.value)
SixchiPvalGWTC <- c( "Six Fold Sites", SixChiRawGWTC$statistic, SixChiRawGWTC$p.value)

OrgPvalCodon <- rbind(OrgPvalCodon, SixchiPvalCodon)
OrgPvalGWTC <- rbind(OrgPvalGWTC, SixchiPvalGWTC)
}
if(sum(observed) == 0)
{
  SixchiPvalCodon <- c("Six Fold Sites", "NA", "NA")
  SixchiPvalGWTC <- c("Six Fold Sites", "NA", "NA")
  OrgPvalCodon <- rbind(OrgPvalCodon, SixchiPvalCodon)
  OrgPvalGWTC <- rbind(OrgPvalGWTC, SixchiPvalGWTC)
  
}

print("Six Fold Complete")
#==========================================#


#============= 4 Fold Analysis=============================#
Atriplet <- c()
Ttriplet <- c()
Ctriplet <- c()
Gtriplet <- c()
Triplet <- FourFold[,2]
fold <- "Fourfold"

for(k in 1:length(Triplet))
{
  if(grepl("[A]", Triplet[k], fixed = TRUE))
  {
    Atriplet <- rbind(Atriplet, FourFold[k,])
  }
  
  if(grepl("[C]", Triplet[k], fixed = TRUE))
  {
    Ctriplet <- rbind(Ctriplet, FourFold[k,])
  }
  
  if(grepl("[G]", Triplet[k], fixed = TRUE))
  {
    Gtriplet <- rbind(Gtriplet, FourFold[k,])
  }
  
  if(grepl("[T]", Triplet[k], fixed = TRUE))
  {
    Ttriplet <- rbind(Ttriplet, FourFold[k,])
  }
}

#observed mutations
observed <- c(Atriplet[,5], Ctriplet[,5], Gtriplet[,5], Ttriplet[,5])
observed <- as.matrix(observed)

#observed Codons
#observed <- c(Atriplet[,3], Ctriplet[,3], Gtriplet[,3], Ttriplet[,3])
#observed <- as.matrix(observed)



CodonExpected <- FourFold[,3]/sum(FourFold[,3])
CodonExpected <- as.matrix(CodonExpected)

GWTCExpected <- FourFold[,4]/sum(FourFold[,4])
GWTCExpected <- as.matrix(CodonExpected)

setwd(Path_to_scripts)
source("GC_Multifold_Analysis_Nfold.R")

if(sum(observed) !=0)
{
FourGCChiRawCodon <- chisq.test(observed, as.numeric(CodonProbMatrix[,4]))
FourChiRawCodon <- chisq.test(observed, p = CodonExpected) #, simulate.p.value = TRUE, B= 5000)
FourChiRawGWTC <- chisq.test(observed, p = GWTCExpected) #, simulate.p.value = TRUE, B = 5000)

FourchiPvalCodon <- c("Four Fold Sites", FourChiRawCodon$statistic, FourChiRawCodon$p.value, FourGCChiRawCodon$statistic, FourGCChiRawCodon$p.value)
FourchiPvalGWTC <- c("Four Fold Sites", FourChiRawGWTC$statistic, FourChiRawGWTC$p.value)

OrgPvalCodon <- rbind(OrgPvalCodon, FourchiPvalCodon)
OrgPvalGWTC <- rbind(OrgPvalGWTC, FourchiPvalGWTC)
}
if(sum(observed) == 0)
{
  FourchiPvalCodon <- c("Four Fold Sites", "NA", "NA")
  FourchiPvalGWTC <- c("Four Fold Sites", "NA", "NA")
  OrgPvalCodon <- rbind(OrgPvalCodon, FourchiPvalCodon)
  OrgPvalGWTC <- rbind(OrgPvalGWTC, FourchiPvalGWTC)
  
}

print("Four Fold Complete")
#==========================================#

#============= 3 Fold Analysis=============================#
Atriplet <- c()
Ttriplet <- c()
Ctriplet <- c()
Gtriplet <- c()
Triplet <- ThreeFold[,2]


for(k in 1:length(Triplet))
{
  if(grepl("[A]", Triplet[k], fixed = TRUE))
  {
    Atriplet <- rbind(Atriplet, ThreeFold[k,])
  }
  
  if(grepl("[C]", Triplet[k], fixed = TRUE))
  {
    Ctriplet <- rbind(Ctriplet, ThreeFold[k,])
  }
  
  if(grepl("[G]", Triplet[k], fixed = TRUE))
  {
    Gtriplet <- rbind(Gtriplet, ThreeFold[k,])
  }
  
  if(grepl("[T]", Triplet[k], fixed = TRUE))
  {
    Ttriplet <- rbind(Ttriplet, ThreeFold[k,])
  }
}

#observed mutations
observed <- c(Atriplet[,5], Ctriplet[,5], Gtriplet[,5], Ttriplet[,5])
observed <- as.matrix(observed)

#observed Codons
#observed <- c(Atriplet[,3], Ctriplet[,3], Gtriplet[,3], Ttriplet[,3])
#observed <- as.matrix(observed)

CodonExpected <- ThreeFold[,3]/sum(ThreeFold[,3])
CodonExpected <- as.matrix(CodonExpected)

GWTCExpected <- ThreeFold[,4]/sum(ThreeFold[,4])
GWTCExpected <- as.matrix(CodonExpected)

setwd(Path_to_scripts)
source("GC_Multifold_Analysis_Nfold.R")

if(sum(observed) != 0)
{

ThreeGCChiRawCodon <- chisq.test(observed, as.numeric(CodonProbMatrix[,4]))
ThreeChiRawCodon <- chisq.test(observed, p = CodonExpected) #, simulate.p.value = TRUE, B= 5000)
ThreeChiRawGWTC <- chisq.test(observed, p = GWTCExpected) #, simulate.p.value = TRUE, B = 5000)

ThreechiPvalCodon <- c("Three Fold Sites", ThreeChiRawCodon$statistic, ThreeChiRawCodon$p.value, ThreeGCChiRawCodon$statistic, ThreeGCChiRawCodon$p.value)
ThreechiPvalGWTC <- c("Three Fold Sites", ThreeChiRawGWTC$statistic, ThreeChiRawGWTC$p.value)

OrgPvalCodon <- rbind(OrgPvalCodon, ThreechiPvalCodon)
OrgPvalGWTC <- rbind(OrgPvalGWTC, ThreechiPvalGWTC)
print("Three Fold Complete")
}

if(sum(observed) == 0)
{
  ThreechiPvalCodon <- c("Three Fold Sites", "NA", "NA")
  ThreechiPvalGWTC <- c("Three Fold Sites", "NA", "NA")
  OrgPvalCodon <- rbind(OrgPvalCodon, ThreechiPvalCodon)
  OrgPvalGWTC <- rbind(OrgPvalGWTC, ThreechiPvalGWTC)
}

print("Three Fold Complete")
#==========================================#

#============= 2 Fold Analysis=============================#
Atriplet <- c()
Ttriplet <- c()
Ctriplet <- c()
Gtriplet <- c()
Triplet <- TwoFold[,2]


for(k in 1:length(Triplet))
{
  if(grepl("[A]", Triplet[k], fixed = TRUE))
  {
    Atriplet <- rbind(Atriplet, TwoFold[k,])
  }

  if(grepl("[C]", Triplet[k], fixed = TRUE))
  {
    Ctriplet <- rbind(Ctriplet, TwoFold[k,])
  }

  if(grepl("[G]", Triplet[k], fixed = TRUE))
  {
    Gtriplet <- rbind(Gtriplet, TwoFold[k,])
  }

  if(grepl("[T]", Triplet[k], fixed = TRUE))
  {
    Ttriplet <- rbind(Ttriplet, TwoFold[k,])
  }
}

#observed mutations
observed <- c(Atriplet[,5], Ctriplet[,5], Gtriplet[,5], Ttriplet[,5])
observed <- as.matrix(observed)

#observed Codons
#observed <- c(Atriplet[,3], Ctriplet[,3], Gtriplet[,3], Ttriplet[,3])
#observed <- as.matrix(observed)

CodonExpected <- TwoFold[,3]/sum(TwoFold[,3])
CodonExpected <- as.matrix(CodonExpected)

GWTCExpected <- TwoFold[,4]/sum(TwoFold[,4])
GWTCExpected <- as.matrix(CodonExpected)

setwd(Path_to_scripts)
source("GC_Multifold_Analysis_Nfold.R")

if(sum(observed) != 0 )
{

TwoGCChiRawCodon <- chisq.test(observed, as.numeric(CodonProbMatrix[,4]))
TwoChiRawCodon <- chisq.test(observed, p = CodonExpected) #, simulate.p.value = TRUE, B= 5000)
TwoChiRawGWTC <- chisq.test(observed, p = GWTCExpected) #, simulate.p.value = TRUE, B = 5000)

TwochiPvalCodon <- c("Two Fold Sites", TwoChiRawCodon$statistic, TwoChiRawCodon$p.value, TwoGCChiRawCodon$statistic, TwoGCChiRawCodon$p.value)
TwochiPvalGWTC <- c("Two Fold Sites", TwoChiRawGWTC$statistic, TwoChiRawGWTC$p.value)

OrgPvalCodon <- rbind(OrgPvalCodon, TwochiPvalCodon)
OrgPvalGWTC <- rbind(OrgPvalGWTC, TwochiPvalGWTC)
}
if(sum(observed) == 0)
{
  TwochiPvalCodon <- c("Two Fold Sites", "NA", "NA")
  TwochiPvalGWTC <- c("Two Fold Sites", "NA", "NA")
  OrgPvalCodon <- rbind(OrgPvalCodon, TwochiPvalCodon)
  OrgPvalGWTC <- rbind(OrgPvalGWTC, TwochiPvalGWTC)
  
}

print("Two Fold Complete")
# #==========================================#
# 
# # #============= 1 Fold Analysis=============================#
# Atriplet <- c()
# Ttriplet <- c()
# Ctriplet <- c()
# Gtriplet <- c()
# Triplet <- OneFold[,2]
# 
# for(k in 1:length(Triplet))
# {
#   if(grepl("[A]", Triplet[k], fixed = TRUE))
#   {
#     Atriplet <- rbind(Atriplet, OneFold[k,])
#   }
# 
#   if(grepl("[C]", Triplet[k], fixed = TRUE))
#   {
#     Ctriplet <- rbind(Ctriplet, OneFold[k,])
#   }
# 
#   if(grepl("[G]", Triplet[k], fixed = TRUE))
#   {
#     Gtriplet <- rbind(Gtriplet, OneFold[k,])
#   }
# 
#   if(grepl("[T]", Triplet[k], fixed = TRUE))
#   {
#     Ttriplet <- rbind(Ttriplet, OneFold[k,])
#   }
# }
# 
# observed <- c(Atriplet[,5], Ctriplet[,5], Gtriplet[,5], Ttriplet[,5])
# observed <- as.matrix(observed)
# 
# 
# 
# CodonExpected <- OneFold[,3]/sum(OneFold[,3])
# CodonExpected <- as.matrix(CodonExpected)
# 
# GWTCExpected <- OneFold[,4]/sum(OneFold[,4])
# GWTCExpected <- as.matrix(CodonExpected)
# 
# OneChiRawCodon <- chisq.test(observed, CodonExpected)
# OneChiRawGWTC <- chisq.test(observed, GWTCExpected)
# OnechiPvalCodon <- c("One Fold Sites", ChiRawCodon$statistic, ChiRawCodon$p.value)
# OnechiPvalGWTC <- c("One Fold Sites", ChiRawGWTC$statistic, ChiRawGWTC$p.value)
# #==========================================#


rownames(OrgPvalCodon) <- NULL
rownames(OrgPvalGWTC) <- NULL
OrgColNames <- c("Fold Site", "Chi Stat", "P-Value", "GC Chi Stat", "GC P-Value")
colnames(OrgPvalCodon) <- OrgColNames
setwd(path_output_organism)

write.csv(OrgPvalCodon, "ChiSquare_Codon_ByFoldSite.csv")
write.csv(OrgPvalGWTC, "ChiSquare_GWTC_ByFoldSite.csv")
#generate the expected value and calculate the Chisquare for each individual subset of each

#generate the ChiSquare of a given organism for all sites

