### GC_MultiFold_Analysis.R###
#================================#
# Objective: 
# 1. (Complete) Conduct analysis of all N-fold sites using genomic GC content to generate the probability P(codon) of each nucleotide triplet
# 2. For each Triplet calculated, calculate the weighted probability of each triplet, W(P(codon)) = P(codon)/SUM(P(codon))
# 3. Calculate the expected # of codons Exp(codon) = W(P(codon))*SUM(OBS_codon)
# 4. Chi Square Test: chisq.test(OBS_codon, EXP_codon)
Codon_colnames <- c("Codon", "Probability", "Weighted Probability", "Expected Value")
CodonProbMatrix <- matrix(nrow = 64, ncol = 4)
colnames(CodonProbMatrix) <- Codon_colnames
CodonTriplets <- gsub("\\[|\\]", "",Triplet)
CodonTriplets <- sort(CodonTriplets)

p <-1
for(p in 1:length(CodonTriplets))
{
  Codon <- CodonTriplets[p]
  chararr <- toupper(unlist(strsplit(Codon, "")))
  
  if(chararr[1] == 'G' | chararr[1] == 'C')
  {
    P1 <- GC_content
  }else{
    P1 <- AT_content
  }
  if(chararr[2] == 'G' | chararr[2] == 'C')
  {
    P2 <- GC_content
  }else{
    P2 <- AT_content
  }
  if(chararr[3] == 'G' | chararr[3] == 'C')
  {
    P3 <- GC_content
  } else{
    P3 <- AT_content
  }
  
  CodonProbability <- P1*P2*P3
  
  CodonProbMatrix[p,1] <- Codon
  CodonProbMatrix[p,2] <- as.numeric(CodonProbability)  
}

SumCodonP <- sum(as.numeric(CodonProbMatrix[,2]))
CodonProbMatrix[,3] <- as.numeric(CodonProbMatrix[,2])/SumCodonP
Sum_observed <- sum(observed)
CodonProbMatrix[,4] <- as.numeric(CodonProbMatrix[,3])*CodonTotal

GC_File_Title <- "GC_Codon_Probability_Matrix.csv"
setwd(path_output_organism)
write.csv(CodonProbMatrix, GC_File_Title)


GCChiRawCodon <- chisq.test(observed, as.numeric(CodonProbMatrix[,4]))
GCchiPvalCodon <- c(organism, GCChiRawCodon$statistic, GCChiRawCodon$p.value)
print(GCchiPvalCodon)
GCChiSquarePvalCodon <- rbind(GCChiSquarePvalCodon, GCchiPvalCodon)

