#MutRate_ChiSq  Analysis.r#
#===========================#

#Objective: to develop expected mutation rates and run a chi square analysis of expected rates per organism

# step 1: mutation rate P(x)
# step 2: inverse (mutation rate P(x))
# step 3: inv(Column F)/sum(inv(column F))
# step 4:  sum of the observed*(step 3)

MutRateProbMatrix <- matrix(nrow = length(Triplet), ncol = 6)
MutRateColnames <- c("codon", "Observed Counts", "Mutation Rate", "Inverse Mutation Rate", "Normalized Inv Mut Rate", "Expected Codon Usage")
colnames(MutRateProbMatrix) <- MutRateColnames
CodonTriplets <- gsub("\\[|\\]", "",Triplet)
#CodonTriplets <- sort(CodonTriplets)
MutRateProbMatrix[,1] <- CodonTriplets
MutRateProbMatrix[,2] <- as.numeric(observedMut)
MutRateProbMatrix[,3] <- as.numeric(observedRate)
#need to trim sites with no observed mutation rate
inds <- c()
for(s in 1:length(MutRateProbMatrix[,1]))
{
  if(MutRateProbMatrix[s,3] == 0)
  {
    #inds <- c(inds, s)
    #replace a 0 mutation rate site with .01 mutations/(site)(generation)(lines)
    MutRateProbMatrix[s,3] <- 0.000000000001
  }
}

for(k in 1:length(MutRateProbMatrix[,1]))
{

    #inds <- c(inds, s)
    #replace a 0 mutation rate site with .01 mutations/(site)(generation)(lines)
    MutRateProbMatrix[k,4] <- 1/as.numeric(MutRateProbMatrix[k,3])
}

SumInvMutRate <- sum(as.numeric(MutRateProbMatrix[,4]))
SumObs <- sum(as.numeric(MutRateProbMatrix[,2]))
MutRateProbMatrix[,5] <- as.numeric(MutRateProbMatrix[,4])/ SumInvMutRate
MutRateProbMatrix[,6] <- sum(observedCodon)*as.numeric(MutRateProbMatrix[,5])

#need to remove row without altering subscript of the matrix!!!!
#for(s in length(MutRateProbMatrix[,1]):1)
#{
#  if(s %in% inds)
#  {
#    tempmatrix <- MutRateProbMatrix[-s,]
#  }
#}

Mut_File_Title <- "Mutation_Probability_Matrix.csv"
setwd(path_output_organism)
write.csv(MutRateProbMatrix, Mut_File_Title)

chiframe <- data.frame(as.numeric(observedCodon), as.numeric(MutRateProbMatrix[,6]))
MutChiRawCodon <- chisq.test(chiframe)
MutchiPvalCodon <- c(organism, MutChiRawCodon$statistic, MutChiRawCodon$p.value)
print(MutchiPvalCodon)
MutChiSquarePvalCodon <- rbind(MutChiSquarePvalCodon, MutchiPvalCodon)

