#MutRate_ChiSq  Analysis.r#
#===========================#

#Objective: to develop expected mutation rates and run a chi square analysis of expected rates per organism

# step 1: mutation rate P(x)
# step 2: inverse (mutation rate P(x))
# step 3: inv(Column F)/sum(inv(column F))
# step 4:  sum of the observed*(step 3)

MutRateProbMatrix <- matrix(nrow = length(MutTriplet), ncol = 6)
MutRateColnames <- c("Codon", "Observed Counts", "Mutation Rate", "Inverse Mutation Rate", "Normalized Inv Mut Rate", "Expected Codon Usage")
colnames(MutRateProbMatrix) <- MutRateColnames
CodonTriplets <- gsub("\\[|\\]", "",MutTriplet)
#CodonTriplets <- sort(CodonTriplets)
MutRateProbMatrix[,1] <- CodonTriplets
MutRateProbMatrix[,2] <- observednum
MutRateProbMatrix[,3] <- observedFoldRate
#need to trim sites with no observed mutation rate

if(length(unique(MutRateProbMatrix) != 1))
{



#Detect zero mutation rate entries and delete/replace
for(s in length(MutRateProbMatrix[,1]):1) #reverse iterate through matrix to avoid subscript out of bound errors
{
  if(MutRateProbMatrix[s,3] == 0) #detect if there is a 0 mutation rate
  {
    MutRateProbMatrix <- MutRateProbMatrix[-s,] #remove row with 0 mutation rate
    #MutRateProbMatrix[s,3] <- 0.000000000001 #option to insert dummy mutation rate instead
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
MutRateProbMatrix[,6] <- SumObs*as.numeric(MutRateProbMatrix[,5])

Mut_File_Title <- "Mutation_Probability_Matrix.csv"
setwd(path_output_organism)
write.csv(MutRateProbMatrix, Mut_File_Title)
}
#need to remove row without altering subscript of the matrix!!!!
#for(s in length(MutRateProbMatrix[,1]):1)
#{
#  if(s %in% inds)
#  {
#    tempmatrix <- MutRateProbMatrix[-s,]
#  }
#}



#MutChiRawCodon <- chisq.test(as.numeric(MutRateProbMatrix[,2]), as.numeric(MutRateProbMatrix[,6]))
#MutchiPvalCodon <- c(organism, MutChiRawCodon$statistic, MutChiRawCodon$p.value)
#print(MutchiPvalCodon)
#MutChiSquarePvalCodon <- rbind(MutChiSquarePvalCodon, MutchiPvalCodon)

