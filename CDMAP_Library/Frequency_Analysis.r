
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


sung_cols <- c("T", "G", "C", "A")
rows <- c("T[X-->Y]T", "T[X-->Y]G","T[X-->Y]C","T[X-->Y]A","G[X-->Y]T","G[X-->Y]G","G[X-->Y]C","G[X-->Y]A",
          "C[X-->Y]T","C[X-->Y]G","C[X-->Y]C","C[X-->Y]A", "A[X-->Y]T","A[X-->Y]G", "A[X-->Y]C","A[X-->Y]A")
FrequencyMatrix <- matrix(0L, nrow =16, ncol =4)
rownames(FrequencyMatrix) <- rows
colnames(FrequencyMatrix) <- sung_cols

path_analyze <- "/Users/triadge/Desktop/PhD_Thesis/Data Analysis/frequency_Analysis/Chromosome"
path_output <- "/Users/triadge/Desktop/PhD_Thesis/Data Analysis/frequency_Analysis/"
#Array to hold organism name
names_array <- c()

#Array to hold the name index of the matrix variable
matrix_index_array <- c()


name_flag <- "Sung"
CCW <- "Counterclockwise"
CCW_output <- FALSE
directory <- list.files(path_analyze)
setwd(path_analyze)


for(i in 1:length(directory))
{
  if(directory[i] == "Correlation_Results")
  {
    next
  }
  if(directory[i] == "ChiSquare_Results")
  {
    next
  }
  
  name <- unlist(strsplit(directory[i], split='.', fixed=TRUE))[1] #Grabs the name of the Output organism file
  names_array <- c(names_array, as.character(name)) #Adds the name to the name indexing array
  
  file_path <- as.character(paste(path_analyze, "/", directory[i], sep = "")) #grabs file path to read in matrix
  mat <- as.matrix(read.csv(file=file_path, header=TRUE, row.names = 1, sep=",")) #assigns matrix to temp variable
  
  assign(paste("matrix", i, sep = ""), mat) # assigns matrix to the ith temp variable, corresponding to the ith entry in the
  # name indexing array to keep track of each matrix
  
  matrix_index_array <- c(matrix_index_array, paste("matrix", i, sep = "")) #stores the the matrix into an array where it's 
  #indice corresponds to the name of the temp matrix
  
}

for(i in 1:length(matrix_index_array))
{
  current <- get(matrix_index_array[i])
  FrequencyMatrix <- FrequencyMatrix+current
}

setwd(path_output)
write.csv(FrequencyMatrix, "Mutation_Frequency_Matrix.csv")


