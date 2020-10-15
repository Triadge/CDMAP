#Chi-Square Analysis Script
#purpose: Compute the Chi-Square test on each organism in the directory to determine If each sample is representative
#of the population


path_chisq_output <- paste(path_analyze, "/ChiSquare_Results", sep = "")
if(!(dir.exists(path_chisq_output)))
{
  dir.create(path_chisq_output)
}

#Array to hold organism name
names_array_rate <- c()
sung_names_array_rate <- c()
#ccw_names_array_rate <- c()
#Array to hold the name index of the matrix variable
matrix_index_array_rate <- c()
sung_matrix_index_array_rate <- c()
#ccw_matrix_index_array_rate <- c()

name_flag <- "Sung"
#CCW <- "Counterclockwise"
#CCW_output <- FALSE
directory <- list.files(path_analyze)
setwd(path_analyze)

#==== Parse the directory and seperate the outputs into respective groups for analysis ===#
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
  

    name <- unlist(strsplit(directory[i], split='.', fixed=TRUE))[1]
    name <- paste(name, "_rate", sep = "")
    names_array_rate <- c(names_array_rate, as.character(name))
    
    file_path <- as.character(paste(path_analyze, "/", directory[i], sep = "")) #grabs file path to read in matrix
    mat <- as.matrix(read.csv(file=file_path, header=TRUE, row.names = 1, sep=",")) #assigns matrix to temp variable
    assign(paste("Ubs_matrix_", i, sep = ""), mat)
    #assigns temp matrix to iterable matrix, which numerically corresponds to names array position
    matrix_index_array_rate <- c(matrix_index_array_rate, paste("Ubs_matrix_", i, sep = "")) 
    #creates index array used to iterated on matrices by numeric designation
  }

#===========#


#generate expected value matrix for each organism, and calculate the Chi Square statistic for each organism ==#

#=== generate probability P(x) for expected value E(xj) = P(x)*A[i,j]

probability_array <- c()
probability_matrix <- matrix(0L, nrow = 0, ncol = 4)
#ccw_probability_array <- c()
#ccw_probability_matrix <- matrix(0L, nrow = 0, ncol = 12)

cols <- c("T", "G", "C", "A")
#cols <- c("T→G", "T→C", "T→A",	"G→T",	"G→C",	"G→A",	"C→T",	"C→G",	"C→A",	"A→T",	"A→G",	"A→C")
rows <- c("T[X-->Y]T", "T[X-->Y]G","T[X-->Y]C","T[X-->Y]A","G[X-->Y]T","G[X-->Y]G","G[X-->Y]C","G[X-->Y]A",
          "C[X-->Y]T","C[X-->Y]G","C[X-->Y]C","C[X-->Y]A", "A[X-->Y]T","A[X-->Y]G", "A[X-->Y]C","A[X-->Y]A")
colnames(probability_matrix) <- cols
# colnames(sung_probability_matrix) <- cols
#colnames(ccw_probability_matrix) <- cols

i <- 1
j <- 1

#========== Calculate the expectation probability for each matrix type=======#
for(i in 1:length(names_array_rate))
{
  mat <- get(matrix_index_array_rate[i])
  for(j in 1:length(mat[1,]))
  {
    testvar <- sum(mat[,j])
    testsum <- testvar/length(mat[,j])
    
    #append(test_array, testsum)
    probability_array <- c(probability_array, testsum)
    j <- j+1
  }
  probability_matrix <- rbind(probability_matrix, probability_array)
  probability_array <-c()
  rownames(probability_matrix)[i] <- as.character(names_array_rate[i])
}

i <- 1
j <- 1




setwd(path_chisq_output)
write.csv(probability_matrix, "probability_Value_Matrix.csv")

#=====================================# 

#======= Calculate the Expected value matrix for each organism ======#
expectation_matrix <- matrix(0L, nrow = 16, ncol = 4)
# sung_expectation_matrix <- matrix(0L, nrow = 16, ncol = 4)
#ccw_expectation_matrix <- matrix(0L, nrow = 16, ncol = 12)

expectation_array_index <- list()
# sung_expectation_array_index <- list()
#ccw_expectation_array_index <- list()

i <- 1
j <- 1
k <- 1

prb_mat <- probability_matrix      #THIS LOOP IS STRUCTURED CORRECTLY! FIX CCW AND SUNG
for(i in 1:length(names_array_rate))
{
  mat <- get(matrix_index_array[i])
  prb_array <- probability_matrix[i,]

  for(j in 1:length(mat[1,]))
  {
    prb <- prb_array[j]
    prb <- round(as.numeric(prb), 20)
    for(k in 1:length(mat[,1]))
    {
      mut <- mat[k,j]
      exp_value <- mut*prb
      exp_value <- round(exp_value, 20)
      if(length(exp_value) == 0)
      {
        exp_value <- 0
      }
      expectation_matrix[k,j] <- as.numeric(exp_value)
      k <- k+1
    }
    j <- j+1
  }
  #assign(paste("exp_matrix", i, sep = ""), expectation_matrix)
  expectation_array_index[[i]] <- expectation_matrix 
  #paste("exp_matrix", i, sep = "")
  i <- i+1
}




#=======================================#

#====== Compute the Chi Square analysis of each organism and save it to an output file ====#
setwd(path_chisq_output)

i <- 1
Report_Name <- "ChiSquare_Summary_Report.txt"
sink(file = Report_Name, type = "output")
for(i in 1:length(names_array))
{
  obs_mat <- get(matrix_index_array[i])
  exp_mat <- expectation_array_index[[i]]
  cat(paste("Post Analysis summary of ", names_array[i], ":", sep = ""), append = TRUE)
  cat("\n=============================\n", append = TRUE)
  print(chisq.test(obs_mat, exp_mat, simulate.p.value = FALSE, correct = FALSE), append = TRUE)
  cat("\n=============================\n", append = TRUE)
  i <- i+1
}
sink()




print("Chi Square analysis of all organisms within database complete!")
#================================