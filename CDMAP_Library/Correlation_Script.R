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


#====================Output Matrix Creation========================#
#p value matrix
cor_pval_matrix <- matrix(0L, nrow = length(names_array), ncol = length(names_array))
colnames(cor_pval_matrix) <- names_array
rownames(cor_pval_matrix) <- names_array

#correlation test statistic matrix
cor_stat_matrix <- matrix(0L, nrow = length(names_array), ncol = length(names_array))
colnames(cor_stat_matrix) <- names_array
rownames(cor_stat_matrix) <- names_array

#Correlation estimate matrix
cor_est_matrix <- matrix(0L, nrow = length(names_array), ncol = length(names_array))
colnames(cor_est_matrix) <- names_array
rownames(cor_est_matrix) <- names_array



#======== Generate the 1:N correlation output matrices of all organisms in the output csv list ====== #

for(i in 1:length(names_array))
{
  if(length(matrix_index_array) != 0)
  {
  object1 <- get(matrix_index_array[i])
  
  for(j in 1:length(names_array))
  {
    if(j < length(names_array))
    {
      object2 <- get(matrix_index_array[j+1])
      #print(paste("corr object 1: ", as.character(names_array[i]), "corr object 2: ", as.character(names_array[j+1]) ))
      corr_object <- cor.test(object1, object2)
      cor_pval_matrix[i, j+1] <- corr_object$p.value
      cor_pval_matrix[j+1, i] <- corr_object$p.value
      cor_stat_matrix[i, j+1] <- corr_object$statistic
      cor_stat_matrix[j+1,i] <- corr_object$statistic
      cor_est_matrix[i, j+1] <- corr_object$estimate
      cor_est_matrix[j+1,i] <- corr_object$estimate
    }
    else
    {
      object2 <- get(matrix_index_array[1])
      #print(paste("corr object 1: ", as.character(names_array[i]), "corr object 2: ", as.character(names_array[1])))
      corr_object <- cor.test(object1, object2)
      cor_pval_matrix[i, 1] <- corr_object$p.value
      cor_pval_matrix[1, i] <- corr_object$p.value
      cor_stat_matrix[i, 1] <- corr_object$statistic
      cor_stat_matrix[1,i] <- corr_object$statistic
      cor_est_matrix[i, 1] <- corr_object$estimate
      cor_est_matrix[1,i] <- corr_object$estimate
    }
  }
  }
  
}

#Generate GC Ccontent based correlation matrix

GC_corr_matrix <- matrix(0L, nrow =dim(GC_output_matrix)[1], ncol = dim(GC_output_matrix)[1])
GC_arr <- as.character(GC_sort_matrix$name)
colnames(GC_corr_matrix) <- GC_arr
rownames(GC_corr_matrix) <- GC_arr


for(i in 1:length(GC_arr))
{
  if(length(cor_est_matrix) != 0)
  {  
    GC_row_flag <- GC_arr[i]
    rownum <- grep(GC_row_flag, rownames(cor_est_matrix), ignore.case=TRUE)
    for(j in 1:length(GC_arr))
    {
      GC_col_flag <- GC_arr[j]
      colnum <- grep(GC_col_flag, colnames(cor_est_matrix), ignore.case=TRUE)
      GC_corr_matrix[i,j] <- cor_est_matrix[rownum, colnum]
    }
  }
}



#============================================#

#=====Create output directory for correlation spreadsheets and output graphics ====#
path_correlate_output <- paste(path_analyze, "/Correlation_Results", sep = "")
if(!(dir.exists(path_correlate_output)))
{
  dir.create(path_correlate_output)
}

#output the text spreadsheets before generating the graphic
setwd(path_correlate_output)
write.csv(cor_est_matrix, "Correlation.csv")
write.csv(GC_corr_matrix, "GC_Correlation.csv")
write.csv(cor_stat_matrix, "Correlation_TestStatistic.csv")
write.csv(cor_pval_matrix, "Correlation_Pvalue.csv")


if(length(cor_est_matrix) != 0)
{
setwd(LibDir)
matrix_name <- "1:N Correlation Graph "
output_data_matrix <- cor_est_matrix
matrix_name_graph <- paste(matrix_name, name_addon, sep = "-")
source("Correlation_Visualizer.r")
}


if(length(GC_corr_matrix) != 0)
{
  setwd(LibDir)
  matrix_name <- paste("GC Correlation Graph ", name_addon, sep = "-")
  output_data_matrix <- GC_corr_matrix
  matrix_name_graph <- matrix_name
  source("Correlation_Visualizer.r")

}
setwd(LibDir)
#=====================================================

