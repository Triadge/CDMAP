#Scratchy boi

#setwd(Path_to_scripts)
#source("4fold_Mutation_Record.r")

Base_Test <- matrix(0L, nrow =16, ncol =12)
cols <- c("T→G", "T→C", "T→A",	"G→T",	"G→C",	"G→A",	"C→T",	"C→G",	"C→A",	"A→T",	"A→G",	"A→C")
rows <- c("T[X-->Y]T", "T[X-->Y]G","T[X-->Y]C","T[X-->Y]A","G[X-->Y]T","G[X-->Y]G","G[X-->Y]C","G[X-->Y]A",
          "C[X-->Y]T","C[X-->Y]G","C[X-->Y]C","C[X-->Y]A", "A[X-->Y]T","A[X-->Y]G", "A[X-->Y]C","A[X-->Y]A")
rownames(Base_Test) <- rows
colnames(Base_Test) <- cols

TestGWTC <- GWTC
TestMutCount <- Mut_Matrix_Total
TestDivisor <- as.matrix(TestGWTC*generations*malines)


num_triplets <- 64
full_col <- 1
for(row in 1:nrow(Sung_Matrix)) 
{
  print(paste("row: ", row, sep = ""))
  for(col in 1:ncol(Sung_Matrix)) 
  {
    #if(col != 1)
    #full_col <- full_col+3
    GWTC_count <- as.numeric(GWTC_Current[row, col])
    print(paste("GWTC: ", GWTC_Current[row,col], sep = ""))
    print(paste("generations: ", generations, sep = ""))
    print(paste("MA lines: ", malines, sep = ""))
    print(paste("Test divisor (product of the 3): ",TestDivisor[row, col]))
    
 #   if((malines*GWTC_Current[row,col]*generations) == as.numeric(TestDivisor[row, col])) 
 #    {
    print(paste("Row: ", row, " Column: ", full_col, sep = ""))
    print(paste( " Mutations: ", TestMutCount[row, col], " Test Divisor: ", TestDivisor[row, col], sep = ""))
    Base_Test[row, full_col] <- TestMutCount[row, full_col]/TestDivisor[row, col]
    print(paste("Mutation Rate: ", Base_Test[row, full_col], sep = "" ))
    
    print(paste("Row: ", row, " Column: ", full_col+1, sep = ""))
    print(paste( "Mutations: ", TestMutCount[row, col+1], " Test Divisor: ", TestDivisor[row, col], sep = ""))
    Base_Test[row, full_col+1] <- TestMutCount[row, full_col+1]/TestDivisor[row, col]
    print(paste("Mutation Rate: ", Base_Test[row, full_col+1], sep = "" ))
    
    print(paste("Row: ", row, " Column: ", full_col+2, sep = ""))
    print(paste( " Mutations: ", TestMutCount[row, col+2], " Test Divisor: ", TestDivisor[row, col], sep = ""))
    Base_Test[row, full_col+2] <- TestMutCount[row, full_col+2]/TestDivisor[row, col]
    print(paste("Mutation Rate: ", Base_Test[row, full_col+2], sep = "" ))
    
    
 #   }
    full_col <- full_col+3
  }
  full_col <- 1 #reset for the next row calculations until all rows calculated.
}

zz <- file("verifyrates.txt","w")
num_triplets <- 64
full_col <- 1
sink(zz)
for(row in 1:nrow(Sung_Matrix)) 
{
  print(paste("row: ", row, sep = ""))
  for(col in 1:ncol(Sung_Matrix)) 
  {
    GWTC_count <- as.numeric(GWTC_Current[row, col])
    print(paste("GWTC: ", GWTC_Current[row,col], sep = ""))
    print(paste("generations: ", generations, sep = ""))
    print(paste("MA lines: ", malines, sep = ""))
    print(paste("Test divisor (product of the 3): ",TestDivisor[row, col]))
    
    if((malines*GWTC_Current[row,col]*generations) == as.numeric(TestDivisor[row, col]))
    {
      print(paste("Row: ", row, " Column: ", full_col, sep = ""))
      print(paste( " Mutations: ", TestMutCount[row, col], " Test Divisor: ", TestDivisor[row, col], sep = ""))
      print(paste("Mutation Rate: ", Base_Test[row, full_col], sep = "" ))
      
      print(paste("Row: ", row, " Column: ", full_col+1, sep = ""))
      print(paste( "Mutations: ", TestMutCount[row, col+1], " Test Divisor: ", TestDivisor[row, col], sep = ""))
      print(paste("Mutation Rate: ", Base_Test[row, full_col+1], sep = "" ))
      
      print(paste("Row: ", row, " Column: ", full_col+2, sep = ""))
      print(paste(" Mutations: ", TestMutCount[row, col+2], " Test Divisor: ", TestDivisor[row, col], sep = ""))
      print(paste("Mutation Rate: ", Base_Test[row, full_col+2], sep = "" ))
      print("")
      print("")
      full_col <- full_col+3
    }
  }
  full_col <- 1 #reset for the next row calculations until all rows calculated.
}

sink()
close(zz)

setwd(Path_MutationRates_output)
write.csv(Base_Test, "Scratch16x12.csv")