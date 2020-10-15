#DirectoryCheck.r

#Script that determines whether I am working at home or in the cubicle and constructs the strings variables to access files and directories

#organism <- "test"
#setwd("/Users/triadge/desktop/CDMAP/Test_Datasets/")
#setwd("/Users/dpatto12/Desktop/CDMAP/Test_Datasets/")

Path_wd <- "/Desktop/CDMAP/Test_Datasets"
#Path_to_scripts <- "/Desktop/CDMAP/R_workspace/Context_Project/Context_Pipeline_0.1.0"
Path_to_scripts <- "/Desktop/CDMAP/CDMAP_Library"

#CHANGE THESE TO SWAP DIRECTORIES BASED ON LOCATION, AND TO OUTPUT OUTSIDE OF Test_Datasets
Path_output <- "/Desktop/CDMAP_Output/"
Path_output_organism <- paste(Path_output, "/", organism, sep = "")
Path_correlate_repo <- "/Desktop/CDMAP_Output/Correlation_Repository"

#4fold mutation debug directories
Path_output_organism <- paste(Path_output, "/", organism, "/", sep = "")
Path_output_triplet <- paste(Path_output_organism, "Triplet", sep = "")
Path_output_upstream <- paste(Path_output_organism, "Upstream", sep = "")
Path_output_downstream <- paste(Path_output_organism, "Downstream", sep = "")

Path_correlate_triplet  <- paste(Path_correlate_repo, "triplet", sep = "/")
Path_correlate_repo_down  <- paste(Path_correlate_repo, "Downstream", sep = "/")
Path_correlate_repo_up  <- paste(Path_correlate_repo, "Upstream", sep = "/")

#### DEBUGGING AND DATASET TESTING DIRECTORIES########
#Path_output <- "/Desktop/CDMAP/Test_Datasets/Scratch_Output"
#Path_correlate_repo <- "/Desktop/CDMAP/Test_Datasets/Scratch_Correlation"
#Path_correlate_triplet <- "/Desktop/CDMAP/Test_Datasets/Scratch_Correlation/triplet"
#Path_correlate_repo_down <- "/Desktop/CDMAP/Test_Datasets/Scratch_Correlation/Downstream"
#Path_correlate_repo_up <- "/Desktop/CDMAP/Test_Datasets/Scratch_Correlation/Upstream"

Path_RefFile <- "/Desktop/CDMAP/Test_Datasets/VF/chr2/CP000021.fasta"
Path_GBFile <- "/Desktop/CDMAP/Test_Datasets/VF/chr2/CP000021.gbk"
Path_InFile <-"/Desktop/CDMAP/Test_Datasets/VF/chr2/mutant_chr2_basecalls.csv"
Path_Refseq <- "/Desktop/CDMAP/Test_Datasets/VF/chr2/CP000021.fasta"

######################################################

if(DirCheck == 'home')
{
  Path_home <- "/Users/triadge"
  Path_wd <- paste(Path_home, Path_wd, sep = "")
  Path_to_scripts <- paste(Path_home, Path_to_scripts, sep = "")
  Path_output <- paste(Path_home, Path_output, sep = "")
  Path_output_organism <- paste(Path_home, Path_output_organism, sep = "")
  Path_output_triplet <- paste(Path_home, Path_output_triplet, sep = "")
  Path_output_upstream <- paste(Path_home, Path_output_upstream, sep = "")
  Path_output_downstream <- paste(Path_home, Path_output_downstream, sep = "")
  Path_RefFile <- paste(Path_home, Path_RefFile, sep = "")
  Path_GBFile <- paste(Path_home, Path_GBFile, sep = "")
  Path_InFile <- paste(Path_home, Path_InFile, sep = "")
  Path_Refseq <- paste(Path_home, Path_Refseq, sep = "")
  Path_correlate_repo <- paste(Path_home, Path_correlate_repo, sep = "")
  Path_correlate_triplet <- paste(Path_home, Path_correlate_triplet, sep = "")
  Path_correlate_repo_up <- paste(Path_home, Path_correlate_repo_up, sep = "")
  Path_correlate_repo_down <- paste(Path_home, Path_correlate_repo_down, sep = "")
}

if(DirCheck == 'work')
{  
  Path_work <- "/Users/dpatto12"
  Path_wd <- paste(Path_work, Path_wd, sep = "")
  Path_to_scripts <- paste(Path_work, Path_to_scripts, sep = "")
  Path_output <- paste(Path_work, Path_output, sep = "")
  Path_output_organism <- paste(Path_work, Path_output_organism, sep = "")
  Path_output_triplet <- paste(Path_work, Path_output_triplet, sep = "")
  Path_output_upstream <- paste(Path_work, Path_output_upstream, sep = "")
  Path_output_downstream <- paste(Path_work, Path_output_downstream, sep = "")
  Path_RefFile <- paste(Path_work, Path_RefFile, sep = "")
  Path_GBFile <- paste(Path_work, Path_GBFile, sep = "")
  Path_InFile <- paste(Path_work, Path_InFile, sep = "")
  Path_Refseq <- paste(Path_work, Path_Refseq, sep = "")
  Path_correlate_repo <- paste(Path_work, Path_correlate_repo, sep = "")
  Path_correlate_triplet <- paste(Path_work, Path_correlate_triplet, sep = "")
  Path_correlate_repo_up <- paste(Path_work, Path_correlate_repo_up, sep = "")
  Path_correlate_repo_down <- paste(Path_work, Path_correlate_repo_down, sep = "")
}

