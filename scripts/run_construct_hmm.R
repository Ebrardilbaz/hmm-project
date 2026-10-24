library(Biostrings)

setwd("/cta/users/ebrar.dilbaz/hmm-project")

# Helper fonksiyonlar
aa_to_num <- function(aa) {
  amino_acids <- c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a)==1, as.numeric(which(amino_acids %in% a)), 21)})
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T")
  aa <- ifelse(num==21, 21, amino_acids[num])
  return(aa)
}

source("scripts/construct_HMM.R")

construct_HMM(
  folder_name = "results/hmm_profiles",
  main_fam = "GRM1",
  close_fam = "GRM5",
  rest_fam = c("GRM2", "GRM4"),
  cfunvals_file = "../identity_scores_grm1.csv",
  bound1_conslevel = 0.98,
  bound2_conslevel = 0.8,
  bound_blosum = 0.5
)

cat("Done! GRM1_UpdatedProfiles.hmm created.\n")
