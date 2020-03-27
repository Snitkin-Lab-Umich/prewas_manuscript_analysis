# This library contains many functions used by seveveral analyses, including
# multiallelic, convergence, and overlap. 

# Libs -------------------------------------------------------------------------
library(tidyverse)
library(reshape2)

# Data needed to run multiallelic stuff ----------------------------------------
one_row <- 1

# Genome sizes by species
Cd_size <- 4300000
Kp_size <- 5400000
Sa_size <- 2800000
Lc_size <- 2040000
Sm_size <- 4850000
Efm_size <- 2900000
Efs_size <- 2900000

project_df <- as.data.frame(matrix(NA, nrow = 9, ncol = 4))
colnames(project_df) <- c("projects", "species", "species_size", "row")

project_df$projects <-
  c("Cd_3", "Cd_4",
    "Kp_1", 
    "Sa_1", "Sa_2",
    "Efm_1", 
    "Efs_1",
    "Sm_1",
    "Lc_1")
project_df$species_size <-
  c(rep(Cd_size, 2),
    rep(Kp_size, 1),
    rep(Sa_size, 2),
    Efm_size,
    Efs_size,
    Sm_size,
    Lc_size)
project_df$species <-
  c(rep("C. difficile", 2),
    rep("K. pneumoniae", 1),
    rep("S. aureus", 2),
    "E. faecium",
    "E. faecalis",
    "S. maltophilia",
    "L. crispatus")
project_df$row <- 1:nrow(project_df)

snp_mat_dir <- "data/snp_matrices_and_parsed_data/"

data_col_names <- c("Project",
                    "Species",
                    "NumIsolates",
                    "AvgGenomeBP",
                    "TotalBP",
                    "NumVariantSite",
                    "NumBiallelicSite",
                    "NumTriallelicSite",
                    "NumQuadallelicSite",
                    "NumMultiallelicSite",
                    "NumVariantLow",
                    "NumBiallelicLow",
                    "NumMultiallelicLow",
                    "NumVariantModerate",
                    "NumBiallelicModerate",
                    "NumMultiallelicModerate",
                    "NumVariantHigh",
                    "NumBiallelicHigh",
                    "NumMultiallelicHigh",
                    "MeanSNPDist",
                    "MedianSNPDist",
                    "MaxSNPDist",
                    "NumBiallelicOverlappingGene")

# Functions --------------------------------------------------------------------
# This function only works for SNP matrices
# Subset a matrix to only variant sites (invariant sites are removed)
keep_only_variant_sites <- function(mat){
  rows_to_keep <- apply(mat, 1, function(row){
    sum(unique(row) %in% c("A", "T", "G", "C"))
  })
  biallelic <-
    multiallelic <-
    triallelic <-
    quadallelic <- rows_to_keep[rows_to_keep > 1]
  biallelic <- as.logical(biallelic == 2)
  triallelic <- as.logical(triallelic == 3)
  quadallelic <- as.logical(quadallelic == 4)
  multiallelic <- as.logical(multiallelic > 2)
  mat <- mat[rows_to_keep > 1, , drop = FALSE]
  
  # Convert N and - to NA
  mat[mat == "N"] <- NA
  mat[mat == "-"] <- NA
  
  print("kept only variant sites")
  results <- list("mat" = mat,
                  "biallelic_vec" = biallelic,
                  "triallelic_vec" = triallelic,
                  "quadallelic_vec" = quadallelic,
                  "multiallelic_vec" = multiallelic)
  return(results)
}

# Calculate  pairwise SNP distance for SNP matrix
calc_snp_dist <- function(mat){
  snp_dist_mat <- apply(mat, 2, function(col){
    colSums(col != mat, na.rm = TRUE)
  })
  print("calculated pairwise snp dist")
  snp_dist <- snp_dist_mat[lower.tri(snp_dist_mat)]
  return(snp_dist)
}  

# Format empty df
generate_storage_df <- function(){
  data <- as_tibble(matrix(NA, ncol = length(data_col_names), nrow = 15))
  colnames(data) <- data_col_names
  return(data)
}

# Format empty df
generate_single_df <- function(){
  data <- as_tibble(matrix(NA, ncol = length(data_col_names), nrow = 1))
  colnames(data) <- data_col_names
  return(data)
}

# Format empty df
generate_nrow_df <- function(num_row){
  data <- as_tibble(matrix(NA, ncol = length(data_col_names), nrow = num_row))
  colnames(data) <- data_col_names
  return(data)
}

# Functions specific to Snitkin Lab made SNP matrices --------------------------
read_snpmat <- function(path){
  mat <- read.table(path,
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    header = TRUE,
                    row.names = 1)
  print("read in snpmat")
  return(mat)
}

# Count the number of overlapping genes given a snpmat
count_overlapping_genes <- function(mat){
  #identify the rows with multiple annotations
  mult_annot <-
    as.logical(sapply(1:nrow(mat),
                      function(x) lengths(regmatches(row.names(mat)[x],
                                                     gregexpr(";[A,C,G,T]",
                                                              row.names(mat)[x])))) > 1)
  
  # identify rows with biallelic variants
  biallelic <- apply(mat, 1, function(row){
    unname(unlist(as.logical(sum(unique(row) %in% c("A", "T", "G", "C")) == 2)))
  })
  
  coding <- as.logical(1:nrow(mat) %in% grep("^Coding SNP", row.names(mat)))
  # remainder of the rows are presumed to be overlapping genes
  num_rows_with_overlapping_genes = sum(mult_annot + biallelic + coding == 3)
  
  print("counted overlapping genes")
  return(num_rows_with_overlapping_genes)
}

# Summarize a lot of multiallelic data
fill_in_data <- function(data_mat,
                         mat,
                         row_num,
                         species_genome_size,
                         bi_vec,
                         tri_vec,
                         quad_vec,
                         multi_vec){
  
  print("started filling out data matrix")
  data_mat$AvgGenomeBP[row_num] <- species_genome_size
  data_mat$NumIsolates[row_num] <- ncol(mat)
  data_mat$TotalBP[row_num] <- ncol(mat) * species_genome_size
  print("filled in basic stats")
  
  data_mat$NumVariantSite[row_num] <- nrow(mat)
  data_mat$NumBiallelicSite[row_num] <- sum(bi_vec)
  data_mat$NumTriallelicSite[row_num] <- sum(tri_vec)
  data_mat$NumQuadallelicSite[row_num] <- sum(quad_vec)
  data_mat$NumMultiallelicSite[row_num] <- sum(multi_vec)
  print("filled in num variant sites")
  
  data_mat$NumVariantLow[row_num] <- length(grep("[|]LOW[|]", row.names(mat)))
  data_mat$NumBiallelicLow[row_num] <-
    length(grep("[|]LOW[|]", row.names(mat)[bi_vec]))
  data_mat$NumMultiallelicLow[row_num] <-
    length(grep("[|]LOW[|]", row.names(mat)[multi_vec]))
  print("filled in LOW")
  
  data_mat$NumVariantModerate[row_num] <-
    length(grep("[|]MODERATE[|]", row.names(mat)))
  data_mat$NumBiallelicModerate[row_num] <-
    length(grep("[|]MODERATE[|]", row.names(mat)[bi_vec]))
  data_mat$NumMultiallelicModerate[row_num] <-
    length(grep("[|]MODERATE[|]", row.names(mat)[multi_vec]))
  print("filled in MODERATE")
  
  data_mat$NumVariantHigh[row_num] <- length(grep("[|]HIGH[|]", row.names(mat)))
  data_mat$NumBiallelicHigh[row_num] <-
    length(grep("[|]HIGH[|]", row.names(mat)[bi_vec]))
  data_mat$NumMultiallelicHigh[row_num] <-
    length(grep("[|]HIGH[|]", row.names(mat)[multi_vec]))
  print("filled in HIGH")
  
  snp_dist <- calc_snp_dist(mat)
  
  data_mat$MeanSNPDist[row_num] <- mean(snp_dist, na.rm = TRUE)
  data_mat$MedianSNPDist[row_num] <- median(snp_dist, na.rm = TRUE)
  data_mat$MaxSNPDist[row_num] <- max(snp_dist, na.rm = TRUE)
  print("filled in snp dist")
  
  data_mat$NumBiallelicOverlappingGene[row_num] <- count_overlapping_genes(mat)
  print("filled in overlapping genes")
  
  print("Done")
  return(data_mat)
}

load_project_snpmat <- function(project_name){
  print("Loading Project:")
  print(project_name)
  mat <- NULL
  temp_file_name <-
    paste0(snp_mat_dir, project_name, "_allele_mat_unsplit_no_ig.tsv")
  if (file.exists(temp_file_name)) {
    mat <- read_snpmat(temp_file_name)
  }
  if (is.null(mat)) {
    print("SNP mat does not exist")
  }
  return(mat)
}

load_parsed_snpmat_by_project <- function(project_name) {
  mat <- NULL
  fname <- paste0('data/hpc/snp_matrices_and_parsed_data/', 
                  project_name,
                  '_parsed.RData')
  if (file.exists(fname)) {
    mat <- local(get(load(fname)))
  } else {
    print("File not found")
  }
  return(mat)
}
