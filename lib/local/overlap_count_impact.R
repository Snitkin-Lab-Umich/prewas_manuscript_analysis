# Goal: find all of the SNPs in overlapping positions per dataset. 

library(tidyverse)
source("lib/local/multiallelic_lib.R")
num_project <- nrow(project_df)

# Get overlap info -------------------------------------------------------------
for (i in 1:num_project) {
  print(project_df$projects[i])
  snpmat <- load_parsed_snpmat_by_project(project_df$projects[i])
  if (!is.null(snpmat)) {
    overlap_logical <- snpmat$allele$annots$rows_with_overlapping_genes_log
    overlap_pos <- snpmat$allele$annots$pos[overlap_logical]
    overlap_snpeff <-
      as.character(snpmat$allele$annots$snpeff_impact[overlap_logical])
    summary_snpeff_df <- cbind(overlap_pos, overlap_snpeff)
    summary_snpeff_df <- as.data.frame(summary_snpeff_df)
    summary_snpeff_df <- dcast(summary_snpeff_df, overlap_pos ~ overlap_snpeff)
    summary_snpeff_df <- summary_snpeff_df %>% 
      select(HIGH, MODERATE, LOW) %>% 
      filter(rowSums(.) > 1)

    if (nrow(summary_snpeff_df) > 0) {
      write_tsv(x = summary_snpeff_df, 
                col_names = TRUE, 
                path = paste0("data/local/overlap/",
                              project_df$projects[i], 
                              "_SNPEFF_annotation_for_overlap_sites.tsv"))
    } else {
      print("No overlapping sites with 2 non-modifier annotations")
    }
  }
}

# Summarize overlap info -------------------------------------------------------
data_dir <- "data/local/overlap/"
names <- c("Project", 
           "NumOverlapSites", 
           "NumIdenticalImpact", 
           "NumDifferentImpact", 
           "NumBadAnnot", 
           "NumHighHigh", 
           "NumHighModerate",
           "NumHighLow", 
           "NumModerateModerate", 
           "NumModerateLow", 
           "NumLowLow")

summary <- as.data.frame(matrix(NA, 
                                ncol = length(names),
                                nrow = nrow(project_df)))
colnames(summary) <- names 
for (i in 1:num_project) {
  filename <- paste0(data_dir, 
                     project_df$projects[i], 
                     "_SNPEFF_annotation_for_overlap_sites.tsv")
  if (file.exists(filename)) {
    temp_mat <- read_tsv(filename)
    temp_mat <- temp_mat > 0
    summary$Project[i] <- project_df$projects[i]
    summary$NumOverlapSites[i] <- nrow(temp_mat)
    num_identical <- num_different <- bad <- 0
    hi_hi <- hi_mo <- hi_lo <- mo_mo <- mo_lo <- lo_lo <- 0
    for (j in 1:nrow(temp_mat)) {
      num_identical <- num_identical + (sum(temp_mat[j, ] == 0) == 2)
      num_different <- num_different + (sum(temp_mat[j, ] == 0) < 2)
      bad <- bad +  (sum(temp_mat[j, ] == c(0, 0, 0)) == 3)
      hi_hi <- hi_hi + (sum(temp_mat[j, ] == c(1, 0, 0)) == 3)
      hi_mo <- hi_mo + (sum(temp_mat[j, ] == c(1, 1, 0)) == 3)
      hi_lo <- hi_lo + (sum(temp_mat[j, ] == c(1, 0, 1)) == 3)
      mo_mo <- mo_mo + (sum(temp_mat[j, ] == c(0, 1, 0)) == 3)
      mo_lo <- mo_lo + (sum(temp_mat[j, ] == c(0, 1, 1)) == 3)
      lo_lo <- lo_lo + (sum(temp_mat[j, ] == c(0, 0, 1)) == 3)
    }
    summary$NumIdenticalImpact[i] <- num_identical
    summary$NumDifferentImpact[i] <- num_different
    summary$NumBadAnnot[i] <- bad
    summary$NumHighHigh[i] <- hi_hi
    summary$NumHighModerate[i] <- hi_mo    
    summary$NumHighLow[i] <- hi_lo   
    summary$NumModerateModerate[i] <- mo_mo    
    summary$NumModerateLow[i] <- mo_lo
    summary$NumLowLow[i] <- lo_lo
  }
  summary <- summary %>% 
    mutate(percent_different = 100 * NumDifferentImpact / NumOverlapSites) 
}

write_tsv(x = summary,
          path = "data/snpeff_overlap_summary.tsv", 
          col_names = TRUE)


## SAVE INFO ON OVERLAPPING POS NUMBER AND OVERLAPPING GENE NUMBER
overlap_stats <- as.data.frame(matrix(NA, nrow = num_project, ncol = 3))
colnames(overlap_stats) <- c("Project", "Overlapping_SNP_loci", "Overlapping_Genes_w_SNPs")
for (i in 1:num_project) {
  print(project_df$projects[i])
  snpmat <- load_parsed_snpmat_by_project(project_df$projects[i])
  if (!is.null(snpmat)) {
    overlap_logical <- snpmat$allele$annots$rows_with_overlapping_genes_log
    num_unique_overlap_SNP_pos <- length(unique(snpmat$allele$annots$pos[overlap_logical]))
    num_unique_overlap_genes_with_SNP <- length(unique(snpmat$allele$annots$locus_tag[overlap_logical]))

    overlap_stats$Project[i] <- project_df$projects[i] 
    overlap_stats$Overlapping_SNP_loci[i] <- num_unique_overlap_SNP_pos
    overlap_stats$Overlapping_Genes_w_SNPs[i] <- num_unique_overlap_genes_with_SNP
  }
}

write_tsv(x = overlap_stats, 
          col_names = TRUE, 
          path = "data/overlap_stats.tsv")

