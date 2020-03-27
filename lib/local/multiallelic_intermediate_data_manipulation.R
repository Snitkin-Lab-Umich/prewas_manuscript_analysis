# This script takes in data generated on the HPC and summarizes and reformats it
# here. The goal of this script is to describe the prevalence of multiallelic
# sites in the bacterial data. Data from these analyses are plotted in Figures
# 2 and S2. 

# Read in all multiallelic summary data (full and subsampled, store nicely here)
source("lib/local/multiallelic_lib.R")
data_dir <- "data/hpc/multiallelic/"
percent_samples <- c(1, 2.5, 5, 7.5, 10, 12.5, 25, 50, 75)
num_col <- 23
project_key <- read_tsv("data/key/project_names_and_colors.tsv")
num_project <- nrow(project_key)
# Read in multiallelic info ----------------------------------------------------
variant_summary <- as.data.frame(matrix(NA, ncol = num_col, nrow = 0))
for (i in 1:length(project_df$projects)) {
  filename <-
    paste0(data_dir,
           "multiallelic_variant_summary_",
           project_df$projects[i],
           ".tsv")
  if (file.exists(filename)) {
    temp_mat <- read_tsv(filename)
    temp_mat <- temp_mat %>% mutate("Percent" = 100)
    variant_summary <- rbind(variant_summary, temp_mat)
  }
}

subsampled_mat <- as.data.frame(matrix(NA, ncol = num_col, nrow = 0))
for (i in 1:length(project_df$projects)) {
  for (j in 1:length(percent_samples)) {
    filename <- paste0(data_dir,
                       "multiallelic_variant_summary_",
                       project_df$projects[i],
                       "_",
                       percent_samples[j],
                       "_percent.tsv")
    if (file.exists(filename)) {
      temp_mat <- read_tsv(filename)
      temp_mat <- temp_mat %>% mutate("Percent" = percent_samples[j])
      subsampled_mat <- rbind(subsampled_mat, temp_mat)
    }
  }
}
colnames(subsampled_mat) <- colnames(temp_mat)

# combine all of the data together here:
snp_mat_summary <- rbind(variant_summary, subsampled_mat)
snp_mat_summary <- snp_mat_summary %>%
  filter(rowSums(is.na(.)) < 10)


# Group by project -------------------------------------------------------------
group_by_project <-
  snp_mat_summary %>%
  group_by(Project, Percent) %>%
  summarise(AvgNumMultiallelicSite = mean(NumMultiallelicSite),
            AvgNumVariantSite = mean(NumVariantSite),
            Species = unique(Species),
            MaxSNPDist = mean(MaxSNPDist),
            MeanSNPDist = mean(MeanSNPDist),
            MedianSNPDist = mean(MedianSNPDist)) %>%
  select(AvgNumVariantSite,
         AvgNumMultiallelicSite,
         Species,
         Percent,
         MeanSNPDist,
         MedianSNPDist,
         MaxSNPDist,
         Project)

group_by_project <- left_join(group_by_project, project_key, by = "Project")

write_tsv(group_by_project,
          path = "data/local/multiallelic/multiallelic_summary_by_project.tsv")

# Multi to any variant ratio calculation ---------------------------------------
variant_summary <- variant_summary %>%
  mutate(Baseline = NumMultiallelicSite / NumVariantSite,
         Low = NumMultiallelicLow / NumVariantLow,
         Moderate = NumMultiallelicModerate / NumVariantModerate,
         High = NumMultiallelicHigh / NumVariantHigh)

variant_summary <- left_join(variant_summary, project_key)
variant_summary <- variant_summary %>%
  mutate("Dataset Size (#)" = NumIsolates)
snp_mat_summary <- snp_mat_summary %>%
  mutate("Dataset Size (#)" = NumIsolates)
snp_mat_summary <- left_join(snp_mat_summary, project_key)
write_tsv(variant_summary, "data/local/multiallelic/multiallelic_summary.tsv")
write_tsv(snp_mat_summary, path = "data/local/multiallelic/multiallelic_summary_subsampled.tsv")

# Melt the variant summary
melted_diversity_variant <- melted_variant_summary <-
  melt(variant_summary,
       id.vars = 'Project',
       measure.vars = c('Baseline', 'Low', 'Moderate', 'High'))

write_tsv(melted_diversity_variant, "data/local/multiallelic/SNP_diversity_long.tsv")

# Calculate significance, functional impact by low vs hi diversity -------------
t_test_results <- matrix(NA, nrow = 3, ncol = 2)
colnames(t_test_results) <- c("pvalue", "adj_pvalue")
row.names(t_test_results) <- c("Low", "Moderate", "High")

for (i in 1:nrow(t_test_results)) {
  t_result <-
    t.test(melted_diversity_variant %>%
             filter(variable == "Baseline") %>%
             select(value) %>%
             unlist() %>%
             unname(),
           melted_diversity_variant %>%
             filter(variable == row.names(t_test_results)[i]) %>%
             select(value) %>%
             unlist() %>%
             unname())
  t_test_results[i, 1] <- t_result$p.value
}

t_test_results[, 2] <-
  p.adjust(t_test_results[, 1], method = "bonf")

t_test_results <-
  rownames_to_column(as.data.frame(t_test_results),
                     var = "functional_impact")

write_tsv(t_test_results,
          "data/local/multiallelic/t_test_multiallelic.tsv")

# Functional impact at multiallelic sites --------------------------------------
for (i in 1:num_project) {
  print(project_key$Project[i])
  snpmat <- load_project_snpmat(project_key$Project[i]) # read in unsplit, no ig
  if (!is.null(snpmat)) {
    var_only_snps <- keep_only_variant_sites(snpmat)
    multi_var_snp_mat <-
      var_only_snps$mat[var_only_snps$multiallelic_vec, , drop = FALSE]
    multi_var_snpeff_annots_mat <-
      matrix(0, nrow = nrow(multi_var_snp_mat), ncol = 3)
    colnames(multi_var_snpeff_annots_mat) <- c("HIGH", "MODERATE", "LOW")
    multi_var_snpeff_annots_mat <- as.data.frame(multi_var_snpeff_annots_mat)
    for (j in 1:nrow(multi_var_snp_mat)) {
      multi_var_snpeff_annots_mat$HIGH[j] <-
        str_count(row.names(multi_var_snp_mat)[j], "HIGH")
      multi_var_snpeff_annots_mat$MODERATE[j] <-
        str_count(row.names(multi_var_snp_mat)[j], "MODERATE")
      multi_var_snpeff_annots_mat$LOW[j] <-
        str_count(row.names(multi_var_snp_mat)[j], "LOW")
    }
    write_tsv(x = multi_var_snpeff_annots_mat,
              col_names = TRUE,
              path = paste0("data/local/multiallelic/",
                            project_key$Project[i],
                            "_SNPEFF_annotation_for_multiallelic_sites.tsv"))

    tri_var_snpeff_annots_mat <-
      multi_var_snpeff_annots_mat[
        rowSums(multi_var_snpeff_annots_mat) == 2, , drop = FALSE]

    print("multi then tri dim")
    print(dim(multi_var_snpeff_annots_mat))
    print(dim(tri_var_snpeff_annots_mat))

    if (ncol(tri_var_snpeff_annots_mat) > 0)
    write_tsv(x = tri_var_snpeff_annots_mat,
              col_names = TRUE,
              path = paste0("data/local/multiallelic/",
                            project_key$Project[i],
                            "_SNPEFF_annotation_for_triallelic_sites.tsv"))
  }
}

data_dir <- "data/local/multiallelic/"
snpeff_summary <- as.data.frame(matrix(NA, ncol = 11, nrow = num_project))
colnames(snpeff_summary) <- c("project",
                              "NumMultiallelicSites",
                              "NumIdenticalImpact",
                              "NumDifferentImpact",
                              "NumHighAndModerate",
                              "NumHighAndLow",
                              "NumModerateAndLow",
                              "NumBadAnnot",
                              "OnlyHigh",
                              "OnlyModerate",
                              "OnlyLow")

for (i in 1:num_project) {
  filename <- paste0(data_dir,
                     project_key$Project[i],
                     "_SNPEFF_annotation_for_multiallelic_sites.tsv")
  if (file.exists(filename)) {
    print(project_key$Project[i])
    temp_mat <- read_tsv(filename)
    temp_mat <- temp_mat > 0
    snpeff_summary$project[i] <- project_key$Project[i]
    snpeff_summary$NumMultiallelicSites[i] <- nrow(temp_mat)
    num_identical <- num_different <- hi_moderate <- hi_low <- 0
     bad <- hi <- moderate <- low <- moderate_low <- 0
    for (j in 1:nrow(temp_mat)) {
      num_identical <- num_identical + (sum(temp_mat[j, ] == 0) == 2)
      num_different <- num_different + (sum(temp_mat[j, ] == 0) < 2)
      hi_moderate <- hi_moderate + (sum(temp_mat[j, ] == c(1, 1, 0)) == 3)
      hi_low <- hi_low + (sum(temp_mat[j, ] == c(1, 0, 1)) == 3)
      moderate_low <- moderate_low + (sum(temp_mat[j, ] == c(0, 1, 1)) == 3)
      bad <- bad +  (sum(temp_mat[j, ] == c(0, 0, 0)) == 3)
      hi <- hi + (sum(temp_mat[j, ] == c(1, 0, 0)) == 3)
      moderate <- moderate + (sum(temp_mat[j, ] == c(0, 1, 0)) == 3)
      low <- low + (sum(temp_mat[j, ] == c(0, 0, 1)) == 3)
    }
    snpeff_summary$NumIdenticalImpact[i] <- num_identical
    snpeff_summary$NumDifferentImpact[i] <- num_different
    snpeff_summary$NumHighAndModerate[i] <- hi_moderate
    snpeff_summary$NumHighAndLow[i] <- hi_low
    snpeff_summary$NumModerateAndLow[i] <- moderate_low
    snpeff_summary$NumBadAnnot[i] <- bad
    snpeff_summary$OnlyHigh[i] <- hi
    snpeff_summary$OnlyModerate[i] <- moderate
    snpeff_summary$OnlyLow[i] <- low
  }
}

snpeff_summary <- snpeff_summary %>%
  mutate(percent_different = 100 * NumDifferentImpact / NumMultiallelicSites)
colnames(snpeff_summary)[colnames(snpeff_summary) == "project"] <- "Project"
write_tsv(x = snpeff_summary,
          path = "data/local/multiallelic/snpeff_multiallelic_summary.tsv",
          col_names = TRUE)

# For the heatmap figure I'm trying to make, make a specific summary of the
# triallelic data
data_dir <- "data/"

tri_summary <- as.data.frame(matrix(NA, ncol = 11, nrow = num_project))
colnames(tri_summary) <- c("project",
                           "NumMultiallelicSites",
                           "NumIdenticalImpact",
                           "NumDifferentImpact",
                           "NumBadAnnot",
                           "NumHighHigh",
                           "NumHighModerate",
                           "NumHighLow",
                           "NumModerateModerate",
                           "NumModerateLow",
                           "NumLowLow")
for (i in 1:num_project) {
  filename <- paste0(data_dir,
                     project_key$Project[i],
                     "_SNPEFF_annotation_for_triallelic_sites.tsv")
  if (file.exists(filename)) {
    print(project_key$Project[i])
    temp_mat <- read_tsv(filename)
    temp_mat <- temp_mat > 0
    tri_summary$project[i] <- project_key$Project[i]
    tri_summary$NumMultiallelicSites[i] <- nrow(temp_mat)
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
    tri_summary$NumIdenticalImpact[i] <- num_identical
    tri_summary$NumDifferentImpact[i] <- num_different
    tri_summary$NumBadAnnot[i] <- bad
    tri_summary$NumHighHigh[i] <- hi_hi
    tri_summary$NumHighModerate[i] <- hi_mo
    tri_summary$NumHighLow[i] <- hi_lo
    tri_summary$NumModerateModerate[i] <- mo_mo
    tri_summary$NumModerateLow[i] <- mo_lo
    tri_summary$NumLowLow[i] <- lo_lo
  }
}

tri_summary <- tri_summary %>%
  mutate(percent_different = 100 * NumDifferentImpact / NumMultiallelicSites)
colnames(tri_summary)[colnames(tri_summary) == "project"] <- "Project"
write_tsv(x = tri_summary,
          path = "data/local/multiallelic/snpeff_triallelic_summary.tsv",
          col_names = TRUE)

tri_htmp_data <- tri_summary %>%
  filter(!Project %in% drop_projects) %>%
  select(NumHighHigh,
         NumHighModerate,
         NumHighLow,
         NumModerateModerate,
         NumModerateLow,
         NumLowLow) %>%
  colSums(na.rm = TRUE)

tri_htmp <- matrix(0, ncol = 3, nrow = 3)
colnames(tri_htmp) <- row.names(tri_htmp) <- c("High", "Moderate", "Low")
tri_htmp[1, 1] <- tri_htmp_data[1]
tri_htmp[2, 2] <- tri_htmp_data[4]
tri_htmp[3, 3] <- tri_htmp_data[6]

tri_htmp[1, 2] <- tri_htmp[2, 1] <- tri_htmp_data[2]
tri_htmp[1, 3] <- tri_htmp[3, 1] <- tri_htmp_data[3]
tri_htmp[2, 3] <- tri_htmp[3, 2] <- tri_htmp_data[5]

write.csv(tri_htmp, file = "data/local/multiallelic/triallelic_heatmap.csv")
