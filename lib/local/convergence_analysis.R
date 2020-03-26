# This script takes in summary data generated on the HPC regarding the number
# of variant transitions  (WT -> Variant or Variant -> WT) for both biallelic 
# and multiallelic sites for each data set. Given this data, it computes the
# number of times convergence occurs. Data analyzed in this script are used for
# Figure S2D. 

data_dir <- "data/hpc/convergence/"
source("lib/local/multiallelic_lib.R")
bi_mat <- "_biallelic_transition_matrix.csv"
multi_mat <- "_multiallelic_transition_matrix.csv"
mat <- "_transition_matrix.csv"
num_projects <- nrow(project_df)
num_columns <- 6
summary <- as.data.frame(matrix(NA, nrow = 0, ncol = num_columns))
sum_col_names <- c("Project",
                   "SNP_Type", 
                   "Total_Num_SNPs", 
                   "Num_Times_SNP_Arises", 
                   "Count_SNP_Arises",
                   "relative_frequency")
colnames(summary) <- sum_col_names

for (i in 1:num_projects) {
  bi_mat_path <- paste0(data_dir, project_df$projects[i], bi_mat)
  multi_mat_path <- paste0(data_dir, project_df$projects[i], multi_mat)
  
  if (file.exists(bi_mat_path)) {
    print(bi_mat_path)
    temp_bi_mat <-
      read.csv(bi_mat_path, stringsAsFactors = FALSE, row.names = 1)
    temp_multi_mat <- 
      read.csv(multi_mat_path, stringsAsFactors = FALSE, row.names = 1)
    for (j in 1:length(table(temp_bi_mat))) {
      temp_summary <- 
        as.data.frame(matrix(c(project_df$projects[i], 
                               "biallelic", 
                               nrow(temp_bi_mat),
                               names(table(temp_bi_mat)[j]), 
                               unname(table(temp_bi_mat)[j]), 
                               100 * (unname(table(temp_bi_mat)[j]) /
                                        nrow(temp_bi_mat))), 
                             nrow = 1, 
                             ncol = num_columns)) 
      colnames(temp_summary) <- sum_col_names
      summary <- rbind(summary, temp_summary)
    }
    
    for (k in 1:length(table(temp_multi_mat))) {
      temp_summary <- 
        as.data.frame(matrix(c(project_df$projects[i], 
                               "multiallelic", 
                               nrow(temp_multi_mat),
                               names(table(temp_multi_mat)[k]), 
                               unname(table(temp_multi_mat)[k]), 
                               100 * (unname(table(temp_multi_mat)[k]) /
                                        nrow(temp_multi_mat))), 
                             nrow = 1, 
                             ncol = num_columns)) 
      colnames(temp_summary) <- sum_col_names
      summary <- rbind(summary, temp_summary)
    }
  }
}

summary <- summary %>% 
  mutate(Num_Convergence_Events = 
           (as.numeric(as.character(Num_Times_SNP_Arises)) - 1), 
         Num_Convergence_Events_8_bin = Num_Convergence_Events)

summary$Num_Convergence_Events_8_bin[
  summary$Num_Convergence_Events_8_bin > 7] <- NA
summary$relative_frequency <-
  as.numeric(as.character(summary$relative_frequency))
summary$Num_Times_SNP_Arises <- 
  as.numeric(as.character(summary$Num_Times_SNP_Arises))
summary <- 
  summary %>% mutate(Num_Convergence_Events = Num_Times_SNP_Arises - 1) 
summary <- 
  summary %>% mutate(Num_Convergence_Events_8_bin = Num_Convergence_Events)
summary$Num_Convergence_Events_8_bin[
  summary$Num_Convergence_Events_8_bin > 7] <- 8
summary$Num_Convergence_Events_8_bin <-
  as.character(summary$Num_Convergence_Events_8_bin)
summary$Num_Convergence_Events_8_bin[
  summary$Num_Convergence_Events_8_bin == "8"] <- ">7"
summary$Num_Convergence_Events_8_bin <-
  as.character(summary$Num_Convergence_Events_8_bin)
summary$Num_Convergence_Events_8_bin[
  is.na(summary$Num_Convergence_Events_8_bin)] <- ">7"

summary$Num_Convergence_Events_8_bin <-
  as.factor(summary$Num_Convergence_Events_8_bin)
summary <- summary %>% mutate("Locus Type" = SNP_Type)

write_tsv(summary, "data/local/convergence/convergence.tsv")


multi_dist <- summary %>% 
  filter(SNP_Type == "multiallelic") %>%
  select(Num_Convergence_Events) %>% 
  filter(Num_Convergence_Events >= 0) %>% 
  unlist() %>% 
  unname()
bi_dist <- summary %>% 
  filter(SNP_Type == "biallelic") %>%
  filter(Num_Convergence_Events >= 0) %>% 
  select(Num_Convergence_Events) %>% 
  unlist() %>% 
  unname()


ks_results <- ks.test(bi_dist, multi_dist)


num_test <- 100
pvals <- rep(NA, num_test)
set.seed(1)
for (i in 1:num_test) {
  sub_bi_dist <- sample(x = bi_dist, size = length(multi_dist), replace = FALSE)
  sub_ks_results <- ks.test(sub_bi_dist, multi_dist)
  pvals[i] <- sub_ks_results$p.value
}

subset_sig_percent <- 100 * sum(pvals < 0.05) / num_test

results_mat <- matrix(NA, nrow = 3)
rownames(results_mat) <- 
  c("KS D-stat", 
    "KS P-value", 
    "% of Subsampled KS tests that are significant P-value")
results_mat[1, 1] <- as.character(round(ks_results$statistic, 3))
results_mat[2, 1] <- as.character(round(ks_results$p.value, 10))
results_mat[3, 1] <- as.character(round(subset_sig_percent, 3))

write.csv(results_mat, file = "data/local/convergence/convergence_ks_test.csv")
# The two  distributions ARE significantly different