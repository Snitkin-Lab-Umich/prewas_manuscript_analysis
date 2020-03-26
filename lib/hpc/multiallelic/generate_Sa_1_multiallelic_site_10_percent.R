# Sa_1 10 % of total samples
source('../lib/multiallelic_lib.R')
source('../lib/subsample_lib.R')
# Initialize data
num_rep <- 10
percent_samples <- 10
current_project <- "Sa_1"
current_species <- "S. aureus"
data <- generate_nrow_df(num_rep)
current_project <- "Sa_1"
mat <- load_project_snpmat(current_project)
num_to_select <- round(ncol(mat) * percent_samples/100, 0)
for (i in 1:num_rep) {
  sub_mat <- subsample_snpmat(mat, num_to_select, seed = i)
  sub_mat_vecs <- keep_only_variant_sites(sub_mat)
  data <- fill_in_data(data_mat = data,
                       mat = sub_mat_vecs$mat,
                       row_num = i,
                       species_genome_size = 2800000,
                       bi_vec = sub_mat_vecs$biallelic_vec,
                       tri_vec = sub_mat_vecs$triallelic_vec,
                       quad_vec = sub_mat_vecs$quadallelic_vec,
                       multi_vec = sub_mat_vecs$multiallelic_vec)
data$Project[i] <- current_project
data$Species[i] <- current_species
write_tsv(data,
          path = paste0('../data/multiallelic_variant_summary_',
          current_project,
          '_',
          percent_samples,
          '_percent',
          '.tsv'))
}
