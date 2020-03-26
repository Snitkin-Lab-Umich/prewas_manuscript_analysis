# Sa_1
source('../lib/multiallelic_lib.R')
# Initialize data
num_rep <- 10
current_project <- "Sa_1"
current_species <- "S. aureus"
data <- generate_single_df()
mat <- load_project_snpmat(current_project)
mat_vecs <- keep_only_variant_sites(mat)
data <- fill_in_data(data_mat = data,
                     mat = mat_vecs$mat,
                     row_num = one_row,
                     species_genome_size = 2800000,
                     bi_vec = mat_vecs$biallelic_vec,
                     tri_vec = mat_vecs$triallelic_vec,
                     quad_vec = mat_vecs$quadallelic_vec,
                     multi_vec = mat_vecs$multiallelic_vec)
data$Project[one_row] <- current_project
data$Species[one_row] <- current_species
write_tsv(data,
          path = paste0('../data/multiallelic_variant_summary_',
          current_project,
          '.tsv'))
