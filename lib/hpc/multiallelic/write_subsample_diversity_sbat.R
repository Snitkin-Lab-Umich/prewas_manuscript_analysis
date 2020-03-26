# This script generates the subsampling of the various SNP matrics to that we
# can compare the relationship between multialleleic sites and isolates number.

source('../lib/multiallelic_lib.R')
percent_samples <- c(1, 2.5, 5, 7.5, 10, 12.5, 25, 50, 75)

for (project in 1:length(project_df$projects)) {
  for (percent in 1:length(percent_samples)) {
    fname <- paste0("../lib/generate_",
                    project_df$projects[project], 
                    "_multiallelic_site_", 
                    percent_samples[percent], 
                    "_percent.R")
    fname <- gsub("[.]5", "p5", x = fname)
    
    project_specific_code <-
      paste(paste0('current_project <- ', '"', project_df$projects[project], '"'), 
            "mat <- load_project_snpmat(current_project)", sep = "\n")
    
    size_line <- paste0("                       species_genome_size = ", 
                        project_df$species_size[project], 
                        ",")
    
    print_statement <- 
      print(paste('#',
                  project_df$projects[project], 
                  percent_samples[percent], 
                  '% of total samples'), sep = " ")
    percent_line <- 
      paste0('percent_samples <- ', percent_samples[percent])
    project_line <-
      paste0('current_project <- ', '"', project_df$projects[project], '"')
    species_line <- 
      paste0('current_species <- ', '"', project_df$species[project], '"')
    
    commands <- 
      paste(print_statement, 
            "source('../lib/multiallelic_lib.R')",
            "source('../lib/subsample_lib.R')", 
            "# Initialize data",  
            "num_rep <- 10",
            percent_line, 
            project_line, 
            species_line, 
            "data <- generate_nrow_df(num_rep)",
            project_specific_code, 
            "num_to_select <- round(ncol(mat) * percent_samples/100, 0)", 
            "for (i in 1:num_rep) {", 
            "  sub_mat <- subsample_snpmat(mat, num_to_select, seed = i)", 
            "  sub_mat_vecs <- keep_only_variant_sites(sub_mat)", 
            "  data <- fill_in_data(data_mat = data,",  
            "                       mat = sub_mat_vecs$mat,",  
            "                       row_num = i,", 
            size_line, 
            "                       bi_vec = sub_mat_vecs$biallelic_vec,", 
            "                       tri_vec = sub_mat_vecs$triallelic_vec,", 
            "                       quad_vec = sub_mat_vecs$quadallelic_vec,", 
            "                       multi_vec = sub_mat_vecs$multiallelic_vec)",
            "data$Project[i] <- current_project", 
            "data$Species[i] <- current_species", 
            "write_tsv(data,", 
            "          path = paste0('../data/multiallelic_variant_summary_',", 
            "          current_project,",
            "          '_',", 
            "          percent_samples,", 
            "          '_percent',",
            "          '.tsv'))",
            "}", 
            sep = "\n")
    write(commands, 
          file = fname,
          append = FALSE)
    
    command <- paste("module load R/3.6.1 Rtidyverse/3.6.1", 
                     paste0("Rscript ", fname), 
                     sep = "\n")
    
    sbat_name <- paste0("../data/run_", 
                       project_df$projects[project], 
                       "_subsampled_to_", 
                       percent_samples[percent], 
                       "_percent", ".sbat")

    writeLines(c("#!/bin/sh",
                 paste0("#SBATCH --job-name=subsample_", 
                        project_df$projects[project], 
                        percent_samples[percent]),
                 "#SBATCH --mail-user=katiephd@umich.edu",  
                 "#SBATCH --mail-type=BEGIN,END,NONE,FAIL,REQUEUE",
                 "#SBATCH --export=ALL",
                 "#SBATCH --partition=standard",
                 "#SBATCH --account=esnitkin",
                 "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=30g --time=10:00",
                 "cd $SLURM_SUBMIT_DIR",
                 "echo $SLURM_SUBMIT_DIR",
                 command), # module load plus Rscript name
               sbat_name,
               sep = "\n")
  }
}