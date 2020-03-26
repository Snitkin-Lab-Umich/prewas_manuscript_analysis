# This script generates the .sbat files needed to calculate the summary of 
# multiallelic sites for each project. 

source('../lib/multiallelic_lib.R')

for (project in 1:nrow(project_df)) {
  fname <- paste0("../lib/generate_multiallelic_site_summary_",
                  project_df$projects[project], 
                  ".R")
  size_line <- paste0("                     species_genome_size = ", 
                      project_df$species_size[project], 
                      ",")
  
  
  print_statement <- print(paste0('# ', project_df$projects[project]))
  project_line <- paste0('current_project <- ', '"', project_df$projects[project], '"')
  species_line <- paste0('current_species <- ', '"', project_df$species[project], '"')
  
  commands <- paste(print_statement, 
                    "source('../lib/multiallelic_lib.R')",
                    "# Initialize data",  
                    "num_rep <- 10",
                    project_line, 
                    species_line, 
                    "data <- generate_single_df()",
                    "mat <- load_project_snpmat(current_project)", 
                    "mat_vecs <- keep_only_variant_sites(mat)", 
                    "data <- fill_in_data(data_mat = data,",  
                    "                     mat = mat_vecs$mat,",  
                    "                     row_num = one_row,", 
                    size_line, 
                    "                     bi_vec = mat_vecs$biallelic_vec,", 
                    "                     tri_vec = mat_vecs$triallelic_vec,", 
                    "                     quad_vec = mat_vecs$quadallelic_vec,", 
                    "                     multi_vec = mat_vecs$multiallelic_vec)",
                    "data$Project[one_row] <- current_project", 
                    "data$Species[one_row] <- current_species", 
                    "write_tsv(data,", 
                    "          path = paste0('../data/multiallelic_variant_summary_',", 
                    "          current_project,",
                    "          '.tsv'))",
                    sep = "\n")
  write(commands, 
        file = fname,
        append = FALSE)
  
  command <- paste("module load R/3.6.1 Rtidyverse/3.6.1", 
                   paste0("Rscript ", fname), 
                   sep = "\n")
  
  sbat_name <- paste0("../data/summarize_multiallelic_in_", 
                     project_df$projects[project], 
                     ".sbat")
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=multi_", project_df$projects[project]),
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
