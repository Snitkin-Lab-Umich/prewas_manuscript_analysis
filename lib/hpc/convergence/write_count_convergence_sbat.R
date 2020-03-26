library(tidyverse)
prewas_inputs <- read_tsv("../../2020-03-13_prewas_resource_usage/data/project_names_and_colors.tsv",col_names = TRUE)                         

prewas_results_dir <- "../../2020-03-13_prewas_resource_usage/data/"
tree_dir <- "../../2019-11-20_clean_data/data/"
for (i in 1:nrow(prewas_inputs)) {
  fname <- paste0("../lib/count_convergence_events_for_",
                  prewas_inputs$Project[i], 
                  ".R")
  temp_prewas_output <- paste0(prewas_results_dir,
                               prewas_inputs$Project[i], 
                               "_prewas_output_with_ancestral_reconstruction_given_tree_10_cores.Rdata")
  
  if (file.exists(temp_prewas_output)) {
    print(temp_prewas_output)
    
    tree_name <- paste0(tree_dir,
                        "/", 
                        prewas_inputs$Project[i], 
                        "_.tree")
    
    commands <- paste('source("../lib/count_convergent_events_lib.R")', 
                      paste0('load(', "'", temp_prewas_output, "'", ')'),
                      'snpmat <- output$bin_mat',
                      paste0('tree <- read.tree(', "'", tree_name, "')"),
                      'tree <- midpoint.root(tree)', 
                      'snpmat <- snpmat[, colnames(snpmat) %in% tree$tip.label, drop = FALSE]', 
                      'tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% colnames(snpmat)])',
                      'tree$edge.length[tree$edge.length == 0] <- 0.0000000001', 
                      'snpmat <- snpmat[, match(tree$tip.label, colnames(snpmat)), drop = FALSE]', 
                      'print("SNPMAT dimensions and Ntip")',
                      'print(dim(snpmat))',
                      'print(ape::Ntip(tree))',
                      'transition_mat <- count_convergent_events(snpmat, tree)',
                      'multiallelic_log <- identify_multiallelic_sites(snpmat)',
                      'biallelic_trans_mat <- make_biallelic_trans_mat(transition_mat, multiallelic_log)',
                      'multiallelic_trans_mat <- make_multiallelic_trans_mat(transition_mat, multiallelic_log)', 
                      paste0('file_name <- "../figures/convergence_histogram_', prewas_inputs$Project[i], '.pdf"'),
                      'plot_num_convergence_histogram(biallelic_trans_mat, multiallelic_trans_mat, file_name)',
                      paste0('save_results(', "'", prewas_inputs$Project[i], "',", ' transition_mat, biallelic_trans_mat, multiallelic_trans_mat)'),
                      sep = "\n")
    write(commands, 
          file = fname,
          append = FALSE)
    
    command <- paste("module load R/3.6.1 Rtidyverse/3.6.1", 
                     paste0("Rscript ", fname), 
                     sep = "\n")
    
    sbat_name <- paste0("../data/count_convergence_events_", 
                        prewas_inputs$Project[i], 
                        ".sbat")
    
    writeLines(c("#!/bin/sh",
                 paste0("#SBATCH --job-name=", prewas_inputs$Project[i], "_convergence"),
                 paste0("#SBATCH --output=", prewas_inputs$Project[i], "_convergence.out"),
                 "#SBATCH --mail-user=katiephd@umich.edu",  
                 "#SBATCH --mail-type=END,FAIL",
                 "#SBATCH --export=ALL",
                 "#SBATCH --partition=standard",
                 "#SBATCH --account=esnitkin1",
                 "#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=10g --time=20:00:00",
                 "cd $SLURM_SUBMIT_DIR",
                 "echo $SLURM_SUBMIT_DIR",
                 "echo $SLURM_JOB_ID",
                 command), # module load plus Rscript name
               sbat_name,
               sep = "\n") 
  }
}
