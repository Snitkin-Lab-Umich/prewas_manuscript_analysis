# To be run from data/
library(tidyverse)
project_names <- read_tsv("../data/project_names_and_colors.tsv", 
                          col_names = TRUE)

# Projects not included in final prewas draft
project_names <- 
  project_names[!project_names$Project %in% 
                  c("Cd_1", "Cd_2", "Kp_2", "Ec_1", "Ece_1", "Lp_1", "Sa_3"), 
                , drop = FALSE]

write_prewas_sbat <- function(project_id, 
                              ancestral_reconstruction_logical,
                              given_tree_log, 
                              num_cores) {

  AR_state <- "major"
  if (ancestral_reconstruction_logical) {
    AR_state <- "w_AR"
  }
  tree_state <- "no_tr_given"
  if (given_tree_log) {
    tree_state <- "w_tr"
  }
  job_name <- paste0(project_id, "_", AR_state, "_", tree_state, "_", num_cores)
  out_name <- paste0(job_name, ".out")
  writeLines(c("#!/bin/sh",
               paste0("#SBATCH --job-name=", job_name),
               paste0("#SBATCH --output=", out_name),
               "#SBATCH --mail-user=katiephd@umich.edu",  
               "#SBATCH --mail-type=END",
               "#SBATCH --export=ALL",
               "#SBATCH --partition=standard",
               "#SBATCH --account=esnitkin1",
               paste0("#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=", 
                      num_cores, 
                      " --mem=40G --time=100:00:00"),
               "cd $SLURM_SUBMIT_DIR",
               "echo $SLURM_SUBMIT_DIR",
               "echo $SLURM_JOB_ID",
               paste(
                 "Rscript ../lib/run_prewas_on_multivcf_with_and_without_AR.R ",
                 project_id, 
                 " ", 
                 ancestral_reconstruction_logical, 
                 " ", 
                 given_tree_log, 
                 " ", 
                 num_cores)),
             paste0(getwd(), "/", job_name, ".sbat"),
             sep = "\n")
}

for (i in 1:nrow(project_names)) {
  prj_name <- project_names$Project[i] 
  vcf_file_name <- paste0("../data/", prj_name, ".vcf")
  if (file.exists(vcf_file_name)) {
    # AR + Given Tree + 1 core
    write_prewas_sbat(prj_name, TRUE, TRUE, 1)
    # AR + Given Tree + 4 core
    write_prewas_sbat(prj_name, TRUE, TRUE, 4)
    # AR + Given Tree + 10 core
    write_prewas_sbat(prj_name, TRUE, TRUE, 10)
    # AR + No Tree Provided + 1 core
    write_prewas_sbat(prj_name, TRUE, FALSE, 1)
    # AR + No Tree Provided + 4 core
    write_prewas_sbat(prj_name, TRUE, FALSE, 4)
    # AR + No Tree Provided + 10 core
    write_prewas_sbat(prj_name, TRUE, FALSE, 10)
    # Major + No Tree Provided + 1 core
    write_prewas_sbat(prj_name, FALSE, FALSE, 1)
  }
}
