# Run each of the datasets for which we have a multiVCF file through prewas.
# The outfile for each run will record the jobID. 
library(prewas)

run_prewas_on_cluster <- function(project_name, 
                                  ancestral_reconstruction_log, 
                                  given_tree_log, 
                                  num_cores) {
  AR_phrase <- "no_ancestral_reconstruction"
  if (ancestral_reconstruction_log) {
    AR_phrase <- "with_ancestral_reconstruction"
  }
  
  tree_phrase <- "no_tree_given"
  if (given_tree_log) {
    tree_phrase <- "given_tree"
  }
  
  vcf_path <- paste0("../data/cleaned_", project_name, ".vcf")
  
  tree_path <- NULL
  if (given_tree_log) {
    tree_path <- paste0("../../2019-11-20_clean_data/data/", 
                        project_name, 
                        "_.tree")
  }

  out_name <- paste0("../data/",
                     project_name, 
                     "_prewas_output_",
                     AR_phrase, 
                     "_", 
                     tree_phrase, 
                     "_", 
                     num_cores, 
                     "_cores.Rdata")
  output <- prewas::prewas(dna = vcf_path, 
                 tree = tree_path,
                 outgroup = NULL, 
                 gff = NULL, 
                 anc = ancestral_reconstruction_log)
  save(output, file = paste0("../data/", out_name))
}

args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the sbatch script
project_name <- args[1]
AR_log <- args[2]
tree_log <- args[3]
cores <- args[4]

run_prewas_on_cluster(project_name, AR_log, tree_log, cores)