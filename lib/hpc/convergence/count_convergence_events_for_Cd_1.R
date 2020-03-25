source("/nfs/esnitkin/Project_Prep_GWAS_Variants/Analysis/Prep_GWAS_SNPs/2019-10-12_count_convergence/lib/count_convergent_events_lib.R")
load('/nfs/esnitkin/Project_Prep_GWAS_Variants/Analysis/Prep_GWAS_SNPs/2019-10-11_prewas_snpmats/data//prewas_output_cdif_st1.RData')
snpmat <- prewas_results$bin_mat
tree <- read.tree('/nfs/esnitkin/Project_Cdiff/Analysis/regional_cdiff_transmission/2019-08-07_format_data_for_gwas/data/trees/2019_07_14_10_25_16_Cdiff_ST1_R20291_whole_genome.treefile')
tree <- ape::root(tree, "cdiff_630")
snpmat <- snpmat[, colnames(snpmat) %in% tree$tip.label, drop = FALSE]
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% colnames(snpmat)])
tree$edge.length[tree$edge.length == 0] <- 0.0000000001
snpmat <- snpmat[, match(tree$tip.label, colnames(snpmat)), drop = FALSE]
transition_mat <- count_convergent_events(snpmat, tree)
multiallelic_log <- identify_multiallelic_sites(snpmat)
biallelic_trans_mat <- make_biallelic_trans_mat(transition_mat, multiallelic_log)
multiallelic_trans_mat <- make_multiallelic_trans_mat(transition_mat, multiallelic_log)
file_name <- "../figures/convergence_histogram_cdif_st1.pdf"
plot_num_convergence_histogram(biallelic_trans_mat, multiallelic_trans_mat, file_name)
save_results('cdif_st1', transition_mat, biallelic_trans_mat, multiallelic_trans_mat)
