source("../lib/count_convergent_events_lib.R")
load('../../2020-03-13_prewas_resource_usage/data/Sa_1_prewas_output_with_ancestral_reconstruction_given_tree_10_cores.Rdata')
snpmat <- output$bin_mat
tree <- read.tree('../../2019-11-20_clean_data/data//Sa_1_.tree')
tree <- midpoint.root(tree)
snpmat <- snpmat[, colnames(snpmat) %in% tree$tip.label, drop = FALSE]
tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% colnames(snpmat)])
tree$edge.length[tree$edge.length == 0] <- 0.0000000001
snpmat <- snpmat[, match(tree$tip.label, colnames(snpmat)), drop = FALSE]
print("SNPMAT dimensions and Ntip")
print(dim(snpmat))
print(ape::Ntip(tree))
transition_mat <- count_convergent_events(snpmat, tree)
multiallelic_log <- identify_multiallelic_sites(snpmat)
biallelic_trans_mat <- make_biallelic_trans_mat(transition_mat, multiallelic_log)
multiallelic_trans_mat <- make_multiallelic_trans_mat(transition_mat, multiallelic_log)
file_name <- "../figures/convergence_histogram_Sa_1.pdf"
plot_num_convergence_histogram(biallelic_trans_mat, multiallelic_trans_mat, file_name)
save_results('Sa_1', transition_mat, biallelic_trans_mat, multiallelic_trans_mat)
