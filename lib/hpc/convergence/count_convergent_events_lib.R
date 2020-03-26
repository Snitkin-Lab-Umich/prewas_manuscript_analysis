library(magrittr)
library(ape)
library(phytools)

count_convergent_events <- function(snpmat, tree) {
  transition_mat <- matrix(NA, nrow = nrow(snpmat), ncol = 1)
  row.names(transition_mat) <- row.names(snpmat)
  colnames(transition_mat) <- "times_snp_appears"
  set.seed(1)
  for (i in 1:nrow(snpmat)) {
    if (sum(snpmat[i, , drop = TRUE]) != 0 & sum(snpmat[i, , drop = TRUE]) != ncol(snpmat)) { # Skip rows with no variants
      current_recon <- ape::ace(x = snpmat[i, ],
                                phy = tree, 
                                type = "discrete",
                                method = "ML", 
                                marginal = FALSE)
      # Extract the mostly likely character state using which.max
      current_recon <-
        as.numeric(colnames(current_recon$lik.anc)[apply(current_recon$lik.anc,
                                                         1,
                                                         which.max)])
      
      transition_direction <- parent_node <- child_node <- 
        integer(ape::Nedge(tree))
      older <- 1 # older node is 1st column in tr$edge
      younger <- 2 # younger node is 2nd column in tr$edge
      parent_0_child_1 <- 1
      
      for (j in 1:Nedge(tree)) {
        parent_node[j] <- current_recon[tree$edge[j, older] - ape::Ntip(tree)]
        if (tree$edge[j, younger] <= ape::Ntip(tree)) {
          # child is a tip
          child_node[j] <- unlist(snpmat[i, ][tree$edge[j, younger]])
        } else {
          # child is internal nodes
          child_node[j] <- 
            unlist(current_recon[tree$edge[j, younger] - ape::Ntip(tree)])
        }
        if (parent_node[j] < child_node[j]) {
          transition_direction[j] <- parent_0_child_1
        }
      }
      
      transition_mat[i, 1] <- sum(transition_direction)
    }
  }
  return(transition_mat)
}

identify_multiallelic_sites <- function(snpmat) {
  multiallelic_log <- 
    row.names(snpmat) %>% 
    gsub("[|].*", "", .) %>% 
    grepl(pattern = "[.]", .)
  return(multiallelic_log)
}


make_biallelic_trans_mat <- function(transition_mat, multiallelic_log){
  biallelic_trans_mat <- transition_mat[!multiallelic_log, , drop = FALSE]
  biallelic_trans_mat <- biallelic_trans_mat[!is.na(biallelic_trans_mat), , drop = FALSE]
  return(biallelic_trans_mat)
}

make_multiallelic_trans_mat <- function(transition_mat, multiallelic_log){
  multiallelic_trans_mat <- transition_mat[multiallelic_log, , drop = FALSE]
  multiallelic_trans_mat <- multiallelic_trans_mat[!is.na(multiallelic_trans_mat), , drop = FALSE]
  return(multiallelic_trans_mat)
}

plot_num_convergence_histogram <- function(biallelic_trans_mat,
                                           multiallelic_trans_mat, 
                                           file_name){
  pdf(file_name)
  bi_color <- rgb(1, 0, 0, 0.25)
  multi_color <- rgb(0, 0, 1, 0.25)
  hist(biallelic_trans_mat[,1], 
       col = bi_color, 
       xlab = "Number of SNP convergence events on tree", 
       ylab = "Count", 
       main = "")
  hist(multiallelic_trans_mat[,1], col = multi_color, add = TRUE)
  legend("topright", c("Biallelic", "Multiallelic"), 
         fill = c(bi_color, multi_color))
  dev.off()
}

save_results <- function(project_name, 
                         transition_mat, 
                         biallelic_trans_mat, 
                         multiallelic_trans_mat) {
  write.csv(transition_mat, 
            file = paste0("../data/", project_name, "_transition_matrix.csv"))
  write.csv(biallelic_trans_mat, 
            file = paste0("../data/",
                          project_name, 
                          "_biallelic_transition_matrix.csv"))
  write.csv(multiallelic_trans_mat, 
            file = paste0("../data/",
                          project_name,
                          "_multiallelic_transition_matrix.csv"))
}