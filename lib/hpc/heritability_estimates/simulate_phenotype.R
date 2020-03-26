# simulate continuous norma phenotypes

# we'll just need to add a while loop to check for normally distrbuted data
make_continuous_phenotypes <- function(tree_list, num_pheno){
  num_trees <- length(tree_list)
  pheno_mat_list_BM <- pheno_mat_list_WN <- rep(list(NULL), num_trees)
  set.seed(1)
  for (i in 1:num_trees) {
    #print(paste('tree',i))
    #print(tree_list[[i]])
    pheno_mat_list_BM[[i]] <- pheno_mat_list_WN[[i]] <- 
      matrix(NA, nrow = ape::Ntip(tree_list[[i]]), ncol = num_pheno)
    colnames(pheno_mat_list_BM[[i]]) <- paste0("BM_pheno_", 1:num_pheno)
    colnames(pheno_mat_list_WN[[i]]) <- paste0("WN_pheno_", 1:num_pheno)
    
    for (j in 1:num_pheno) {
      
      bm_not_from_normal = TRUE
      while(bm_not_from_normal){
        lamdba_not_close_to_1 <- TRUE
        while (lamdba_not_close_to_1) {
          continuous_BM_pheno <- phytools::fastBM(tree = tree_list[[i]])
          # Check that lambda is close to 1 for BM phenotype
          BM_lambda <- phytools::phylosig(tree = tree_list[[i]],
                                          x = continuous_BM_pheno,
                                          method = "lambda")
          lamdba_not_close_to_1 <- BM_lambda$lambda < 0.95 & BM_lambda$lambda > 1.05
        }
        #print(shapiro.test(continuous_BM_pheno)$p.value)
        bm_not_from_normal = shapiro.test(continuous_BM_pheno)$p.value < 0.05
      }
      
      wn_not_from_normal = TRUE
      while(wn_not_from_normal){
        # When lambda is close to zero, the phylogenetic signal is low (White Noise)
        lambda_has_high_signal <- TRUE
        while (lambda_has_high_signal) {
          jumbled_pheno <- sample(unname(continuous_BM_pheno), 
                                  size = ape::Ntip(tree_list[[i]]), 
                                  replace = FALSE)
          names(jumbled_pheno) <- tree_list[[i]]$tip.label
          jumbled_lambda <- phytools::phylosig(tree = tree_list[[i]],
                                               x = jumbled_pheno,
                                               method = "lambda")
          lambda_has_high_signal <- jumbled_lambda$lambda < -0.05 & jumbled_lambda$lambda > 0.05
        }
        cont_low_signal_pheno <- jumbled_pheno
        #print(shapiro.test(cont_low_signal_pheno)$p.value)
        wn_not_from_normal = shapiro.test(cont_low_signal_pheno)$p.value < 0.05
      }
      
      
      continuous_BM_pheno <- 
        matrix(c(names(continuous_BM_pheno), continuous_BM_pheno), ncol = 2)  
      row.names(pheno_mat_list_BM[[i]]) <- continuous_BM_pheno[, 1]
      continuous_BM_pheno <- continuous_BM_pheno[, 2, drop = FALSE]
      
      cont_low_signal_pheno <- matrix(c(names(cont_low_signal_pheno), cont_low_signal_pheno), ncol = 2) 
      row.names(pheno_mat_list_WN[[i]]) <- cont_low_signal_pheno[, 1]
      cont_low_signal_pheno <- cont_low_signal_pheno[, 2, drop = FALSE]
      
      pheno_mat_list_BM[[i]][, j] <- continuous_BM_pheno
      pheno_mat_list_WN[[i]][, j] <- cont_low_signal_pheno
    }
  }
  results <- list("cont_pheno_BM_mat_list" = pheno_mat_list_BM, 
                  "cont_pheno_WN_mat_list" = pheno_mat_list_WN)
  return(results)
}