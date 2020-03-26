subsample_snpmat <- function(snpmat, num_samples, seed){
  set.seed(seed)
  selected_cols <- sample(x = 1:ncol(snpmat),
                          size = num_samples, 
                          replace = FALSE)
  subsampled_mat <- snpmat[, selected_cols, drop = FALSE]
  return(subsampled_mat)
}
