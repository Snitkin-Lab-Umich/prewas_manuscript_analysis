# compare heritability based on snps pre- and post-prewas

# so what we need to do is:
# - simulate phenotype data for each tree in dataset
# - get binary snp matrix without multiallelic sites
# - get post-prewas snp matrix with multiallelic sites
# - find pairwise snp distances
# - run limix

# load libraries
library(ape)
library(tidyverse)
source('lib/simulate_phenotype.R')

# rdata file
#args[[1]] = rdas[15]
rda = snakemake@input[[1]]
out_phen = snakemake@output[[1]]
out_multi = snakemake@output[[2]]
out_no_multi = snakemake@output[[3]]

# load prewas output
load(rda)

# get project
basename = gsub('../2020-03-13_prewas_resource_usage/data/|_1_cores.Rdata','',rda)
project_name = gsub('_prewas_output.*','',basename)

# get tree
tree = output$tree
if(is.null(tree)){
  tree = read.tree(paste0('../2019-11-20_clean_data/data/',project_name,'_.tree'))
}

# get phenotypes
phenotypes = make_continuous_phenotypes(list(tree),1)
phenotypes = cbind(phenotypes$cont_pheno_BM_mat_list[[1]],phenotypes$cont_pheno_WN_mat_list[[1]])

# binmat with multiallelic sites
mat_multi = output$bin_mat

# binmat with multiallelic sites
multi = gsub('\\.1','',rownames(mat_multi)[grep('\\.1',rownames(mat_multi))])
no_multi = !sapply(rownames(mat_multi),function(x) gsub('\\.1','',x) %in% multi)
mat_no_multi = mat_multi[no_multi,]

#minimat = matrix(sample(c(0,1),16,replace=T),ncol=4)

# find pariwise snp distances
pw_dists_multi = apply(mat_multi, 2, function(x){
  colSums(mat_multi != x)
})
pw_dists_no_multi = apply(mat_no_multi, 2, function(x){
  colSums(mat_no_multi != x)
})

# save output in format for limix
# phenotypes
phenotypes = phenotypes[rownames(pw_dists_multi),]
write.csv(phenotypes,out_phen,quote=F)
# pairwise snp distances (correlation matrix)
write.csv(pw_dists_multi,out_multi,quote=F,row.names = F)
write.csv(pw_dists_no_multi,out_no_multi,quote=F,row.names = F)




# # phenotypes
# write.csv(phenotypes,paste0('data/',basename,'_phenotype.csv'),quote=F)
# # pairwise snp distances (correlation matrix)
# write.csv(pw_dists_multi,paste0('data/',basename,'_pw_dists_multi.csv'),quote=F,row.names = F)
# write.csv(pw_dists_no_multi,paste0('data/',basename,'_pw_dists_no_multi.csv'),quote=F,row.names = F)
# 
# 
