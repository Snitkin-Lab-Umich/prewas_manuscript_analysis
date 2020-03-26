# get shared variants
library(tidyverse)

# function to get shared variants from binary variant matrix
get_shared_variants = function(mat){
  return(apply(mat,2,function(x) colSums(mat == 1 & x == 1,na.rm=T)))
}

# inputs and outputs
f = snakemake@input[[1]]
phen_path = snakemake@input[[2]]
outfile_anc_multi = snakemake@output[[1]]
outfile_maj_multi = snakemake@output[[2]]
outfile_ref_multi = snakemake@output[[3]]
outfile_anc_nomulti = snakemake@output[[4]]
outfile_maj_nomulti = snakemake@output[[5]]
outfile_ref_nomulti = snakemake@output[[6]]

# get dataset name
dataset = gsub('../2019-11-20_clean_data/data//|_parsed.RData','',f)

# load dataset
load(f)

# get phenotype data
phen = read.csv(phen_path,row.names = 1)
# get annotations
annots = parsed$bin$annots
# get matrix referenced to ancestral allele
mat_anc = parsed$bin$mat
# rename to match phenotype names
names(mat_anc) = gsub('^X','',names(mat_anc)) %>% gsub('\\.','-',.)
# make sample order match phenotype sample order
mat_anc = mat_anc[,rownames(phen)]
# remove duplicated genes b/c of overlapping
rm_ogs = !(duplicated(paste0(annots$split_rows_flag, annots$rows_with_overlapping_genes_log)) & annots$rows_with_overlapping_genes_log)
# no multiallelic
rm_multi = !(duplicated(paste0(annots$split_rows_flag, annots$rows_with_mult_var_allele_log)) & annots$rows_with_mult_var_allele_log)
# multiallelic mat
mat_anc_multi = mat_anc[rm_ogs,]
# no multiallelic mat
mat_anc_nomulti = mat_anc[rm_ogs & rm_multi,]

# make major allele mat
mat_maj = mat_anc
anc_maj_diff = which(annots$anc != annots$maj)
for(i in anc_maj_diff){
  mat_maj[i,] = ifelse(unlist(mat_maj[i,])==1,0,1)
}
mat_maj_multi = mat_maj[rm_ogs,]
mat_maj_nomulti = mat_maj[rm_ogs & rm_multi,]

# make reference genome allele mat
mat_ref = mat_anc
anc_ref_diff = which(annots$anc != annots$ref)
for(i in anc_ref_diff){
  mat_ref[i,] = ifelse(unlist(mat_ref[i,])==1,0,1)
}
mat_ref_multi = mat_ref[rm_ogs,]
mat_ref_nomulti = mat_ref[rm_ogs & rm_multi,]

# get shared variant matrices (i.e. kinship matrices)
shared_anc_multi = get_shared_variants(mat_anc_multi)
shared_anc_nomulti = get_shared_variants(mat_anc_nomulti)
shared_maj_multi = get_shared_variants(mat_maj_multi)
shared_maj_nomulti = get_shared_variants(mat_maj_nomulti)
shared_ref_multi = get_shared_variants(mat_ref_multi)
shared_ref_nomulti = get_shared_variants(mat_ref_nomulti)

# write kinship matrices
write.csv(shared_anc_multi,outfile_anc_multi,quote=F,row.names = T)
write.csv(shared_maj_multi,outfile_maj_multi,quote=F,row.names = T)
write.csv(shared_ref_multi,outfile_ref_multi,quote=F,row.names = T)
write.csv(shared_anc_nomulti,outfile_anc_nomulti,quote=F,row.names = T)
write.csv(shared_maj_nomulti,outfile_maj_nomulti,quote=F,row.names = T)
write.csv(shared_ref_nomulti,outfile_ref_nomulti,quote=F,row.names = T)

# clean up
rm(parsed)

