library(ggplot2)
library(reshape2)

files = list.files('../../2019-11-20_clean_data/data/',pattern = 'annot',full.names = T)

frac_mismatch = matrix(NA,nrow=length(files),ncol=3)
colnames(frac_mismatch) = c('ref_maj','ref_anc','maj_anc')
rownames(frac_mismatch) = gsub('.*//|_annot_split_no_overlap_no_ig.tsv','',files)

num_mismatch = matrix(NA,nrow=length(files),ncol=6)
colnames(num_mismatch) = c('ref_maj','ref_anc','maj_anc','ref_ref','maj_maj','anc_anc')
rownames(num_mismatch) = gsub('.*//|_annot_split_no_overlap_no_ig.tsv','',files)

anc_probs = data.frame(dat=NA,anc_prob=NA)

for(f in files){
  pref = gsub('.*//|_annot_split_no_overlap_no_ig.tsv','',f)
  print(pref)
  annots = read.delim(f,stringsAsFactors = F)
  anc_prob = as.numeric(annots$anc_prob)
  annots = annots[anc_prob >= 0.875,]
  #anc_prob_high_conf = as.numeric(annots$anc_prob)
  ref_allele = annots$ref
  maj_allele = annots$maj
  anc_allele = annots$anc
  
  alleles = data.frame(reference=ref_allele,major=maj_allele,ancestral=anc_allele)
  rn_alleles = rownames(alleles)
  rownames(alleles) = 1:nrow(alleles)
  
  dat = data.frame(dat=rep(pref,length(anc_prob)),anc_prob)
  anc_probs = rbind(anc_probs,dat)
  
  alleles = alleles[alleles$major %in% c('A','C','G','T'),]
  
  frac_mismatch[pref,'ref_maj'] = sum(maj_allele != ref_allele)/nrow(alleles)
  frac_mismatch[pref,'ref_anc'] = sum(anc_allele != ref_allele)/nrow(alleles)
  frac_mismatch[pref,'maj_anc'] = sum(anc_allele != maj_allele)/nrow(alleles)
  
  num_mismatch[pref,'ref_maj'] = sum(maj_allele != ref_allele)
  num_mismatch[pref,'ref_anc'] = sum(anc_allele != ref_allele)
  num_mismatch[pref,'maj_anc'] = sum(anc_allele != maj_allele)
  num_mismatch[pref,'ref_ref'] = sum(ref_allele != ref_allele)
  num_mismatch[pref,'maj_maj'] = sum(maj_allele != maj_allele)
  num_mismatch[pref,'anc_anc'] = sum(anc_allele != anc_allele)
  
  rm(annots)
  
}

write.table(frac_mismatch,file = '../data/frac_mismatch.tsv',sep = '\t')
write.table(num_mismatch,file = '../data/num_mismatch.tsv',sep = '\t')

anc_probs = data.frame(anc_probs[2:nrow(anc_probs),])
write.table(anc_probs,file = '../data/anc_probs.tsv',sep='\t')






