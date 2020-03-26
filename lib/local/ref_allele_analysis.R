# This script finds the number and fraction of mismatches between genome reference, ancestral, and major allele variant matrices
# Input: 
#    data/hpc/snp_matrices_and_parsed_data/*annot*: Parsed annotations for each dataset
# Output: 
#    data/local/ref_allele/frac_mismatch.tsv: file with fraction of mismatches for each combination of reference alleles for each dataset
#    data/local/ref_allele/num_mismatch.tsv: file with number of mismatches for each combination of reference alleles for each dataset
#    data/local/ref_alleleanc_probs.tsv: file with information about ancestral reconstruction probabilities for each dataset 

library(ggplot2)
library(reshape2)

# get list of annotation files
files = list.files('data/hpc/snp_matrices_and_parsed_data/',pattern = 'annot',full.names = T)

# initialize fraction matrix
frac_mismatch = matrix(NA,nrow=length(files),ncol=3)
colnames(frac_mismatch) = c('ref_maj','ref_anc','maj_anc')
rownames(frac_mismatch) = gsub('.*//|_annot_split_no_overlap_no_ig.tsv|.gz','',files)

# initialize number matrix
num_mismatch = matrix(NA,nrow=length(files),ncol=6)
colnames(num_mismatch) = c('ref_maj','ref_anc','maj_anc','ref_ref','maj_maj','anc_anc')
rownames(num_mismatch) = gsub('.*//|_annot_split_no_overlap_no_ig.tsv|.gz','',files)

# initialize ancestral probabilities matrix
anc_probs = data.frame(dat=NA,anc_prob=NA)

# for each dataset
for(f in files){
  pref = gsub('.*//|_annot_split_no_overlap_no_ig.tsv|.gz','',f)
  print(pref)
  # get annotations
  annots = read.delim(f,stringsAsFactors = F)
  # get ancestral probabilities
  anc_prob = as.numeric(annots$anc_prob)
  # keep annotations with high ancestral state probability
  annots = annots[anc_prob >= 0.875,]
  # get different allele types
  ref_allele = annots$ref
  maj_allele = annots$maj
  anc_allele = annots$anc
  
  alleles = data.frame(reference=ref_allele,major=maj_allele,ancestral=anc_allele)
  rn_alleles = rownames(alleles)
  rownames(alleles) = 1:nrow(alleles)
  
  dat = data.frame(dat=rep(pref,length(anc_prob)),anc_prob)
  anc_probs = rbind(anc_probs,dat)
  
  alleles = alleles[alleles$major %in% c('A','C','G','T'),]
  
  # calculate mismatches (fractions)
  frac_mismatch[pref,'ref_maj'] = sum(maj_allele != ref_allele)/nrow(alleles)
  frac_mismatch[pref,'ref_anc'] = sum(anc_allele != ref_allele)/nrow(alleles)
  frac_mismatch[pref,'maj_anc'] = sum(anc_allele != maj_allele)/nrow(alleles)
  
  ## calculate mismatches (counts)
  num_mismatch[pref,'ref_maj'] = sum(maj_allele != ref_allele)
  num_mismatch[pref,'ref_anc'] = sum(anc_allele != ref_allele)
  num_mismatch[pref,'maj_anc'] = sum(anc_allele != maj_allele)
  num_mismatch[pref,'ref_ref'] = sum(ref_allele != ref_allele)
  num_mismatch[pref,'maj_maj'] = sum(maj_allele != maj_allele)
  num_mismatch[pref,'anc_anc'] = sum(anc_allele != anc_allele)
  
  rm(annots)
  
}

# write files
write.table(frac_mismatch,file = 'data/local/ref_allele/frac_mismatch.tsv',sep = '\t')
write.table(num_mismatch,file = 'data/local/ref_allele/num_mismatch.tsv',sep = '\t')

anc_probs = data.frame(anc_probs[2:nrow(anc_probs),])
write.table(anc_probs,file = 'data/local/ref_allele/anc_probs.tsv',sep='\t')






