# SOURCE PARSER 
source('../../../../../Github/snitkitr/R/parse_snps.R')
library(ape)

# NAME OF DATASET 
dataset_name = 'Sa_2'

# READ IN ALLELE MAT PATH LOOKUP TABLE 
lookup = read.csv('../../data/project_paths.txt', header = TRUE, stringsAsFactors = FALSE)

# GET PATH TO ALLELE MAT FOR DATASET OF INTEREST 
abs_path = lookup$path[lookup$dataset == dataset_name]
rel_path = paste0('../../../../../', gsub('/nfs/esnitkin/', '', abs_path))

path_to_snpmat_code =
  paste0(rel_path, '/data_matrix/matrices/SNP_matrix_code.csv')
path_to_snpmat_allele = 
  paste0(rel_path, '/data_matrix/matrices/SNP_matrix_allele_new.csv')
path_to_tree = 
  list.files(path = paste0(rel_path, '/gubbins/iqtree_results/'),
             pattern = '.*_aln_w_alt_allele_unmapped.*.treefile',
             full.names = TRUE)
  
# READ IN ALLELE MAT 
allele_mat = read.table(path_to_snpmat_allele,
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t",
                        quote = "", 
                        row.names = 1,
                        flush = TRUE)

# READ IN CODE MAT
code_mat = read.table(path_to_snpmat_code,
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      sep = "\t",
                      quote = "", 
                      row.names = 1)
# READ IN TREE 
tree = read.tree(path_to_tree)

# CHANGE NAMES IN MATRICES AND TREE TO BE THE SAME
names(allele_mat) = gsub('^X|_R1_001.fastq.gz','',names(allele_mat))
names(code_mat) = gsub('^X|_R1_001.fastq.gz','',names(code_mat))
names(allele_mat) = gsub('\\.','-',names(allele_mat))
names(code_mat) = gsub('\\.','-',names(code_mat))
names(allele_mat) = gsub('_S.*|_L.*','',names(allele_mat))
names(code_mat) = gsub('_S.*|_L.*','',names(code_mat))
tree$tip.label = gsub('_$','',tree$tip.label)

# REMOVE BAD SAMPLES & SAMPLES WE'RE NOT ALLOWED TO USE FROM ALLELE MAT,
#   CODE MAT, TREE (AND REFERENCE GENOME FROM TREE)
to_drop = c("MRSA_USA_100_DAR3548")
tree_subset = drop.tip(tree, to_drop)
allele_mat = allele_mat[, !names(allele_mat) %in% to_drop]
code_mat = code_mat[, !names(code_mat) %in% to_drop]

# CHECK THAT THE SNPMAT AND TREE HAVE SAME LABELS 
if (length(setdiff(tree_subset$tip.label, colnames(allele_mat))) > 0 | 
    length(setdiff(colnames(allele_mat), tree_subset$tip.label)) > 0){
  stop('tree tip labels do not match snpmat column names')
}

# RUN PARSER 
# SPECIFY OUTGROUP IF ONE EXISTS 
parsed = parse_snps(code_mat, 
                    allele_mat, 
                    tree = tree_subset, 
                    og = NULL, 
                    remove_multi_annots = FALSE, 
                    return_binary_matrix = TRUE, 
                    ref_to_anc = T)
#parsed = parse_snps(code_mat, allele_mat, tree = NULL, og = NULL, 
# remove_multi_annots = FALSE, return_binary_matrix = FALSE, ref_to_anc = F)


#### AFTER PARSER ####
# REMOVE INTERGENIC SITES 
# COLLAPSE SPLIT ALLELE MAT 

intergenic_unsplit = 
  parsed$allele$annots$intergenic[!duplicated(parsed$allele$annots$split_rows_flag)]
allele_mat_unsplit = 
  parsed$allele$mat[!duplicated(parsed$allele$annots$split_rows_flag), ]
raw_rownames_unsplit = 
  parsed$allele$annots$raw_rownames[!duplicated(parsed$allele$annots$split_rows_flag)]
  
allele_mat_unsplit_no_ig = allele_mat_unsplit[!intergenic_unsplit,] 
row.names(allele_mat_unsplit_no_ig) = raw_rownames_unsplit[!intergenic_unsplit]

# GET ANNOTATIONS (NO OVERLAPPING GENES, NO INTERGENIC REGIONS)

annots_split_no_ol_no_ig <- parsed$allele$annots[!parsed$allele$annots$rows_with_overlapping_genes_log &
                                                   !parsed$allele$annots$intergenic,]

# SAVE CLEANED MATRIX, ANNOTATIONS, AND TREE
prefix <- paste0("../data/", dataset_name, "_")

write.table(x = allele_mat_unsplit_no_ig, 
            file = paste0(prefix, "allele_mat_unsplit_no_ig.tsv"), 
            append = FALSE, 
            sep = "\t")

write.table(x = annots_split_no_ol_no_ig, 
            file = paste0(prefix, "annot_split_no_overlap_no_ig.tsv"), 
            append = FALSE, 
            sep = "\t")

write.tree(phy = tree_subset, file = paste0(prefix, ".tree"))

save(parsed, file = paste0(prefix, "parsed.RData"))



