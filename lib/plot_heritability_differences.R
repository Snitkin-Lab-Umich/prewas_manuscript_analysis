# plot heritability estimates
library(tidyverse)

cols = read_delim('../data/project_names_and_colors.tsv', delim='\t')
her_ests = read.csv('data/her_ests.csv', header = F, stringsAsFactors = F)
colors = cols$hex_color
names(colors) = cols$Project
proj_name = cols$Dataset
names(proj_name) = cols$Project

names(her_ests) = c('datname','phen_type','pw_type','her_est')

her_ests$dataset = gsub('_prewas.*','',her_ests$datname)
her_ests$prewas_run = gsub('.*prewas_output_','',her_ests$datname)
her_ests$prewas_run[her_ests$prewas_run == 'no_ancestral_reconstruction_no_tree_given'] = 'Major Allele'
her_ests$prewas_run[her_ests$prewas_run == 'with_ancestral_reconstruction_no_tree_given'] = 'AR + NJ tree'
her_ests$prewas_run[her_ests$prewas_run == 'with_ancestral_reconstruction_given_tree'] = 'AR + ML tree'
her_ests$phen_type[her_ests$phen_type == 'bm'] = 'BM'
her_ests$phen_type[her_ests$phen_type == 'wn'] = 'WN'

he = her_ests[her_ests$pw_type == 'multi',]
he$diff = (her_ests$her_est[her_ests$pw_type == 'multi'] - her_ests$her_est[her_ests$pw_type == 'no multi'])*100

pdf('figures/Figure_S3_heritability_estimate_differences.pdf',width = 6,height = 4)
ggplot(he, aes(x=phen_type,y=diff)) + geom_boxplot() + geom_jitter(aes(col=dataset),size=3) + facet_grid(cols = vars(prewas_run)) + 
  xlab('Phenotype Evolutionary Model') + ylab('Heritability with multiallelic -\nHeritability without multiallelic (%)') + 
  scale_color_manual(name = 'Dataset',values = colors, labels = proj_name) + 
  theme_bw() +
  theme(text = element_text(size=10))
dev.off()
