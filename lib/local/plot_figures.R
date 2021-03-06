# Libraries --------------------------------------------------------------------
library(magrittr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
source("lib/local/multiallelic_lib.R")

# Color palettes and project names ---------------------------------------------
color_key <- read_tsv("data/key/plot_colors.tsv")
project_key <- read_tsv("data/key/project_names_and_colors.tsv")
palette <-  project_key$hex_color
names(palette) <- project_key$Dataset

# Save function & Figure defaults ----------------------------------------------
default_width <- 8
default_height <- 6
default_font_size <- 20
rsquared_font_size <- 7
default_dot_size <- 8
default_line_size <- 2

save_as_pdf_eps_png <- function(file_prefix, width_size, height_size){
  ggsave(filename = paste0("figures/", file_prefix, ".pdf"),
         width = width_size,
         height = height_size,
         units = "in",
         device = "pdf")
  ggsave(filename = paste0("figures/", file_prefix, ".eps"),
         width = width_size,
         height = height_size,
         units = "in",
         device = "eps")
  ggsave(filename = paste0("figures/", file_prefix, ".png"),
         width = width_size,
         height = height_size,
         units = "in",
         device = "png")
}

# Update Dataset name for final publication
rename_datasets <- function(df) {
  df$Dataset[df$Dataset == "C. difficile #3"] <- "C. difficile no. 1"
  df$Dataset[df$Dataset == "C. difficile #4"] <- "C. difficile no. 2"
  df$Dataset[df$Dataset == "E. faecium #1"] <- "E. faecium"
  df$Dataset[df$Dataset == "E. faecalis #1" ] <- "E. faecalis" 
  df$Dataset[df$Dataset == "K. pneumoniae #1"] <- "K. pneumoniae"
  df$Dataset[df$Dataset == "L. crispatus #1" ] <- "L. crispatus" 
  df$Dataset[df$Dataset == "S. maltophilia #1"] <- "S. maltophilia"
  df$Dataset[df$Dataset == "S. aureus #1"] <- "S. aureus no. 1"
  df$Dataset[df$Dataset == "S. aureus #2"] <- "S. aureus no. 2"
  return(df)
}


# Read in Multiallelic Data ----------------------------------------------------
variant_summary <- read_tsv("data/local/multiallelic/multiallelic_summary.tsv")
snp_mat_summary <- read_tsv("data/local/multiallelic/multiallelic_summary_subsampled.tsv")
snp_mat_grouped <- read_tsv("data/local/multiallelic/multiallelic_summary_by_project.tsv")


functional_impact_long <- read_tsv("data/local/multiallelic/SNP_diversity_long.tsv")
functional_variant_summary <- read_tsv(file = "data/local/multiallelic/snpeff_multiallelic_summary.tsv")
tri_functional_summary <- read_tsv(file = "data/local/multiallelic/snpeff_triallelic_summary.tsv")
tri_functional_heatmap <- read.csv(file = "data/local/multiallelic/triallelic_heatmap.csv", row.names = 1)
convergence <- read_tsv(file = "data/local/convergence/convergence.tsv")

# Filter out datasets to exclude:
drop_projects <- c("Kp_2", "Cd_1", "Cd_2")
variant_summary <- variant_summary %>% filter(!Project %in% drop_projects)
snp_mat_summary <- snp_mat_summary %>% filter(!Project %in% drop_projects)
snp_mat_grouped <- snp_mat_grouped %>% filter(!Project %in% drop_projects)
functional_impact_long <- functional_impact_long %>%
  filter(!Project %in% drop_projects)
functional_variant_summary <- functional_variant_summary %>%
  filter(!Project %in% drop_projects)
tri_functional_summary <- tri_functional_summary %>%
  filter(!Project %in% drop_projects)
convergence <- convergence %>% filter(!Project %in% drop_projects)

# Figure 2: Multiallelic Plots -------------------------------------------------
## 2A: % Multiallelic while subsampling datasets --------------------------------
snp_mat_grouped <- rename_datasets(snp_mat_grouped)
snp_mat_grouped %>%
  ggplot(aes(x = Percent / 100,
             y = AvgNumMultiallelicSite / AvgNumVariantSite,
             color = Dataset)) +
  geom_line(size = default_line_size) +
  geom_point(size = 4) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent_format(accuracy = .1)) +
  scale_x_continuous(labels = scales::percent) +
  #scale_color_manual(values = palette) +
  ylab("All multiallelic sites / all variant sites") +
  xlab("Dataset subset") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  theme(text = element_text(size = default_font_size)) + 
  scale_color_manual(values = palette, 
                     labels = c(expression(paste(italic("C. difficile "), "no. 1")), 
                                expression(paste(italic("C. difficile "), "no. 2")), 
                                expression(italic("E. faecalis")), 
                                expression(italic("E. faecium")), 
                                expression(italic("K. pneumoniae")), 
                                expression(italic("L. crispatus")), 
                                expression(paste(italic("S. aureus "), "no. 1")), 
                                expression(paste(italic("S. aureus "), "no. 2")), 
                                expression(italic("S. maltophila")))) +
  theme(legend.text.align = 0)

save_as_pdf_eps_png("Figure_2A_multiallelic_vs_subsample",
                    default_width,
                    default_height)

## 2B: % Multiallelic vs Maximum Pairwise SNP distance -------------------------
multi_vs_MEAN_model <-
  lm(formula = 100 * NumMultiallelicSite/NumVariantSite ~ MeanSNPDist,
     data = variant_summary)
r2 <- round(summary(multi_vs_MEAN_model)$r.squared, 2)

r2_mean <- r2

r_string <- paste("~R^2==~", r2)
x_pos <- max(variant_summary$MeanSNPDist) * 0.85
y_pos <-
  min(variant_summary$NumMultiallelicSite / variant_summary$NumVariantSite)
variant_summary <- rename_datasets(variant_summary)

variant_summary %>%
  mutate(`Dataset size (#)` = `Dataset Size (#)`) %>% 
  ggplot(aes(x = MeanSNPDist, y = NumMultiallelicSite/NumVariantSite)) +
  geom_point(mapping = aes(color = Dataset, size = `Dataset size (#)`)) +
  theme_bw() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "black",
              size = default_line_size) +
  scale_color_manual(values = palette) +
  ylab("All multiallelic sites / all variant sites") +
  xlab("Mean pairwise SNP distance (BP)") +
  guides(color = FALSE) +
  theme(text = element_text(size = default_font_size)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  annotate("text",
           size = rsquared_font_size,
           x = x_pos,
           y = y_pos,
           label = r_string,
           parse = TRUE) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

save_as_pdf_eps_png("Figure_2B_multiallelic_vs_SNP_distance",
                    default_width,
                    default_height)

## 2C: Heatmap of triallelic site SNPEFF impact groupings ----------------------
whiteToOrange = colorRampPalette(c("white", "orange"))
# We want a unique shade of red for each level of the gene
numColors = length(min(tri_functional_heatmap):max(tri_functional_heatmap))
# Generate a character vector of colors with each shade of red
heatMapCols = whiteToOrange(numColors)
pheatmap::pheatmap(tri_functional_heatmap,
                   color = heatMapCols,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE,
                   border_color = "black",
                   fontsize = default_font_size,
                   angle_col = 0,
                   number_format = "%.0f",
                   height = default_height,
                   width = default_width,
                   filename = "figures/Figure_2C_multiallelic_impact.pdf")
pheatmap::pheatmap(tri_functional_heatmap,
                   color = heatMapCols,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE,
                   border_color = "black",
                   fontsize = default_font_size,
                   angle_col = 0,
                   number_format = "%.0f",
                   height = default_height,
                   width = default_width,
                   filename = "figures/Figure_2C_multiallelic_impact.png")

# Figure S2: Multiallelic Sites
## S2A: % Multiallelic vs Number of Isolates -----------------------------
multi_vs_sample_model <-
  lm(formula = 100 * NumMultiallelicSite/NumVariantSite ~ NumIsolates,
     data = variant_summary)
r2 <- round(summary(multi_vs_sample_model)$r.squared, 2)
r_string <- paste("~R^2==~", r2)

x_pos <- max(variant_summary$NumIsolates) * .8
y_pos <- max(variant_summary$NumMultiallelicSite /
               variant_summary$NumVariantSite)
variant_summary %>%
  ggplot(aes(x = NumIsolates,
             y = NumMultiallelicSite/NumVariantSite)) +
  geom_point(mapping = aes(color = Dataset), size = default_dot_size) +
  theme_bw() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "black",
              size = default_line_size) +
  ylab("All multiallelic sites / all variant sites") +
  xlab("Samples") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme(text = element_text(size = default_font_size)) +
  annotate("text",
           size = rsquared_font_size,
           x = x_pos,
           y = y_pos,
           label = r_string,
           parse = TRUE) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_color_manual(values = palette, 
                     labels = c(expression(paste(italic("C. difficile "), "no. 1")), 
                                expression(paste(italic("C. difficile "), "no. 2")), 
                                expression(italic("E. faecalis")), 
                                expression(italic("E. faecium")), 
                                expression(italic("K. pneumoniae")), 
                                expression(italic("L. crispatus")), 
                                expression(paste(italic("S. aureus "), "no. 1")), 
                                expression(paste(italic("S. aureus "), "no. 2")), 
                                expression(italic("S. maltophila")))) +
  theme(legend.text.align = 0)

save_as_pdf_eps_png("Figure_S2A_multiallelic_vs_num_sample",
                    default_width,
                    default_height)


## S2B: SNPEFF impact (multi vs any) -------------------------------------------
impact_level_order <- c("Low", "Moderate", "High")

functional_impact_long <-
  functional_impact_long %>%
  filter(variable != "Baseline")

functional_impact_long$variable2 <-
  factor(functional_impact_long$variable, levels = impact_level_order)

functional_impact_long %>%
  left_join(project_key, by = "Project") %>%
  ggplot(aes(x = variable2,
             y = value)) +
  geom_boxplot(color = "black",
               alpha = 1.0,
               size = default_line_size) +
  ylab("Multiallelic sites / variant sites") +
  xlab("Predicted functional impact") +
  theme_bw() +
  geom_jitter(aes(color = Dataset),
              size = default_dot_size) +
  theme(text = element_text(size = default_font_size)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_color_manual(values = palette, 
                     labels = c(expression(paste(italic("C. difficile "), "no. 1")), 
                                expression(paste(italic("C. difficile "), "no. 2")), 
                                expression(italic("E. faecalis")), 
                                expression(italic("E. faecium")), 
                                expression(italic("K. pneumoniae")), 
                                expression(italic("L. crispatus")), 
                                expression(paste(italic("S. aureus "), "no. 1")), 
                                expression(paste(italic("S. aureus "), "no. 2")), 
                                expression(italic("S. maltophila")))) +
  theme(legend.text.align = 0)

save_as_pdf_eps_png("Figure_S2B_multiallelic_functional_impact",
                    default_width,
                    default_height)

## S2C: Barplot of % of multiallelic sites with mismatching SNPEFF impact ------
functional_variant_summary_df <- as.data.frame(functional_variant_summary)
project_key_df <- as.data.frame(project_key)
functional_variant_summary_key_df <-
  left_join(functional_variant_summary_df, project_key_df, by = "Project")
functional_variant_summary_key_df$Dataset <- as.factor(functional_variant_summary_key_df$Dataset)
functional_variant_summary_key_df %>%
  filter(!is.na(`percent_different`)) %>%
  ggplot() +
  geom_bar(mapping = aes(x = reorder(Dataset, percent_different),
                         y = percent_different / 100,
                         fill = Dataset),
           stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = default_font_size), 
        axis.text.x = element_text(angle = 90, color = "black", hjust = 1), 
        legend.position = "none", 
        axis.text.y = element_text(color = "black")) +
  ylab("Functional impact mismatches") +
  xlab("") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = palette) +
  scale_x_discrete(labels = c(
                      expression(italic("K. pneumoniae")), 
                      expression(italic("S. maltophila")), 
                      expression(italic("E. faecalis")), 
                      expression(paste(italic("S. aureus "), "no. 1")), 
                      expression(italic("L. crispatus")), 
                      expression(italic("E. faecium")), 
                      expression(paste(italic("C. difficile "), "no. 1")),
                      expression(paste(italic("S. aureus "), "no. 2")), 
                      expression(paste(italic("C. difficile "), "no. 2"))))

save_as_pdf_eps_png("Figure_S2C_multiallelic_functional_impact",
                    default_width,
                    default_height)

## S2D: Convergence on tree -----------------------------------------------------
convergence %>%
  mutate(`Locus type` = `Locus Type`) %>% 
  filter(Num_Convergence_Events_8_bin %in%
           c("0", "1", "2", "3", "4", "5", "6", "7", ">7")) %>%
  ggplot(aes(x = Num_Convergence_Events_8_bin,
             y = relative_frequency / 100,
             col = `Locus type`)) +
  geom_boxplot(size = default_line_size / 2) +
  theme_bw() +
  scale_color_manual(values = c("grey","black"), labels = c("Biallelic", "Multiallelic")) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Convergence events per SNP") +
  ylab("Relative frequency") +
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5", "6", "7", ">7")) +
  theme(text = element_text(size = default_font_size)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
save_as_pdf_eps_png("Figure_S2D_convergence",
                    default_width,
                    default_height)

# Read in Reference Allele Data ------------------------------------------------
frac_mismatch = read.table('data/local/ref_allele/frac_mismatch.tsv')
num_mismatch = read.table('data/local/ref_allele/num_mismatch.tsv')
anc_probs = read.table('data/local/ref_allele/anc_probs.tsv')

anc_probs <- anc_probs %>% filter(!dat %in% drop_projects)

# Figure 3: Reference Allele Plots ---------------------------------------------
# 3A: Mismatch counts for each dataset------------------------------------------
fm = melt(frac_mismatch, variable.name = 'object')
fm$Project = rep(rownames(frac_mismatch), ncol(frac_mismatch))
fm$value = fm$value

fm <- fm %>% filter(!Project %in% drop_projects)


fm = fm %>% left_join(project_key, by = 'Project') %>% left_join(color_key, by = 'object')
fm %>%
  ggplot(aes(x = reorder(Dataset, value),
             y = value,
             fill = object)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  ylab('Reference allele mismatch') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(name = "Mismatched allele", 
                    labels = c(expression(paste("Ref. genome" != "major"), 
                               paste("Ref. genome" != "ancestral"), 
                               paste("Major" != "ancestral"))),
                    values = c(color_key$hex_color[color_key$object == "ref_maj"],
                               color_key$hex_color[color_key$object == "ref_anc"],
                               color_key$hex_color[color_key$object == "maj_anc"])) +
  theme_bw() +
  theme(legend.text.align = 0) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size = default_font_size * .8)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  xlab('') +
  scale_x_discrete(labels = c(
    expression(paste(italic("S. aureus "), "no. 1")), 
    expression(italic("K. pneumoniae")), 
    expression(paste(italic("S. aureus "), "no. 2")), 
    expression(italic("E. faecalis")), 
    expression(paste(italic("C. difficile "), "no. 2")), 
    expression(paste(italic("C. difficile "), "no. 1")),
    expression(italic("S. maltophila")), 
    expression(italic("E. faecium")), 
    expression(italic("L. crispatus"))))
save_as_pdf_eps_png("Figure_3A_allele_mismatch_correct_equals_sign",
                    default_width,
                    default_height)

# 3B: Ancestral state confidence for each dataset-------------------------------
mean_mismatches = fm %>% group_by(Dataset) %>% summarise(mean(value))
dataset_order = mean_mismatches$Dataset[order(mean_mismatches$`mean(value)`)]


names(anc_probs)[1] = 'Project'
anc_probs = anc_probs %>% left_join(project_key, by = 'Project')

anc_probs = anc_probs %>%
  group_by(Dataset,low_conf = anc_prob < 0.875) %>%
  summarise(count = n()) %>%
  mutate(perc = count/sum(count))

anc_probs$hex_color = sapply(anc_probs$Dataset,function(x) project_key$hex_color[project_key$Dataset == x])

ggplot(anc_probs[anc_probs$low_conf == TRUE,],
       aes(x = reorder(Dataset, perc),y = perc, fill = Dataset)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  #scale_fill_manual(values = palette) +
  #scale_x_discrete(limits = dataset_order) +
  xlab('') +
  ylab('Low-confidence ancestral reconstruction') +
  theme_bw() +
  theme(text = element_text(size = default_font_size * .8)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_fill_manual(values = palette) + 
  scale_x_discrete(limits = dataset_order, 
                   labels = c(
    expression(paste(italic("S. aureus "), "no. 1")), 
    expression(italic("K. pneumoniae")), 
    expression(paste(italic("S. aureus "), "no. 2")), 
    expression(italic("E. faecalis")), 
    expression(paste(italic("C. difficile "), "no. 2")), 
    expression(paste(italic("C. difficile "), "no. 1")),
    expression(italic("S. maltophila")), 
    expression(italic("E. faecium")), 
    expression(italic("L. crispatus"))))
save_as_pdf_eps_png("Figure_3B_ancestral_reconstruction",
                    default_width * .5,
                    default_height)

# Read in Overlap Data ---------------------------------------------------------
overlap_func_impact_long <- read_tsv("data/local/overlap/overlap_diversity_long.tsv")
overlap_snpeff_summary <- read_tsv("data/local/overlap/overlap_snpeff_summary.tsv")
overlap_heatmap <- read_tsv("data/local/overlap/overlap_heatmap.tsv")
overlap_stats <- read_tsv("data/local/overlap/overlap_stats.tsv")

overlap_stats <- overlap_stats %>% filter(!Project %in% drop_projects)
overlap_func_impact_long <- overlap_func_impact_long %>%
  filter(!Project %in% drop_projects)
overlap_snpeff_summary <- overlap_snpeff_summary %>%
  filter(!Project %in% drop_projects)
overlap_heatmap <- NULL

# Figure 4: Overlap Plot -------------------------------------------------------
long_overlap_snpeff_summary <- overlap_snpeff_summary %>%
  select(NumHighModerate, NumHighLow, NumModerateLow, Project, NumOverlapSites)
long_overlap_snpeff_summary <- long_overlap_snpeff_summary %>%
  mutate("PercHighModerate" = NumHighModerate / NumOverlapSites,
         "PercHighLow" = NumHighLow / NumOverlapSites,
         "PercModerateLow" = NumModerateLow / NumOverlapSites) %>%
  select(-c(NumHighModerate, NumHighLow, NumModerateLow, NumOverlapSites))
long_overlap_snpeff_summary <- melt(long_overlap_snpeff_summary, variable.name = 'object')
long_overlap_snpeff_summary <- left_join(long_overlap_snpeff_summary, project_key, by = "Project")
long_overlap_snpeff_summary <-
  long_overlap_snpeff_summary %>% filter(!Project %in% drop_projects)

long_overlap_snpeff_summary %>%
  ggplot(aes(x = reorder(Dataset, value),
             y = value,
             fill = object)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  xlab('') +
  ylab('Functional impact mismatch') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(name = "Predicted functional impact",
                    labels = c("High & moderate", "High & low", "Moderate & low"),
                    values = c(color_key$hex_color[color_key$object == "mod_hi"],
                               color_key$hex_color[color_key$object == "low_hi"],
                               color_key$hex_color[color_key$object == "low_mod"])) +
  theme_bw() +
  theme(text = element_text(size = default_font_size)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  theme(legend.title = element_text(size = default_font_size * .6),
        legend.text = element_text(size = default_font_size * .6)) + 
  scale_x_discrete(labels = c(
    expression(paste(italic("C. difficile "), "no. 1")),
    expression(italic("S. maltophila")), 
    expression(paste(italic("C. difficile "), "no. 2")),
    expression(paste(italic("S. aureus "), "no. 2")), 
    expression(paste(italic("S. aureus "), "no. 1")), 
    expression(italic("E. faecium")), 
    expression(italic("K. pneumoniae")), 
    expression(italic("L. crispatus")),
    expression(italic("E. faecalis"))))
    
save_as_pdf_eps_png("Figure_4A_overlap_functional_impact",
                    default_width,
                    default_height)


# Figure S4 Heritability Estimate Differences
cols = project_key 
her_ests = read.csv('data/local/heritability/her_ests_shared.csv', header = F, stringsAsFactors = F)
colors = cols$hex_color
names(colors) = cols$Project
proj_name = cols$Dataset
names(proj_name) = cols$Project

names(her_ests) = c('datname','phen_type','ref_type','multi_type','her_est')

her_ests$dataset = gsub('_prewas.*','',her_ests$datname)
her_ests$phen_type[her_ests$phen_type == 'bm'] = 'BM'
her_ests$phen_type[her_ests$phen_type == 'wn'] = 'WN'
her_ests$ref_type[her_ests$ref_type == 'anc'] = 'Ancestral allele'
her_ests$ref_type[her_ests$ref_type == 'maj'] = 'Major allele'
her_ests$ref_type[her_ests$ref_type == 'ref'] = 'Ref. genome allele'
her_ests$multi_type[her_ests$multi_type == 'multi'] = 'With multiallelic sites'
her_ests$multi_type[her_ests$multi_type == 'nomulti'] = 'Without multiallelic sites'

her_ests %>% 
  ggplot(aes(x = phen_type, y = her_est * 100)) + 
  geom_jitter(aes(col = dataset), size = 3) + 
  facet_grid(multi_type ~ ref_type) +
  xlab('Phenotype evolutionary model') + 
  ylab('Heritability estimate (%)') +
  scale_color_manual(name = 'Dataset',
                     values = colors,
                     #labels = proj_name) +
                     labels = c(expression(paste(italic("C. difficile "), "no. 1")), 
                                expression(paste(italic("C. difficile "), "no. 2")), 
                                expression(italic("E. faecium")), 
                                expression(italic("E. faecalis")), 
                                expression(italic("L. crispatus")), 
                                expression(paste(italic("S. aureus "), "no. 1")), 
                                expression(paste(italic("S. aureus "), "no. 2")))) + 
  theme_bw() + 
  theme(text = element_text(size = 10)) + 
  theme(legend.text.align = 0)


save_as_pdf_eps_png("Figure_S4_heritability_estimate_differences", 6, 4)

# Figure S5 --------------------------------------------------------------------
## S5A Number of Overlapping Genes ----------------------------------------------
overlap_stats <- left_join(overlap_stats, project_key, by = "Project")
overlap_stats <- overlap_stats %>% filter(!is.na(Dataset))
overlap_stats <- left_join(overlap_stats, project_key,
                           by = c("Project", "Dataset", "hex_color", "color"))

overlap_stats %>%
  ggplot() +
  geom_bar(mapping = aes(x = reorder(Dataset, Overlapping_SNP_loci),
                         y = Overlapping_SNP_loci,
                         color = Dataset),
           fill = "white",
           size = 2,
           stat = "identity") +
  theme_bw() +
  scale_color_manual(values = palette) +
  theme(text = element_text(size = default_font_size)) +
  ylab("SNPs found in overlapping genes") +
  xlab("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_x_discrete(labels = c(
    expression(italic("E. faecium")), 
    expression(italic("E. faecalis")),
    expression(italic("L. crispatus")),
    expression(paste(italic("S. aureus "), "no. 2")), 
    expression(italic("K. pneumoniae")), 
    expression(paste(italic("S. aureus "), "no. 1")), 
    expression(italic("S. maltophila")), 
    expression(paste(italic("C. difficile "), "no. 2")),
    expression(paste(italic("C. difficile "), "no. 1"))))

save_as_pdf_eps_png("Figure_S5A_overlap_SNP_count",
                    default_width,
                    default_height)

## S5B Number of Overlapping Positions ------------------------------------------
overlap_stats %>%
  select(c(Overlapping_Genes_w_SNPs, Dataset, hex_color)) %>%
  ggplot() +
  geom_bar(mapping = aes(x = reorder(Dataset, Overlapping_Genes_w_SNPs),
                         y = Overlapping_Genes_w_SNPs,
                         color = Dataset),
           fill = "white",
           size = 2,
           stat = "identity") +
  theme_bw() +
  theme(text = element_text(size = default_font_size)) +
  scale_color_manual(values = palette)  +
  ylab("Overlapping genes with SNPs") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_x_discrete(labels = c(
    expression(italic("E. faecium")), 
    expression(italic("E. faecalis")),
    expression(italic("L. crispatus")),
    expression(paste(italic("S. aureus "), "no. 2")), 
    expression(paste(italic("S. aureus "), "no. 1")), 
    expression(italic("K. pneumoniae")), 
    expression(italic("S. maltophila")), 
    expression(paste(italic("C. difficile "), "no. 1")),
    expression(paste(italic("C. difficile "), "no. 2"))))

save_as_pdf_eps_png("Figure_S5B_overlap_gene_count",
                    default_width,
                    default_height)

## S6 Resource Usage -----------------------------------------------------------
usage_df <- read_csv("data/hpc/resource_usage/prewas_resource_usage.csv")
colnames(usage_df)[3] <- "Dataset"
for (i in 1:nrow(usage_df)) {
  for  (j in 1:nrow(project_key)) {
    if (usage_df$Dataset[i] == project_key$Project[j]) {
      usage_df$Dataset[i] <- project_key$Dataset[j]
    }
  }
}


resource_palette <- palette[names(palette) %in% usage_df$Dataset] 
storage.mode(usage_df$`Cores (#)`) <- "character"
get_legend <- function(myggplot){
  # This function is from http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

usage_df$Tree[usage_df$Tree == "No Tree Provided"] <- "No tree provided"
usage_df$Tree[usage_df$Tree == "Tree Provided"] <- "Tree provided"
usage_df$Method[usage_df$Method == "Major Allele"] <- "Major allele"
usage_df$Method[usage_df$Method == "Ancestral Reconstruction"] <- "Ancestral reconstruction"

# Remove cores == 4
usage_df <- usage_df %>% filter(`Cores (#)` != "4")

memory_plot <- usage_df %>% 
  ggplot(aes(x = `Variants (#)`, y = `Memory (GB)`)) + 
  geom_point(aes(shape = Tree, 
                 size = `Cores (#)`, 
                 color = Dataset)) + 
  theme_bw() + 
  ylim(0, 30) + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 15),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black")) + 
  scale_color_manual(values = resource_palette) +
  scale_size_discrete(range = c(1.5, 3)) + 
  facet_wrap(~ Method) +
  scale_y_continuous(minor_breaks = seq(1, 30, 1))

time_plot <- usage_df %>% 
  mutate(`Time (hours)` = `Time (Hours)`) %>% 
  ggplot(aes(x = `Variants (#)`, y = `Time (hours)`)) + 
  geom_point(aes(shape = Tree, 
                 size = `Cores (#)`, 
                 color = Dataset)) + 
  theme_bw() + 
  ylim(0, max(usage_df$`Time (Hours)`)) + 
  theme(text = element_text(size = 15),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black")) + 
  scale_size_discrete(range = c(1.5, 3)) + 
  facet_wrap(~ Method) +
  scale_y_continuous(minor_breaks = seq(1, 40, 1)) + 
  scale_color_manual(values = resource_palette, 
                     labels = c(expression(paste(italic("C. difficile "), "no. 1")), 
                                expression(paste(italic("C. difficile "), "no. 2")), 
                                expression(italic("E. faecalis")), 
                                expression(italic("E. faecium")), 
                                expression(italic("K. pneumoniae")), 
                                expression(italic("L. crispatus")), 
                                expression(paste(italic("S. aureus "), "no. 1")), 
                                expression(paste(italic("S. aureus "), "no. 2")), 
                                expression(italic("S. maltophila")))) +
  theme(legend.text.align = 0)

legend <- get_legend(time_plot)

time_plot <- time_plot + theme(legend.position = "none")
pdf("figures/Figure_S6_prewas_resource_usage.pdf", width = 16, height = 6)
grid.arrange(memory_plot, time_plot, legend, ncol = 3)
dev.off()

png("figures/Figure_S6_prewas_reource_usage.png", width = 16, height = 6, units = "in", res = 500)
grid.arrange(memory_plot, time_plot, legend, ncol = 3)
dev.off()

sink("data/local/table/numbers_for_paper.txt")
'Number of datasets:'
nrow(variant_summary)
'Number of species:'
length(unique(gsub('#.*','',variant_summary$Dataset)))
'Total sample size of all datasets:'
sum(variant_summary$NumIsolates)
'Multiallelic variant sites:'
'Range:'
'Count:'
range(variant_summary$NumMultiallelicSite)
'Percent:'
paste0(round(range(variant_summary$NumMultiallelicSite/variant_summary$NumVariantSite)*100,1),'%')
'R^2 of mean pairwise SNP distance vs. multiallelic variant sites:'
r2_mean
'Mismatch between reference allele choice methods:'
'Range:'
'Count:'
range(num_mismatch)
'Percent:'
range(frac_mismatch)
'Mismatched predicted functional impact for overlapping genes:'
'Range:'
'Count:'
range(overlap_snpeff_summary$NumDifferentImpact)
'Percent:'
range(overlap_snpeff_summary$NumDifferentImpact/overlap_snpeff_summary$NumOverlapSites)
'Mismatched predicted functional impact for multiallelic sites:'
'Range:'
'Count:'
range(tri_functional_summary$NumDifferentImpact,na.rm = T)
'Percent:'
range(tri_functional_summary$NumDifferentImpact/tri_functional_summary$NumMultiallelicSites,na.rm = T)
sink()
