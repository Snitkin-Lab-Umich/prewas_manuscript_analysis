# Note: Kp_2, Cd_1, and Cd_2 must be excluded from draft 1 of the paper
drop_projects <- c("Kp_2", "Cd_1", "Cd_2")

# Libraries --------------------------------------------------------------------
library(magrittr)
library(ggplot2)
source("lib/multiallelic_lib.R")

# Color palettes and project names ---------------------------------------------
color_key <- read_tsv("data/plot_colors.tsv")
project_key <- read_tsv("data/project_names_and_colors.tsv")
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
  ggsave(filename = paste0("../figures/", file_prefix, ".pdf"),
         width = width_size,
         height = height_size,
         units = "in",
         device = "pdf")
  ggsave(filename = paste0("../figures/", file_prefix, ".eps"),
         width = width_size,
         height = height_size,
         units = "in",
         device = "eps")
  ggsave(filename = paste0("../figures/", file_prefix, ".png"),
         width = width_size,
         height = height_size,
         units = "in",
         device = "png")
}

# Read in Multiallelic Data ----------------------------------------------------
variant_summary <- read_tsv("../data/multiallelic_summary.tsv")
snp_mat_summary <- read_tsv("../data/multiallelic_summary_subsampled.tsv")
snp_mat_grouped <- read_tsv("../data/multiallelic_summary_by_project.tsv")


functional_impact_long <- read_tsv("../data/SNP_diversity_long.tsv")
functional_variant_summary <- read_tsv(file = "../data/snpeff_multiallelic_summary.tsv")
tri_functional_summary <- read_tsv(file = "../data/snpeff_triallelic_summary.tsv")
tri_functional_heatmap <- read.csv(file = "../data/triallelic_heatmap.csv", row.names = 1)
convergence <- read_tsv(file = "../data/convergence.tsv")

# Filter out datasets to exclude:
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
snp_mat_grouped %>%
  ggplot(aes(x = Percent / 100,
             y = AvgNumMultiallelicSite / AvgNumVariantSite,
             color = Dataset)) +
  geom_line(size = default_line_size) +
  geom_point(size = 4) +
  theme_bw() +
  scale_y_continuous(labels=scales::percent_format(accuracy = .1)) +
  scale_x_continuous(labels=scales::percent) +
  scale_color_manual(values = palette) +
  ylab("All Multiallelic Sites / All Variant Sites") +
  xlab("Dataset Subset") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  theme(text = element_text(size = default_font_size))

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
variant_summary %>%
  ggplot(aes(x = MeanSNPDist, y = NumMultiallelicSite/NumVariantSite)) +
  geom_point(mapping = aes(color = Dataset, size = `Dataset Size (#)`)) +
  theme_bw() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "black",
              size = default_line_size) +
  scale_color_manual(values = palette) +
  ylab("All Multiallelic Sites / All Variant Sites") +
  xlab("Mean Pairwise SNP Distance (BP)") +
  scale_x_continuous(labels = scales::comma) +
  guides(color = FALSE) +
  theme(text = element_text(size = default_font_size)) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 0.1)) +
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
                   filename = "../figures/Figure_2C_multiallelic_impact.pdf")
#dev.off()
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
                   filename = "../figures/Figure_2C_multiallelic_impact.png")
#dev.off()

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
  ylab("All Multiallelic Sites / All Variant Sites") +
  xlab("Samples") +
  scale_color_manual(values = palette) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 0.1)) +
  theme(text = element_text(size = default_font_size)) +
  annotate("text",
           size = rsquared_font_size,
           x = x_pos,
           y = y_pos,
           label = r_string,
           parse = TRUE) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

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
  ylab("Multiallelic Sites / Variant Sites") +
  xlab("Predicted Functional Impact") +
  theme_bw() +
  geom_jitter(aes(color = Dataset),
              size = default_dot_size) +
  scale_color_manual(values = palette) +
  theme(text = element_text(size = default_font_size)) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

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
  theme(text = element_text(size = default_font_size)) +
  ylab("Functional impact mismatches") +
  xlab("") +
  scale_fill_manual(values = palette)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none") +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

save_as_pdf_eps_png("Figure_S2C_multiallelic_functional_impact",
                    default_width,
                    default_height)

## S2D: Convergence on tree -----------------------------------------------------
convergence %>%
  filter(Num_Convergence_Events_8_bin %in%
           c("0", "1", "2", "3", "4", "5", "6", "7", ">7")) %>%
  ggplot(aes(x = Num_Convergence_Events_8_bin,
             y = relative_frequency / 100,
             col = `Locus Type`)) +
  geom_boxplot(size = default_line_size / 2) +
  theme_bw() +
  scale_color_manual(values = c("grey","black")) +
  scale_y_continuous(labels=scales::percent) +
  xlab("Convergence Events Per SNP") +
  ylab("Relative Frequency") +
  scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5", "6", "7", ">7")) +
  theme(text = element_text(size = default_font_size)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
save_as_pdf_eps_png("Figure_S2D_convergence",
                    default_width,
                    default_height)

# Read in Reference Allele Data ------------------------------------------------
frac_mismatch = read.table('../data/frac_mismatch.tsv')
num_mismatch = read.table('../data/num_mismatch.tsv')
anc_probs = read.table('../data/anc_probs.tsv')

anc_probs <- anc_probs %>% filter(!dat %in% drop_projects)

# Figure 3: Reference Allele Plots ---------------------------------------------

# 3A: Mismatch counts for each dataset------------------------------------------
fm = melt(frac_mismatch, variable.name = 'object')
fm$Project = rep(rownames(frac_mismatch), ncol(frac_mismatch))
fm$value = fm$value

fm <- fm %>% filter(!Project %in% drop_projects)


fm = fm %>% left_join(project_key, by='Project') %>% left_join(color_key,by='object')
cairo_pdf(paste0('../figures/Figure_3A_allele_mismatch_correct_equals_sign.pdf'),width = 8, height = 6)#,device=cairo_pdf)
fm %>%
  ggplot(aes(x = reorder(Dataset, value),
             y = value,
             fill = object)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  xlab('') +
  ylab('Reference allele mismatch') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(name = "Mismatched Allele",labels = c("Ref. Genome \u2260 Major", "Ref. Genome \u2260 Ancestral", "Major \u2260 Ancestral"),
                    values = c(color_key$hex_color[color_key$object == "ref_maj"],
                               color_key$hex_color[color_key$object == "ref_anc"],
                               color_key$hex_color[color_key$object == "maj_anc"])) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size = default_font_size * .8)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

# 3B: Ancestral state confidence for each dataset-------------------------------
mean_mismatches = fm %>% group_by(Dataset) %>% summarise(mean(value))
dataset_order = mean_mismatches$Dataset[order(mean_mismatches$`mean(value)`)]


names(anc_probs)[1] = 'Project'
anc_probs = anc_probs %>% left_join(project_key, by='Project')

anc_probs = anc_probs %>%
  group_by(Dataset,low_conf=anc_prob < 0.875) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

anc_probs$hex_color = sapply(anc_probs$Dataset,function(x) project_key$hex_color[project_key$Dataset == x])

ggplot(anc_probs[anc_probs$low_conf == TRUE,],
       aes(x = reorder(Dataset, perc),y = perc, fill = Dataset)) +
  geom_bar(stat ='identity') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = palette) +
  scale_x_discrete(limits = dataset_order) +
  xlab('') +
  ylab('Low confidence ancestral reconstruction') +
  theme_bw() +
  theme(text = element_text(size = default_font_size * .8)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
save_as_pdf_eps_png("Figure_3B_ancestral_reconstruction",
                    default_width * .5,
                    default_height)

# Read in Overlap Data ---------------------------------------------------------
overlap_func_impact_long <- read_tsv("../data/overlap_diversity_long.tsv")
overlap_snpeff_summary <- read_tsv("../data/overlap_snpeff_sumary.tsv")
overlap_heatmap <- read_tsv("../data/overlap_heatmap.tsv")
overlap_stats <- read_tsv("../data/overlap_stats.tsv")

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
  scale_fill_manual(name = "Predicted Functional Impact",
                    labels = c("High & Moderate", "High & Low", "Moderate & Low"),
                    values = c(color_key$hex_color[color_key$object == "mod_hi"],
                               color_key$hex_color[color_key$object == "low_hi"],
                               color_key$hex_color[color_key$object == "low_mod"])) +
  theme_bw() +
  theme(text = element_text(size = default_font_size)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  theme(legend.title = element_text(size = default_font_size * .6),
        legend.text = element_text(size = default_font_size * .6))
save_as_pdf_eps_png("Figure_4A_overlap_functional_impact",
                    default_width,
                    default_height)

## No longer included: SNPEFF impact (Overlap vs any variant) split by diveristy (hi vs low)----
# impact_level_order <- c("Low", "Moderate", "High")
# overlap_func_impact_long <-
#   overlap_func_impact_long %>%
#   filter(variable != "Baseline")
#
# overlap_func_impact_long$variable2 <-
#   factor(overlap_func_impact_long$variable, levels = impact_level_order)
#
# overlap_func_impact_long %>%
#   left_join(project_key, by = "Project") %>%
#   ggplot(aes(x = variable2, y = value * 100)) +
#   geom_boxplot(color = "black",
#                alpha = 1.0) +
#   ylab("Frequency at Overlapping Sites (%)") +
#   xlab("Predicted Functional Impact") +
#   theme_bw() +
#   geom_jitter(aes(color = Dataset)) +
#   scale_color_manual(values = palette) +
#   theme(axis.text.x = element_text(color = "black"),
#         axis.text.y = element_text(color = "black"))


## No longer included: Barplot of % of overlapping sites with mismatching SNPEFF impact-------------
# overlap_snpeff_summary %>%
#   left_join(project_key, by = "Project") %>%
#   select(c(percent_different, Dataset, hex_color)) %>%
#   ggplot() +
#   geom_bar(mapping = aes(x = reorder(Dataset, percent_different),
#                          y = percent_different,
#                          color = Dataset),
#            fill = "white",
#            size = 2,
#            stat = "identity") +
#   scale_color_manual(values = palette)  +
#   theme_bw() +
#   ylab("Functional impact mismatches (% overlaping sites)") +
#   xlab("") +
#   ylim(c(0.0, 100)) +
#   theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(axis.text.x = element_text(color = "black"),
#         axis.text.y = element_text(color = "black"))
## Heatmap of all overlapping site SNPEFF impact groupings
# overlap_heatmap <- as.matrix(overlap_heatmap[1:3, 2:4])
# row.names(overlap_heatmap) <- c("High", "Moderate", "Low")
#
# whiteToRed = colorRampPalette(c("white", "purple"))
# # We want a unique shade of red for each level of the gene
# numColors = length(min(overlap_heatmap):max(overlap_heatmap))
# # Generate a character vector of colors with each shade of red
# heatMapCols = whiteToRed(numColors)
#
# pheatmap::pheatmap((overlap_heatmap),
#                    color = heatMapCols,
#                    cluster_rows = FALSE,
#                    cluster_cols = FALSE,
#                    display_numbers = TRUE,
#                    border_color = "black",
#                    number_format = "%.0f",
#                    fontsize = 20,
#                    angle_col = 0,
#                    filename = "../figures/overlap_mismatch_impact_heatmap.pdf")
#
# dev.off()


# Figure S3 --------------------------------------------------------------------
## S3A Number of Overlapping Genes ----------------------------------------------
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
  scale_color_manual(values = palette)  +
  theme(text = element_text(size = default_font_size)) +
  ylab("SNPs found in overlapping genes") +
  xlab("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"))

save_as_pdf_eps_png("Figure_S4A_overlap_SNP_count",
                    default_width,
                    default_height)

## S3B Number of Overlapping Positions ------------------------------------------
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
        axis.text.y = element_text(color = "black"))
save_as_pdf_eps_png("Figure_S4B_overlap_gene_count",
                    default_width,
                    default_height)


sink("../data/numbers_for_paper.txt")
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

