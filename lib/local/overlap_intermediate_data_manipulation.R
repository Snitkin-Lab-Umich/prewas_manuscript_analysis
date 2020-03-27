# Goal: find all of the SNPs in overlapping positions per dataset.
# Data summarized in this analysis is plotted in Figure 4 and Figure S5. 
library(tidyverse)
source("lib/local/multiallelic_lib.R")

drop_projects <- c("Kp_2", "Cd_1", "Cd_2")
variant_summary <- read_tsv("data/local/multiallelic/multiallelic_summary.tsv")

variant_summary <- variant_summary %>%
  select(Project, Species, Dataset, hex_color, color,
         NumVariantSite, NumBiallelicOverlappingGene)

# Reformat data for plotting ---------------------------------------------------
overlap_snpeff_summary <-
  read_tsv(file = "data/local/overlap/snpeff_overlap_summary.tsv")
overlap_snpeff_htmp_data <-
  overlap_snpeff_summary %>%
  filter(!Project %in% drop_projects) %>%
  select(NumHighHigh,
         NumHighModerate,
         NumHighLow,
         NumModerateModerate,
         NumModerateLow,
         NumLowLow) %>%
  colSums(na.rm = TRUE)
overlap_snpeff_htmp <- matrix(0, ncol = 3, nrow = 3)
colnames(overlap_snpeff_htmp) <-
  row.names(overlap_snpeff_htmp) <- c("High", "Moderate", "Low")
overlap_snpeff_htmp[1, 1] <- overlap_snpeff_htmp_data[1]
overlap_snpeff_htmp[2, 2] <- overlap_snpeff_htmp_data[4]
overlap_snpeff_htmp[3, 3] <- overlap_snpeff_htmp_data[6]

overlap_snpeff_htmp[1, 2] <-
  overlap_snpeff_htmp[2, 1] <-
  overlap_snpeff_htmp_data[2]
overlap_snpeff_htmp[1, 3] <-
  overlap_snpeff_htmp[3, 1] <-
  overlap_snpeff_htmp_data[3]

overlap_snpeff_htmp[2, 3] <-
  overlap_snpeff_htmp[3, 2] <-
  overlap_snpeff_htmp_data[5]

overlap_snpeff_summary <- overlap_snpeff_summary %>% filter(!is.na(Project))

overlap_variant_summary <-
  left_join(overlap_snpeff_summary,
            variant_summary,
            by = c("Project"))
overlap_variant_summary <- overlap_variant_summary %>% filter(!is.na(Dataset))

overlap_variant_summary <- overlap_variant_summary %>%
  mutate(# Baseline = NumOverlapSites / NumVariantSite,
         Low = (NumHighLow + NumModerateLow + NumLowLow) / NumOverlapSites,
         Moderate = (NumHighModerate + NumModerateModerate + NumModerateLow) / 
           NumOverlapSites,
         High = (NumHighHigh + NumHighLow + NumHighModerate) / NumOverlapSites)

melted_overlap_variant_summary <- melt(overlap_variant_summary,
                                       id.vars = 'Project',
                                       measure.vars = c( #'Baseline',
                                                        'Low',
                                                        'Moderate',
                                                        'High'))

melted_overlap_variant_summary <- melted_overlap_variant_summary %>%
  filter(!is.na(Project))

bad_ones <- melted_overlap_variant_summary %>% filter(is.infinite(value))
melted_overlap_both_variant <-
  melted_overlap_variant_summary %>% filter(!is.infinite(value))

overlap_snpeff_htmp <- rownames_to_column(as.data.frame(overlap_snpeff_htmp),
                                          var = "functional_impact")

write_tsv(melted_overlap_both_variant, "data/local/overlap/overlap_diversity_long.tsv")
write_tsv(overlap_snpeff_summary, "data/local/overlap/overlap_snpeff_summary.tsv")
write_tsv(overlap_snpeff_htmp, "data/local/overlap/overlap_heatmap.tsv")

