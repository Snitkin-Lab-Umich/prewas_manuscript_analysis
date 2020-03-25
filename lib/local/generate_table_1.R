# This script aggregates data together for both Table_1 and Table_S1 for the 
# paper. 

library(tidyverse)
table1 <- read_tsv("data/local/multiallelic/multiallelic_summary.tsv")
table1 <- table1 %>%
  select(Dataset, 
         NumIsolates,
         NumMultiallelicSite,
         MeanSNPDist, 
         NumBiallelicOverlappingGene)
colnames(table1) <- c("Name", 
                      "Samples (Count)", 
                      "Multiallelic Sites (Count)", 
                      "Mean SNP Distance (BP)", 
                      "SNPs in overlapping genes (Count)")

table1_part2 <- read_tsv("data/key/handwritten_Table_1.tsv")
table1_part2 <-
  table1_part2 %>% select("Name",
                          "Bioproject",
                          "Reference_genome_biosample",
                          "Paper",
                          "Dataset_description")
colnames(table1_part2) <- c("Name",
                            "Bioproject",
                            "Reference Genome Biosample",
                            "Paper",
                            "Dataset Description")
table1 <- left_join(table1_part2, table1, by = c("Name"))

write_tsv(table1, path = "data/local/table/Table_1.tsv")


