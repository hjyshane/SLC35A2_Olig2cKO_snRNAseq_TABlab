library(tidyverse)
library(Seurat)
library(ggplot2)

# Cell count data frame
cell_counts <-
    table(oligsub_man2[["ref_short"]][,1], oligsub_man2[["TimePoint"]][,1]) %>%
    as.data.frame.matrix() %>%
    tibble::rownames_to_column(var = "Cell_Type")

# Filter
clusters_to_keep <- cell_counts %>%
    filter(P00_cKO >= 10 & P00_Control >= 10 |
               P08_cKO >= 10 & P08_Control >= 10 |
               P23_cKO >= 10 & P23_Control >= 10 |
               P56_cKO >= 10 & P56_Control >= 10) %>%
    pull(Cell_Type)

oligsub_man2_filtered <- subset(oligsub_man2, ref_short %in% clusters_to_keep)
cell_counts2 <-
    table(oligsub_man2_filtered[["ref_short"]][,1], oligsub_man2_filtered[["TimeGeno"]][,1]) %>%
    as.data.frame.matrix() %>%
    tibble::rownames_to_column(var = "Cell_Type")
