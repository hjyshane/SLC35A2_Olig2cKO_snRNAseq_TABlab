#### Library load ####
# Analysis
library(Seurat)
library(SeuratExtend)
library(glmGamPoi)
library(presto)

# File load and save
library(qs)

# Data manipulation
library(tidyverse)

source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/07_clustering_oligsub.R')
set.seed(0827)

save_dir <-
    '~/SLC35A2_Olig2cKO_snRNA/WIP'
qsave_dir <-
    file.path(save_dir, "QSAVE")
csv_dir <-
    file.path(save_dir, 'CSV')

# Read previous object if needed.
oligsub <- qs::qread(file.path(qsave_dir, "12_olig_processed.qs"))

oligsub <- clustering_oligsub(oligsub,
                                 range_start = 0.1,
                                 range_end = 1.2,
                                 range_step = 0.1,
                                 save = TRUE,
                                 qsave_dir = qsave_dir,
                                 markers_dir = file.path(csv_dir, "Markers"))