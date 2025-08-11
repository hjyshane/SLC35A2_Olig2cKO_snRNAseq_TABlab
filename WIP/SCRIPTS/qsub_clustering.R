#### Library load ####
# Analysis
library(Seurat)
library(SeuratExtend)
library(glmGamPoi)
library(presto)

# File load and save
library(qs)

# Load HDF5 library from IGM
dyn.load("/igm/apps/hdf5/hdf5-1.12.1/lib/libhdf5_hl.so.200")
library(hdf5r)

# Data manipulation
library(tidyverse)

source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/07__clustering.R')
set.seed(0827)

save_dir <-
    '~/SLC35A2_Olig2cKO_snRNA/WIP'
qsave_dir <-
    file.path(save_dir, "QSAVE")
csv_dir <-
    file.path(save_dir, 'CSV')

# Read previous object if needed.
integrated <- qs::qread(file.path(qsave_dir, "06_integrated_rna.qs"))

# Integration
integrated <- clustering(integrated,
                         range_start = 0.4,
                         range_end = 1.2,
                         range_step = 0.1,
                         save = TRUE,
                         qsave_dir = qsave_dir,
                         markers_dir = file.path(csv_dir, "Markers"))