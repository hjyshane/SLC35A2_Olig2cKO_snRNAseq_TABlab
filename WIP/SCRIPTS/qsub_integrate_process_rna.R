#### Library load ####
# Analysis
library(Seurat)
library(SeuratExtend)
library(glmGamPoi)

# File load and save
library(qs)

# Load HDF5 library from IGM
dyn.load("/igm/apps/hdf5/hdf5-1.12.1/lib/libhdf5_hl.so.200")
library(hdf5r)

# Data manipulation
library(tidyverse)

source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/06__process_integrated.R')
set.seed(0827)

save_dir <-
    '~/SLC35A2_Olig2cKO_snRNA/WIP'
qsave_dir <-
    file.path(save_dir, "QSAVE")

# Read previous object if needed.
integrated <- qs::qread(file.path(qsave_dir, "05_integrated_rna.qs"))

# Integration
integrated <- process_integrated(integrated,
                                 save = TRUE,
                                 qsave_dir = qsave_dir)
