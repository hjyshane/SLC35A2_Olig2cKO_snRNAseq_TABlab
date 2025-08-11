#### Library load ####
# Analysis
library(Seurat)

# File load and save
library(qs)

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
integrated <- qs::qread(file.path(qsave_dir, "08_clustered_metadata_added.qs"))

# Read reference object for cell typing
ref_obj <- readRDS("~/SLC35A2_Olig2cKO_snRNA/WIP/mouseatlas_cut_processed.rds")

Seurat::DefaultAssay(integrated) <- "SCT"
Seurat::DefaultAssay(ref_obj) <- "SCT"

# Find anchors between sobj and ref_obj
anchors <- Seurat::FindTransferAnchors(
    reference = ref_obj,
    query = integrated,
    dims = 1:30,
    normalization.method = 'SCT')

# Predict cell types
predictions <- Seurat::TransferData(
    anchorset = anchors,
    refdata = ref_obj$Description,   # or whichever column contains cell type labels
    dims = 1:30
)

# Add metadata
integrated <- Seurat::AddMetaData(
    object = integrated,
    metadata = predictions)

qs::qsave(integrated, file = file.path(qsave_dir, "09_ref_annotated_obj.qs"))