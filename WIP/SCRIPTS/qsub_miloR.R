#### Library load ####
BiocManager::install(c("miloR", "scater", "scran"))

# library load
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(tidyverse)
library(patchwork)

# File load and save
library(qs)

set.seed(0827)

save_dir <-
    '~/SLC35A2_Olig2cKO_snRNA/WIP'
qsave_dir <-
    file.path(save_dir, "QSAVE")

# Read previous object if needed.
integrated <- qs::qread(file.path(qsave_dir, "WIP/QSAVE/10_metadata_edited_obj.qs"))

integrated_sce <- as.SingleCellExperiment(integrated)
message("sce")

integrated_milo <- Milo(integrated_sce)
message("milo object created")

integrated_milo <- buildGraph(integrated_milo, k = 10, d = 30)
message("knn constructed")

integrated_milo <- makeNhoods(integrated_milo, prop = 0.1, k = 10, d=30, refined = TRUE)
message("neighbourhoods defined")

integrated_milo <- calcNhoodDistance(integrated_milo, d=30)
message("neighborhood distance calculated")

integrated_milo <- countCells(integrated_milo, meta.data = data.frame(colData(integrated_milo)), samples="MouseID")
message("cell counted")

qs::qsave(integrated_milo, file = "~/SLC35A2_Olig2cKO_snRNA/WIP/QSAVE/99_miloRobject.qs")