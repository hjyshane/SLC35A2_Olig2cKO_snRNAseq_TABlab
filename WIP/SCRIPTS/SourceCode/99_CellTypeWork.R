library(Seurat)
library(SeuratExtend)
library(qs)
library(tidyverse)

# Read previous object if needed.
# integrated <- qs::qread(file.path(qsave_dir, "09_ref_annotated_obj.qs"))

# ID cell clutsteres with FeaturePlot.
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "RNA_cluster_0.7"
Idents(integrated) <- "ref_short"

FeaturePlot3.grid(integrated, features = genes_to_plot)
FeaturePlot3.grid(integrated, features = gene_list)
FeaturePlot(integrated, features = "Bcl11b", label = T) |
    FeaturePlot(integrated, features = "Tle4", label = T) |
    FeaturePlot(integrated, features = "Foxp2", label = T)

DimPlot(integrated, reduction= 'umap', label = T, repel = T) + NoLegend() |
    DimPlot(integrated, reduction= 'umap', label = T, repel = T, group.by = "ref_short") + NoLegend()

cells <- colnames(integrated)[integrated$ref_short == "Ex_Cortex"]
DimPlot2(integrated, cells.highlight = cells)


# Set idents
Idents(integrated) <- integrated$RNA_cluster_0.7

# Check clusters
DimPlot(integrated, reduction = "umap", label = T) + NoLegend()

# Finalized Clusters
cell_type <- c(
    "0" = "ExNeurons",
    "1" = "ExNeurons",
    "2" = "ExNeurons",
    "6" = "ExNeurons",
    "7" = "ExNeurons",
    "10" = "ExNeurons",
    "13" = "ExNeurons",
    "14" = "ExNeurons",
    "15" = "ExNeurons",
    "18" = "ExNeurons",
    "19" = "ExNeurons",
    "22" = "ExNeurons",
    "26" = "ExNeurons",
    "28" = "ExNeurons",
    "32" = "ExNeurons",
    "33" = "ExNeurons",
    "34" = "ExNeurons",
    "36" = "ExNeurons",
    "38" = "ExNeurons",
    "48" = "ExNeurons",

    "4" = "InNeurons",
    "5" = "InNeurons",
    "8" = "InNeurons",
    "9" = "InNeurons",
    "12" = "InNeurons",
    "25" = "InNeurons",
    "30" = "InNeurons",
    "44" = "InNeurons",
    "47" = "InNeurons",
    "49" = "InNeurons",

    "27"  = "NPC",

    "3" = "Astrocytes",
    "11" = "Astrocytes",
    "20" = "Astrocytes",
    "35" = "Astrocytes",
    "39" = "Astrocytes",
    "45" = "Astrocytes",
    "46" = "Astrocytes",
    "50" = "Astrocytes",

    "17" = "Microglia",
    "24" = "Microglia",

    "16" = "OPC",
    "23" = "OPC",

    "42" = "NFOL",

    "21" = "MFOL",
    "43" = "MFOL",

    "37" = "Ependymal",
    "40" = "Ependymal",

    "29" = "Pericytes",
    "31" = "Pericytes",
    "41" = "Pericytes"
)

# Assign the cell types to the metadata based on cluster IDs
integrated@meta.data$CellTypes <- cell_type[as.character(integrated$RNA_cluster_0.7)]

# Save
qs::qsave(integrated, file = file.path(qsave_dir, "10_metadata_edited_obj.qs"))
