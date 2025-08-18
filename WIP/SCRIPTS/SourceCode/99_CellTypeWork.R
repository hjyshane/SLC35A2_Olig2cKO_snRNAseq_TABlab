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

DimPlot(integrated, reduction= 'umap', label = T, repel = T, group.by = "RNA_cluster_0.7") + NoLegend() |
    DimPlot(integrated, reduction= 'umap', label = T, repel = T, group.by = "ref_short") + NoLegend()


# Set idents
Idents(integrated) <- integrated$RNA_cluster_0.7

# Check clusters
DimPlot(integrated, reduction = "umap", label = T) + NoLegend()

# Finalized Clusters
cell_type <- c(
    "0" = "ExN",
    "1" = "ExN",
    "2" = "ExN",
    "6" = "ExN",
    "7" = "ExN",
    "10" = "ExN",
    "13" = "ExN",
    "14" = "ExN",
    "15" = "ExN",
    "18" = "ExN",
    "19" = "ExN",
    "22" = "ExN",
    "26" = "ExN",
    "28" = "ExN",
    "32" = "ExN",
    "33" = "ExN",
    "34" = "ExN",
    "36" = "ExN",
    "38" = "ExN",
    "46" = "ExN",
    "48" = "ExN",

    "4" = "InhN",
    "5" = "InhN",
    "8" = "InhN",
    "9" = "InhN",
    "12" = "InhN",
    "25" = "InhN",
    "30" = "InhN",
    "44" = "InhN",
    "47" = "InhN",
    "49" = "InhN",

    "27"  = "NPC",

    "3" = "Astrocytes",
    "11" = "Astrocytes",
    "20" = "Astrocytes",
    "35" = "Astrocytes",
    "39" = "Astrocytes",
    "45" = "Astrocytes",
    "50" = "Astrocytes",

    "17" = "Microglia",
    "24" = "Microglia",

    "16" = "OPC_1",
    
    "23" = "OPC_2",

    "42" = "COP/NFOL",

    "21" = "MFOL",
    
    "43" = "Mature_OL",

    "37" = "Ependymal",
    
    "40" = "Ependymal",

    "29" = "Vascular",
    
    "31" = "Pericytes",
    "41" = "Pericytes"
)

# Assign the cell types to the metadata based on cluster IDs
integrated@meta.data$CellTypes <- cell_type[as.character(integrated$RNA_cluster_0.7)]
unique(integrated@meta.data$CellTypes)

# Save
qs::qsave(integrated, file = file.path(qsave_dir, "10_metadata_edited_obj.qs"))
