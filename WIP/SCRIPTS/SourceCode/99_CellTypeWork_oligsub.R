library(Seurat)
library(SeuratExtend)
library(qs)
library(tidyverse)

# Read previous object if needed.
# integrated_ref_npc <- qs::qread(file.path(qsave_dir, "OligoSubset", "10__ref_npc.qs"))
# integrated_man_npc <- qs::qread(file.path(qsave_dir, "OligoSubset", "10__man_npc.qs"))

# Markers
OligMarkers <- c(
    'Pdgfra', 'Cspg4', 'Olig1', 'Olig2', 'Sox10', 'Sox6',            # OPC markers
    'Ptprz1', 'Lhfpl3', 'Bmp4', 'Neu4', 'Id2', 'Id4', 'Tcf7l2',      # COP markers, Sox10 as well.
    'Enpp6', 'Mog', 'Plp1', 'Prom1', 'Tspan2',                       # NFOL markers 'Ctps1', 'Grp17',   
    'Mobp', 'Mag', 'Mal', 'Opalin', 'Serinc5', 'Cnp', 'Mbp', 'Myrf', # MFOL markers, Mog, Plp1 as well.
    'Cnp', 'Klk6', 'Apod', 'Trf', 'Pmp22'                            # Mature_OL markers, Mobp, Mbp, Mag, Opalin,  Plp1 as well
)

OligMarkers_selected <- c(
    'Pdgfra', 'Olig2', 'Sox6',
    'Ptprz1', 'Neu4', 'Sox10',
    'Enpp6', 'Mog', 'Plp1', 
    'Mag', 'Mbp', 'Cnp'
)

# ID cell clutsteres with FeaturePlot.
DefaultAssay(integrated_man_npc) <- "RNA"
DefaultAssay(integrated_ref_npc) <- "RNA"


# Check clusters
FeaturePlot3.grid(integrated_man_npc, features = OligMarkers_selected) |
    FeaturePlot3.grid(integrated_ref_npc, features = OligMarkers_selected)

FeaturePlot(integrated_man_npc, features = OligMarkers_selected) |
    FeaturePlot(integrated_ref_npc, features = OligMarkers_selected)

DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "Olig_cluster_man_npc_0.4") + NoLegend() |
    DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "ref_short2") + NoLegend()
FeaturePlot(integrated_man_npc, features = OligMarkers_selected)

DimPlot(integrated_ref_npc, reduction= 'umap', label = T, repel = T, group.by = "Olig_cluster_ref_npc_0.4") + NoLegend() |
    DimPlot(integrated_ref_npc, reduction= 'umap', label = T, repel = T, group.by = "ref_short2") + NoLegend()
FeaturePlot(integrated_ref_npc, features = OligMarkers_selected)

#### clusters ####
# Finalized Clusters
cell_type_man <- c(
    "0" = "OPC_1",
    "12" = "OPC_1",
    
    "2" = "OPC_2",
    "5" = "OPC_2",
    "8" = "OPC_2",
    "13" = "OPC_2",
    "14" = "OPC_2",
    "15" = "OPC_2",
    
    "6" = "NPC",
    "9" = "NPC",
    
    "4" = "COP",
    
    "11" = "NFOL",
    
    "7" = "MFOL",
    
    "1" = "Mature_OL",
    "3" = "Mature_OL",
    "10" = "Mature_OL"
)

cell_type_ref <- c(
    "0" = "OPC_1",
    "14" = "OPC_1",
    "15" = "OPC_1",
    
    "2" = "OPC_2",
    "5" = "OPC_2",
    "11" = "OPC_2",
    "18" = "OPC_2",
    "19" = "OPC_2",
    "21" = "OPC_2",
    
    "3" = "NPC",
    "6" = "NPC",
    "9" = "NPC",
    "10" = "NPC",
    
    "4" = "Vascular",
    "16" = "Microglia",
    
    "8" = "COP",
    
    "13" = "NFOL",
    
    "7" = "MFOL",
    
    "1" = "Mature_OL",
    "12" = "Mature_OL",
    "17" = "Mature_OL",
    "20" = "Mature_OL"
)

#### Assign the cell types to the metadata based on cluster IDs ####
integrated_man_npc@meta.data$CellTypes2 <- cell_type_man[as.character(integrated_man_npc$Olig_cluster_man_npc_0.4)]
unique(integrated_man_npc@meta.data$CellTypes2)

integrated_ref_npc@meta.data$CellTypes2 <- cell_type_ref[as.character(integrated_ref_npc$Olig_cluster_ref_npc_0.4)]
unique(integrated_ref_npc@meta.data$CellTypes2)

# Check clusters
DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "Olig_cluster_man_npc_0.4") + NoLegend() |
    DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "CellTypes2") + NoLegend()
Idents(integrated_man_npc) <- "CellTypes2"
FeaturePlot(integrated_man_npc, features = OligMarkers_selected, label = T)

DimPlot(integrated_ref_npc, reduction= 'umap', label = T, repel = T, group.by = "Olig_cluster_ref_npc_0.4") + NoLegend() |
    DimPlot(integrated_ref_npc, reduction= 'umap', label = T, repel = T, group.by = "CellTypes2") + NoLegend()
Idents(integrated_ref_npc) <- "CellTypes2"
FeaturePlot(integrated_ref_npc, features = OligMarkers_selected, label = T)


#### Save ####
qs::qsave(integrated_man_npc, file = file.path(qsave_dir, "OligoSubset", "10__man_npc.qs"))
qs::qsave(integrated_ref_npc, file = file.path(qsave_dir, "OligoSubset", "10__ref_npc.qs"))

#### adjustment- man #### 
# Check
DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "Olig_cluster_man_npc_1") + NoLegend() |
    DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "CellTypes") + NoLegend()
FeaturePlot(integrated_man_npc, features = OligMarkers_selected)

# Finalized Clusters
cell_type_man <- c(
    "0" = "OPC_1",
    "1" = "OPC_1",
    "17" = "OPC_1",
    "18" = "OPC_1",
    "20" = "OPC_1",
    
    
    "16" = "OPC_2",
    "19" = "OPC_2",

    "4" = "OPC_3",
    "6" = "OPC_3",
    "7" = "OPC_3",
    "10" = "OPC_3",
    "21" = "OPC_3",
    "22" = "OPC_3",
    "23" = "OPC_3",
    "24" = "OPC_3",
    
    "5" = "NPC",
    "15" = "NPC",
    
    "9" = "COP/NFOL",
    "12" = "COP/NFOL",
    "11" = "COP/NFOL",
    
    "2" = "MFOL/Mature_OL",
    "8" = "MFOL/Mature_OL",
    "13" = "MFOL/Mature_OL",
    "3" = "MFOL/Mature_OL",
    "14" = "MFOL/Mature_OL"
)

# Assign the cell types to the metadata based on cluster IDs 
integrated_man_npc@meta.data$CellTypes2 <- cell_type_man[as.character(integrated_man_npc$Olig_cluster_man_npc_1.2)]
unique(integrated_man_npc@meta.data$CellTypes2)

# check
DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "Olig_cluster_man_npc_1.2") + NoLegend() |
DimPlot(integrated_man_npc, reduction= 'umap', label = T, repel = T, group.by = "CellTypes2") + NoLegend()
FeaturePlot(integrated_man_npc, features = OligMarkers_selected)

# Save 
qs::qsave(integrated_man_npc, file = file.path(qsave_dir, "OligoSubset", "10__man_npc.qs"))
