# (wrong approach- subsettied obj should be processed from integration) Work only with Olig population
## Subset Olig only and re-process
# ```{r [O] Subset, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# integrated <- qs::qread(file.path(qsave_dir, "10_metadata_edited_obj.qs"))

# Subset based on manual
oligsub_man <- subset(integrated, integrated$CellTypes %in% c("OPC", "NFOL", "MFOL"))

# Subset based on manual
oligsub_ref <- subset(integrated, integrated$ref_short %in% c("OPC", "COP", "NFOL", "MFOL", "Mature_OL"))

# Save
qs::qsave(oligsub_man, file = file.path(qsave_dir, "a1_OligOnlyInit_manual.qs"))
qs::qsave(oligsub_ref, file = file.path(qsave_dir, "a1_OligOnlyInit_refshort.qs"))
# ```   

# ```{r [O] Re-process, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# oligsub_man <- qs::qread(file.path(qsave_dir, "a1_OligOnlyInit_refshort.qs"))
# oligsub_ref <- qs::qread(file.path(qsave_dir, "a1_OligOnlyInit_manual.qs"))

# Normalize/process 
# source('~/SLC35A2_Olig2cKO_snRNA_JY/WIP/SCRIPTS/SourceCode/06_process_oligsub.R')

# Process integrated object
oligsub_man <- process_oligsub(oligsub_man,
                               save = TRUE,
                               suffix = "_manual",
                               qsave_dir = qsave_dir)

oligsub_ref <- process_oligsub(oligsub_ref,
                               save = TRUE,
                               suffix = "_refshort",
                               qsave_dir = qsave_dir)
# ```

## Initial visualization
### Manual cluster
# ```{r [V] DimPlot manual, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# oligsub_man <- qs::qread(file.path(qsave_dir, "a2_olig_processed_manual.qs"))

# set assay and Idents
Seurat::DefaultAssay(oligsub_man) <- "RNA"
Seurat::Idents(oligsub_man) <- "MouseID"

# DimPLot
dp <- Seurat::DimPlot(oligsub_man, reduction = 'umap', group.by = 'TimePoint', alpha = 0.3) |
  Seurat::DimPlot(oligsub_man, reduction = 'umap', group.by = 'GenoType', alpha = 0.3) |
  Seurat::DimPlot(oligsub_man, reduction = 'umap', group.by = "MouseID", alpha = 0.3) |
  Seurat::DimPlot(oligsub_man, reduction = 'umap', group.by = "SampleID", alpha = 0.3)

dp_time <- Seurat::DimPlot(oligsub_man, reduction = 'umap', group.by = 'TimePoint', alpha = 0.3, split.by = "GenoType")

dp_geno <- Seurat::DimPlot(oligsub_man, reduction = 'umap', group.by = 'GenoType', alpha = 0.3, split.by = "TimePoint")

dp_tg <- Seurat::DimPlot(oligsub_man, reduction = 'umap', split.by = 'TimeGeno', alpha = 0.3)

# Save
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_all_manual.png"), dp, width = 18, height = 4, dpi = 300)
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_timepoint_manual.png"), dp_time, width = 9, height = 4, dpi = 300)
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_genotype_manual.png"), dp_geno, width = 9, height = 4, dpi = 300)
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_TimeGeno_manual.png"), dp_tg, width = 18, height = 4, dpi = 300)
#```

### ref_short cluster
#```{r [V] DimPlot refshort, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# oligsub_ref <- qs::qread(file.path(qsave_dir, "a2_olig_processed_refshort.qs"))

# set assay and Idents
Seurat::DefaultAssay(oligsub_ref) <- "RNA"
Seurat::Idents(oligsub_ref) <- "MouseID"

# DimPLot
dp <- Seurat::DimPlot(oligsub_ref, reduction = 'umap', group.by = 'TimePoint', alpha = 0.3) |
  Seurat::DimPlot(oligsub_ref, reduction = 'umap', group.by = 'GenoType', alpha = 0.3) |
  Seurat::DimPlot(oligsub_ref, reduction = 'umap', group.by = "MouseID", alpha = 0.3) |
  Seurat::DimPlot(oligsub_ref, reduction = 'umap', group.by = "SampleID", alpha = 0.3)

dp_time <- Seurat::DimPlot(oligsub_ref, reduction = 'umap', group.by = 'TimePoint', alpha = 0.3, split.by = "GenoType")

dp_geno <- Seurat::DimPlot(oligsub_ref, reduction = 'umap', group.by = 'GenoType', alpha = 0.3, split.by = "TimePoint")

dp_tg <- Seurat::DimPlot(oligsub_ref, reduction = 'umap', split.by = 'TimeGeno', alpha = 0.3)

# Save
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_all_refshort.png"), dp, width = 18, height = 4, dpi = 300)
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_timepoint_refshort.png"), dp_time, width = 9, height = 4, dpi = 300)
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_genotype_refshort.png"), dp_geno, width = 9, height = 4, dpi = 300)
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_TimeGeno_refshort.png"), dp_tg, width = 18, height = 4, dpi = 300)
#```

## Clustering oligsub
#```{r [O] Clustering, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# source('~/SLC35A2_Olig2cKO_snRNA_JY/WIP/SCRIPTS/SourceCode/07_clustering_oligsub.R')

# Read previous object if needed.
# oligsub_man <- qs::qread(file.path(qsave_dir, "a2_olig_processed_manual.qs"))
# oligsub_ref <- qs::qread(file.path(qsave_dir, "a2_olig_processed_refshort.qs"))

oligsub_man <- clustering_oligsub(oligsub_man,
                                  range_start = 0.1,
                                  range_end = 1.2,
                                  range_step = 0.1,
                                  save = TRUE,
                                  suffix = "_manual",
                                  qsave_dir = qsave_dir,
                                  markers_dir = file.path(csv_dir, "Markers"))

oligsub_ref <- clustering_oligsub(oligsub_ref,
                                  range_start = 0.1,
                                  range_end = 1.2,
                                  range_step = 0.1,
                                  save = TRUE,
                                  suffix = "_refshort",
                                  qsave_dir = qsave_dir,
                                  markers_dir = file.path(csv_dir, "Markers"))

#```

### Manual clustering
#```{r [V] Check clusters - manual, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object if needed.
# oligsub_man <- qs::qread(file.path(qsave_dir, "a3_clustered_oligsub_manual.qs"))

# Olig markers
olig_markers <- c(
  'Olig1', 'Olig2', 'Sox10', 'Cd9', 'Smarca4', 'Serinc5',
  'Pdgfra', 'Cspg4', 'Lhfpl3', 'Ptprz1', 'Pmp22',
  'Bmp4', 'Enpp6', 'Neu4', 'Prom1',
  'Cnp',
  'Plp1', 'Mbp', 'Mag', 'Enpp4',
  'Mog', 'Cnp', 'Cd82', 'Opalin', 'Klk6', 'Apod', 'Trf',  'Mobp',
  "Tbx18", 'Vtn', 'Lum', 'Col1a2'
)

# Set assay and Idents
Seurat::DefaultAssay(oligsub_man) <- 'RNA'
Seurat::Idents(oligsub_man) <- oligsub_man$Olig_cluster_0.3

# Glimpse 
Seurat::DimPlot(oligsub_man, reduction = "umap")


# Feature plots
fp_olig <- Seurat::FeaturePlot(oligsub_man, features = olig_markers, label = T)

# Save
ggplot2::ggsave(fp_olig, file = file.path(plot_dir, "FeaturePlot_Olig_markers_manual.png"), width = 16, height= 28, dpi = 300)
#```

### ref_short clustering
#```{r [V] Check clusters - refshort, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object if needed.
# oligsub_ref <- qs::qread(file.path(qsave_dir, "a3_clustered_oligsub_refshort.qs"))

# Olig markers
olig_markers <- c(
  'Olig1', 'Olig2', 'Sox10', 'Cd9', 'Smarca4', 'Serinc5',
  'Pdgfra', 'Cspg4', 'Lhfpl3', 'Ptprz1', 'Pmp22',
  'Bmp4', 'Enpp6', 'Neu4', 'Prom1',
  'Cnp',
  'Plp1', 'Mbp', 'Mag', 'Enpp4',
  'Mog', 'Cnp', 'Cd82', 'Opalin', 'Klk6', 'Apod', 'Trf',  'Mobp',
  "Tbx18", 'Vtn', 'Lum', 'Col1a2')

# Set assay and Idents
DefaultAssay(oligsub_ref) <- 'RNA'
Idents(oligsub_ref) <- oligsub_ref$Olig_cluster_0.3

# Glimpse 
Seurat::DimPlot(oligsub_ref, reduction = "umap")

# Feature plots
fp_olig <- Seurat::FeaturePlot(oligsub_ref, features = olig_markers, label = T)

# Save
ggplot2::ggsave(fp_olig, file = file.path(plot_dir, "FeaturePlot_Olig_markers_refshort.png"), width = 16, height= 28, dpi = 300)
#```

#```{r [V] Visualization from paper, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# https://www.sciencedirect.com/science/article/pii/S0959438817301150
# https://www.sciencedirect.com/science/article/pii/S2213671124000778?via%3Dihub
# https://www.sciencedirect.com/science/article/pii/S1534580718305586?via%3Dihub
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5221728/pdf/emss-70326.pdf
# https://onlinelibrary.wiley.com/doi/10.1111/epi.18413
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4691794/

# Set assay and Idents
Seurat::DefaultAssay(oligsub_man) <- "RNA"
Seurat::DefaultAssay(oligsub_ref) <- "RNA"

Seurat::Idents(oligsub_man) <- oligsub_man$Olig_cluster_0.3
Seurat::Idents(oligsub_ref) <- oligsub_ref$Olig_cluster_0.3

# Markers list
markers <- list(
  OligMarkers = c('Olig1', 'Olig2', 'Sox10', 'Cd9', 'Smarca4', 'Serinc5',
                  'Pdgfra', 'Cspg4', 'Lhfpl3', 'Ptprz1', 'Pmp22',
                  'Bmp4', 'Enpp6', 'Neu4', 'Prom1',
                  'Cnp',
                  'Plp1', 'Mbp', 'Mag', 'Enpp4',
                  'Mog', 'Cnp', 'Cd82', 'Opalin', 'Klk6', 'Apod', 'Trf',  'Mobp',
                  "Tbx18", 'Vtn', 'Lum', 'Col1a2'),
  paper1 = c('Ascl1', 'Nfix', 'Dll3', 'Slc38a1', 'Rplp0', 'Rpl4', 'Eef2', 'Npas1', 'Hes', 'Epas1',
             'Pnlip', 'Pcp4', 'Ptprn', 'Nrarp', 'Clu', 'Gcp5', 'Ttyh1'),
  paper2 = c("Top2a", "Fyn", "Etv5", "Ednrb"),
  paper3 = c("Adora1", "Brinp3", "Resp18", "Matn4", "Ccnb1", "Bche", "Cpm", "Elovl7", "Parvb", "Rspo2",
             "Wnt8b", "Clic6", "Trpm3", "Pax2", "Arx", "Lhx6", "Col6a3", "Ddr2", "Abcc9", "Kcnj8",
             "Gjb6", "Ranbp3l", "Itih2", "Il33"),
  paper4 = c("Nr2e1", "Foxg1", "Neurod6", "Otx2", "Nkx2-1", "Eomes", "Hoxa13", "Hoxc11", "Hoxd9",
             "Hoxa11", "Hoxc6", "Hoxb8", "Hoxa3", "Hoxb6", "Hoxc8", "Hoxc9", "Hoxc10", "Hoxa5", "Hoxa6",
             "Hoxd10", "Hoxb2", "Hoxd11", "Hoxa9", "Hoxa7", "Hoxb7", "Hoxa10", "Hoxd9", "Hoxb9", "Pax8",
             "Tlx3", "Nhlh2"),
  paper5 = c("Gli1", "Ebf3", "Scrt2", "Mecom", "Lin28b", "Plag1", "Zic4", "Tcfap2b", "Osr1", "Tead2",
             "Gli3", "Prdm6", "Bnc2", "Gli2", "Smad6", "Hoxb9", "Pax3", "Neurog2", "Lhx9", "Isl1",
             "Neurod6", "Neurod1", "Neurod2", "Dlx1", "Satb2", "Sox8", "Nkx2-2", "Dbx2", "Barx2",
             "Npas1", "Nkx6-2", "Cux1", "Cux2", "Bcl11b"),
  IEGs = c("Jun", "Fos", "Egr1", "Junb")
)

# FeaturePlot
Seurat::FeaturePlot(oligsub_man, features = markers[['OligMarkers']])
Seurat::FeaturePlot(oligsub_ref, features = markers[['OligMarkers']])

Seurat::FeaturePlot(oligsub_man, features = markers[['paper1']])
Seurat::FeaturePlot(oligsub_ref, features = markers[['paper1']])

Seurat::FeaturePlot(oligsub_man, features = markers[['paper2']])
Seurat::FeaturePlot(oligsub_ref, features = markers[['paper2']])

Seurat::FeaturePlot(oligsub_man, features = markers[['paper3']])
Seurat::FeaturePlot(oligsub_ref, features = markers[['paper3']])

Seurat::FeaturePlot(oligsub_man, features = markers[['paper4']])
Seurat::FeaturePlot(oligsub_ref, features = markers[['paper4']])

Seurat::FeaturePlot(oligsub_man, features = markers[['paper5']])
Seurat::FeaturePlot(oligsub_ref, features = markers[['paper5']])

Seurat::FeaturePlot(oligsub_man, features = markers[['IEGs']])
Seurat::FeaturePlot(oligsub_ref, features = markers[['IEGs']])
# ```

# ```{r [O] Name cluster, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Set Idents
Seurat::Idents(oligsub_man) <- oligsub_man$ref_short

# Glimpse
Seurat::DimPlot(oligsub_man, reduction = 'umap', label = T)

# Map
map <- c(
  "0" = "OPC_1", "3" = "OPC_1",
  "5" = "COP",
  "9" = "NFOL",
  "1" = "MFOL",
  "6" = "Mature_OL", "4" = "Mature_OL",
  "2" = "OPC_2", "7" = "OPC_2", "8" = "OPC_2","10" = "OPC_2"
)

# Add to meta data
oligsub_man@meta.data$OligSub <- map[as.character(oligsub_man$Olig_cluster_0.3)]

# Set Idents
Seurat::Idents(oligsub_man) <- oligsub_man$OligSub

# Glimpse
Seurat::DimPlot(oligsub_man, reduction = 'umap', label = T)

# FeaturePlot
SeuratExtend::FeaturePlot3.grid(oligsub_man, features = c('Pdgfra', 'Cspg4', 'Sox10', 'Enpp6', 'Mog', 'Plp1', 'Mag', 'Mbp', 'Opalin', 'Mobp', 'Serinc5', 'Apod'))

Seurat::FeaturePlot(oligsub_man, features = c('Pdgfra', 'Enpp6', 'Plp1', 'Mag', 'Mbp', 'Apod'), label = T)

# Re-Subset based on manual
oligsub_man2 <- subset(oligsub_man, oligsub_man$OligSub %in% c("OPC_1", "COP", "NFOL", "MFOL", "Mature_OL"))

# Glimpse
Seurat::DimPlot(oligsub_man2, reduction = 'umap', label = T, split.by = "TimeGeno", group.by = "OligSub")
Seurat::FeaturePlot(oligsub_man2, features = c('Pdgfra', 'Enpp6', 'Plp1', 'Mag', 'Mbp', 'Apod'), label = T)

###############################################################
# Set Idents
Seurat::Idents(oligsub_ref) <- oligsub_ref$ref_short

# Glimpse
Seurat::DimPlot(oligsub_ref, reduction = 'umap', label = T)

# FeaturePlot
SeuratExtend::FeaturePlot3.grid(oligsub_ref, features = c('Pdgfra', 'Cspg4', 'Sox10', 'Enpp6', 'Mog', 'Plp1', 'Mag', 'Mbp', 'Opalin', 'Mobp', 'Serinc5', 'Apod'))

Seurat::FeaturePlot(oligsub_ref, features = c('Pdgfra', 'Enpp6', 'Plp1', 'Mag', 'Mbp', 'Apod'), label = T)


# ```


# (wrong approach- subsettied obj should be processed from integration) Work only with Olig population - Doublet ver.
## Subset Olig only and re-process
# ```{r [O] Subset, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# integrated <- qs::qread(file.path(qsave_dir, "11_dbl_added.qs"))

# Subset based on manual
oligsub_man <- subset(integrated, integrated$CellTypes %in% c("OPC", "NFOL", "MFOL"))

# Subset based on manual
oligsub_ref <- subset(integrated, integrated$ref_short %in% c("OPC", "COP", "NFOL", "MFOL", "Mature_OL"))

# Save
qs::qsave(oligsub_man, file = file.path(qsave_dir, "b1_OligOnlyInit_manual.qs"))
qs::qsave(oligsub_ref, file = file.path(qsave_dir, "b1_OligOnlyInit_refshort.qs"))
# ```   

# ```{r [O] Re-process, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# oligsub_man <- qs::qread(file.path(qsave_dir, "b1_OligOnlyInit_refshort.qs"))
# oligsub_ref <- qs::qread(file.path(qsave_dir, "b1_OligOnlyInit_manual.qs"))

# Normalize/process 
# check save file name in source code, if you put save = T)
# source("~/SLC35A2_Olig2cKO_snRNA_JY/WIP/SCRIPTS/SourceCode/06_process_oligsub.R")

# Process integrated object
oligsub_man <- process_oligsub(oligsub_man, save = F)

oligsub_ref <- process_oligsub(oligsub_ref, save = F)

# Save
qs::qsave(oligsub_man, file = file.path(qsave_dir, "b2_olig_processed_manual.qs"))
qs::qsave(oligsub_ref, file = file.path(qsave_dir, "b2_olig_processed_refshort.qs"))
# ```

## Initial visualization
### Manual cluster
# ```{r [V] DimPlot manual, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# oligsub_man <- qs::qread(file.path(qsave_dir, "b2_olig_processed_manual.qs"))

# set assay and Idents
Seurat::DefaultAssay(oligsub_man) <- "RNA"
Seurat::Idents(oligsub_man) <- "MouseID"

# DimPLot
dp_dbl <- Seurat::DimPlot(oligsub_man, reduction = 'umap', split.by = 'Doublet_Status', alpha = 0.3)

# Save
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_doublet_manual.png"), dp_dbl, width = 18, height = 4, dpi = 300)
# ```

### ref_short cluster
# ```{r [V] DimPlot refshort, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# Read previous object
# oligsub_ref <- qs::qread(file.path(qsave_dir, "b2_olig_processed_refshort.qs"))

# set assay and Idents
Seurat::DefaultAssay(oligsub_ref) <- "RNA"
Seurat::Idents(oligsub_ref) <- "MouseID"

# DimPLot
dp_dbl <- Seurat::DimPlot(oligsub_ref, reduction = 'umap', split.by = 'Doublet_Status', alpha = 0.3)

# Save
ggplot2::ggsave(filename = file.path(plot_dir, "DimPlot_umap_oligsub_doublet_refshort.png"), dp_dbl, width = 18, height = 4, dpi = 300)
# ```

## Clustering oligsub
# ```{r [O] Clustering, echo=TRUE, message=TRUE, warning=TRUE, paged.print=TRUE}
# check save file name in source code, if you put save = T)
# source('~/SLC35A2_Olig2cKO_snRNA_JY/WIP/SCRIPTS/SourceCode/07_clustering_oligsub.R')

# Read previous object if needed.
# oligsub_man <- qs::qread(file.path(qsave_dir, "b2_olig_processed_manual.qs"))
# oligsub_ref <- qs::qread(file.path(qsave_dir, "b2_olig_processed_refshort.qs"))

oligsub_man <- clustering_oligsub(oligsub_man,
                                  range_start = 0.1,
                                  range_end = 1.2,
                                  range_step = 0.1,
                                  save = F)

oligsub_ref <- clustering_oligsub(oligsub_ref,
                                  range_start = 0.1,
                                  range_end = 1.2,
                                  range_step = 0.1,
                                  save = F)

# Save
qs::qsave(oligsub_man, file = file.path(qsave_dir, "b3_clustered_oligsub_manual.qs"))
qs::qsave(oligsub_ref, file = file.path(qsave_dir, "b3_clustered_oligsub_refshort.qs"))
# ```

