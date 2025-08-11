library(Seurat)
library(sctransform)
library(ggplot2)
library(clustree)
library(SeuratDisk)
dyn.load("/igm/apps/hdf5/hdf5-1.12.1/lib/libhdf5_hl.so.200")
library(hdf5r)


mouse.loom <- Connect(filename = "~/PTZ_ATAC_scRNA_072024/File/l5_all.loom")
mouseatlas <- as.Seurat(mouse.loom)
saveRDS(mouseatlas, file = "~/PTZ_ATAC_scRNA_072024/WIP/File/mouseatlas.rds")

Idents(mouseatlas) <- mouseatlas$ClusterName
# Idents(mouseatlas) <- mouseatlas$TaxonomyRank1
# Idents(mouseatlas) <- mouseatlas$TaxonomyRank2
# Idents(mouseatlas) <- mouseatlas$TaxonomyRank3
# Idents(mouseatlas) <- mouseatlas$TaxonomyRank4
# Idents(mouseatlas) <- mouseatlas$Taxonomy_group

# hippocampus/ctx
mouseatlas.cut <- subset(mouseatlas, idents = c("TEGLU1","TEGLU3","TEGLU20","TEGLU2",
                                                "TEGLU13","TEGLU14","TEGLU5","TEGLU21",
                                                "TEGLU24","TEGLU23","DGGRC1","DGGRC2",
                                                "DGNBL2","DGNBL1","SZNBL","TEINH17",
                                                "TEINH18","TEINH19","TEINH21","TEINH16",
                                                "TEINH15","TEINH14","TEINH20","TEINH13",
                                                "TEINH12","TEINH9","TEINH10","TEINH11",
                                                "TEINH4","TEINH5","TEINH8","TEINH7",
                                                "TEINH6","CR","RGDG","ACTE1","ACTE2",
                                                "PER3","PER1","VECA","VECC","VECV",
                                                "PVM1","PVM2","MGL3","MGL2","MGL1",
                                                "COP1","NFOL1","MFOL2","MOL1","EPEN",
                                                "OPC","VLMC1","VSMCA"))

saveRDS(mouseatlas.cut, file = "~/PTZ_ATAC_scRNA_072024/WIP/File/mouseatlas_cut_hippocampus.rds")

Idents(mouseatlas.cut) <- mouseatlas.cut$ClusterName
# Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank1
# Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank2
# Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank3
# Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank4
# Idents(mouseatlas.cut) <- mouseatlas.cut$Taxonomy_group
# Idents(mouseatlas.cut) <- mouseatlas.cut$Description

levels(Idents(mouseatlas.cut))

# process
mouseatlas.cut <- SCTransform(mouseatlas.cut, method = "glmGamPoi",verbose = T)
DefaultAssay(mouseatlas.cut) <- "SCT"
mouseatlas.cut <- RunPCA(mouseatlas.cut)

# Let's find neighbors, clusters, markers and identify our clusters

mouseatlas.cut <- FindNeighbors(mouseatlas.cut, dims = 1:30, verbose = T)
mouseatlas.cut <- RunUMAP(mouseatlas.cut, dims = 1:30, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.2, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.3, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.4, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.5, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.6, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.7, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.8, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 0.9, verbose = T)
mouseatlas.cut <- FindClusters(mouseatlas.cut,resolution = 1, verbose = T)

clustree(mouseatlas.cut)

mouseatlas.cut$seurat_clusters <- mouseatlas.cut$SCT_snn_res.0.8
mouseatlas.cut@active.ident <- mouseatlas.cut$TaxonomyRank4

Idents(mouseatlas.cut) <- mouseatlas.cut$ClusterName
Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank1
Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank2
Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank3
Idents(mouseatlas.cut) <- mouseatlas.cut$TaxonomyRank4
Idents(mouseatlas.cut) <- mouseatlas.cut$Taxonomy_group
Idents(mouseatlas.cut) <- mouseatlas.cut$Description

unique(mouseatlas.cut$TaxonomyRank4)

p1 <- DimPlot(mouseatlas.cut, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
# ggsave("/igm/home/hxy008/PTZ_ATAC_scRNA_072024/WIP/plot/mouseatlas_cut_umap.png", p1, width = 12, height = 10)
saveRDS(mouseatlas.cut, file = "~/SLC35A2_Olig2cKO_snRNA/WIP/mouseatlas_cut_processed.rds")