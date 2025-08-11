library(Seurat)
library(SoupX)
library(ggplot2)
library(knitr)

# Make the soup
toc = Seurat::Read10X(file.path("/igm/home/afh005/Bedrosian/Collaborations/switzerland_071025/P1_250116G/P1_250116G/outs/per_sample_outs/F-0M_28/count/", "sample_filtered_feature_bc_matrix"))

tod = Seurat::Read10X(file.path("/igm/home/afh005/Bedrosian/Collaborations/switzerland_071025/P1_250116G/P1_250116G/outs/per_sample_outs/F-0M_28/count/", "sample_raw_feature_bc_matrix"))

obj  <- CreateSeuratObject(counts = toc)

sc = SoupChannel(tod, toc)

obj    <- NormalizeData(obj, verbose = F)
obj    <- FindVariableFeatures(obj, verbose = F)
obj    <- ScaleData(obj, verbose = F)
obj    <- RunPCA(obj, verbose = F)
obj    <- RunUMAP(obj, dims = 1:30, verbose = F)
obj    <- FindNeighbors(obj, dims = 1:30, verbose = F)
obj    <- FindClusters(obj, resolution = 0.1, verbose = T)

meta    <- obj@meta.data
umap    <- obj@reductions$umap@cell.embeddings
soup.channel  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
head(meta)

soup.channel  <- autoEstCont(soup.channel)
head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)
plotMarkerDistribution(soup.channel)

adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)

obj <- CreateSeuratObject(adj.matrix)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rpl")

VlnPlot(obj, features = c("percent.mt", "percent.ribo"), pt.size = 0, y.max = 10)

obj <- subset(obj, percent.ribo < 2)

saveRDS(obj, file = "/igm/home/tab013/Fabio_finaldata/F0M_28_SoupX.rds")

obj    <- NormalizeData(obj, verbose = F)
obj    <- FindVariableFeatures(obj, verbose = F)
obj    <- ScaleData(obj, verbose = F)
obj    <- RunPCA(obj, verbose = F)
obj    <- RunUMAP(obj, dims = 1:30, verbose = F)
obj    <- RunTSNE(obj, dims = 1:30, verbose = F)
obj    <- FindNeighbors(obj, dims = 1:30, verbose = F)
obj    <- FindClusters(obj, resolution = 0.1, verbose = T)
DimPlot(obj, reduction = "tsne")
FeaturePlot(obj, reduction = "tsne", features = c("Prox1", "Bcl11b", "Satb2", "Pdgfra", "Mki67", "Gad1"))
FeaturePlot(obj, reduction = "tsne", features = c("Tle4", "Ar", "Esr1", "Sox2", "Aqp4", "P2ry12"))
markers <- FindAllMarkers(obj, only.pos = T)
