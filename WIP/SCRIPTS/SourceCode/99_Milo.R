# library load
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(tidyverse)
library(patchwork)

# Set seed
set.seed(0827)

# Create miloR object
integrated <- qs::qread("~/SLC35A2_Olig2cKO_snRNA/WIP/QSAVE/10_metadata_edited_obj.qs")
integrated_sce <- as.SingleCellExperiment(integrated)
integrated_milo <- Milo(integrated_sce)

qs::qsave(integrated_milo, file = "~/SLC35A2_Olig2cKO_snRNA/WIP/QSAVE/99_miloRobject.qs")

# Construct KNN graph
integrated_milo <- buildGraph(integrated_milo, k = 10, d = 30)

# Defining representative neighbourhoods
# prop: the proportion of cells to randomly sample to start with (usually 0.1 - 0.2 is sufficient)
# k: the k to use for KNN refinement (we recommend using the same k used for KNN graph building)
# d: the number of reduced dimensions to use for KNN refinement (we recommend using the same d used for KNN graph building)
# refined indicated whether you want to use the sampling refinement algorithm, or
# just pick cells at random. The default and recommended way to go is to use refinement.
# The only situation in which you might consider using random instead, is
# if you have batch corrected your data with a graph based correction algorithm, such as BBKNN,
# but the results of DA testing will be suboptimal.

integrated_milo <- makeNhoods(integrated_milo, prop = 0.1, k = 10, d=30, refined = TRUE)

# Check
plotNhoodSizeHist(integrated_milo)

#Counting cells in neighbourhoods
integrated_milo <- countCells(integrated_milo, meta.data = data.frame(colData(integrated_milo)), samples="MouseID")
head(nhoodCounts(integrated_milo))

# Differential abundance testing
design <- data.frame(colData(integrated_milo))[,c("MouseID", "TimeGeno")]
design <- distinct(design)
rownames(design) <- design$Sample

## Reorder rownames to match columns of nhoodCounts(milo)
design <- design[colnames(nhoodCounts(integrated_milo)), , drop=FALSE]
integrated_milo <- calcNhoodDistance(integrated_milo, d=30)

rownames(design) <- design$Sample
da_results <- testNhoods(integrated_milo, design = ~ Condition, design.df = design)
da_results %>%
    arrange(- SpatialFDR) %>%
    head()

# Visualize neighbourhoods displaying DA
integrated_milo <- buildNhoodGraph(integrated_milo)

plotUMAP(integrated_milo) + plotNhoodGraphDA(integrated_milo, da_results, alpha=0.05) +
    plot_layout(guides="collect")

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

integrated_milo <- buildNhoodGraph(integrated_milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(integrated_milo, dimred = "umap", colour_by="stage", text_by = "celltype",
                          text_size = 3, point_size=0.5) +
    guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(integrated_milo, da_results, layout="umap",alpha=0.1)

umap_pl + nh_graph_pl +
    plot_layout(guides="collect")

da_results <- annotateNhoods(integrated_milo, da_results, coldata_col = "celltype")
head(da_results)

ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=50)

plotDAbeeswarm(da_results, group.by = "celltype")

# Finding markers of DA populations
## Add log normalized count to Milo object
integrated_milo <- logNormCounts(integrated_milo)

da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
da_nhood_markers <- findNhoodGroupMarkers(integrated_milo, da_results, subset.row = rownames(integrated_milo)[1:10])

da_nhood_markers <- findNhoodGroupMarkers(integrated_milo, da_results, subset.row = rownames(integrated_milo)[1:10],
                                          aggregate.samples = TRUE, sample_col = "sample")

## Run buildNhoodGraph to store nhood adjacency matrix
integrated_milo <- buildNhoodGraph(integrated_milo)

## Find groups
da_results <- groupNhoods(integrated_milo, da_results, max.lfc.delta = 2)
head(da_results)

plotNhoodGroups(integrated_milo, da_results, layout="umap")

plotDAbeeswarm(da_results, "NhoodGroup")
plotDAbeeswarm(groupNhoods(integrated_milo, da_results, max.lfc.delta = 0.5), group.by = "NhoodGroup") + ggtitle("max LFC delta=0.5")
plotDAbeeswarm(groupNhoods(integrated_milo, da_results, max.lfc.delta = 1) , group.by = "NhoodGroup") + ggtitle("max LFC delta=1")
plotDAbeeswarm(groupNhoods(integrated_milo, da_results, max.lfc.delta = 2)   , group.by = "NhoodGroup") + ggtitle("max LFC delta=2")

da_results <- groupNhoods(integrated_milo, da_results, max.lfc.delta = 3, overlap=5)
plotNhoodGroups(integrated_milo, da_results, layout="umap")
plotDAbeeswarm(da_results, group.by = "NhoodGroup")

## Exclude zero counts genes
keep.rows <- rowSums(logcounts(integrated_milo)) != 0
integrated_milo <- integrated_milo[keep.rows, ]

## Find HVGs
dec <- modelGeneVar(integrated_milo)
hvgs <- getTopHVGs(dec, n=2000)
head(hvgs)

nhood_markers <- findNhoodGroupMarkers(integrated_milo, da_results, subset.row = hvgs,
                                       aggregate.samples = TRUE, sample_col = "sample")

head(nhood_markers)

gr2_markers <- nhood_markers[c("logFC_2", "adj.P.Val_2")]
colnames(gr2_markers) <- c("logFC", "adj.P.Val")

head(gr2_markers[order(gr2_markers$adj.P.Val), ])

nhood_markers <- findNhoodGroupMarkers(integrated_milo, da_results, subset.row = hvgs,
                                       aggregate.samples = TRUE, sample_col = "sample",
                                       subset.groups = c("2")
)

head(nhood_markers)

nhood_markers <- findNhoodGroupMarkers(integrated_milo, da_results, subset.row = hvgs,
                                       subset.nhoods = da_results$NhoodGroup %in% c('11','2'),
                                       aggregate.samples = TRUE, sample_col = "sample")

head(nhood_markers)

ggplot(nhood_markers, aes(logFC_2,-log10(adj.P.Val_2 ))) +
    geom_point(alpha=0.5, size=0.5) +
    geom_hline(yintercept = 2)

markers <- rownames(nhood_markers)[nhood_markers$adj.P.Val_2 < 0.01 & nhood_markers$logFC_2 > 0]

plotNhoodExpressionGroups(integrated_milo, da_results, features=markers,
                          subset.nhoods = da_results$NhoodGroup %in% c('11','2'), scale_to_1 = TRUE,
                          grid.space = "fixed")

dge_9 <- testDiffExp(integrated_milo, da_results, design = ~ stage, meta.data = data.frame(colData(integrated_milo)),
                     subset.row = rownames(integrated_milo)[1:10], subset.nhoods=da_results$NhoodGroup=="9")