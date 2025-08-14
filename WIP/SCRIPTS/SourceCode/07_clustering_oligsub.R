# read this function with
# source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/07__clustering_oligsub.R')

#' Process integrated Seurat object (PCA, clustering, markers, UMAP/TSNE)
#'
#' @param integrated_obj Integrated Seurat object
#' @param assay Name of assay to use (default: "integrated")
#' @param graph_name Character vector with kNN and SNN graph names (default: c("RNA_nn", "RNA_snn"))
#' @param cluster_name Prefix for clustering results
#' @param umap_name Name to assign to UMAP embedding
#' @param range_start Start of resolution range for clustering
#' @param range_end End of resolution range for clustering
#' @param range_step Step size for resolution
#' @param save Logical, whether to save outputs
#' @param qsave_dir Directory to save processed object
#' @param markers_dir Directory to save marker CSVs
#' @return Processed Seurat object
#' @export
clustering_oligsub <- function(
        integrated_rna,
        range_start = 0.4,
        range_end = 1.2,
        range_step = 0.1,
        save = TRUE,
        suffix = NULL,
        qsave_dir = NULL,
        markers_dir = NULL) {
    # Save check
    if (save) {if (is.null(qsave_dir)) {stop("You must provide 'qsave_dir' when save = TRUE.")}
        if (is.null(markers_dir)) {stop("You must provide 'markers_dir' when save = TRUE.")}
        if (!dir.exists(qsave_dir)) {dir.create(qsave_dir, recursive = TRUE)}
        if (!dir.exists(markers_dir)) {dir.create(markers_dir, recursive = TRUE)}
    }

    # Process
    Seurat::DefaultAssay(integrated_rna) <- 'integrated'

    # Clustering
    cluster_range = seq(range_start, range_end, by = range_step)
    for (i in cluster_range) {
        integrated_rna <- Seurat::FindClusters(
            integrated_rna,
            graph.name = 'snn',
            resolution = i,
            cluster.name = paste0("Olig_cluster", suffix, "_", i))
    }

    # Find markers
    # The "SCT" assay contains normalized and variance-stabilized data suitable for identifying differentially expressed genes between clusters.
    # The "integrated" assay is batch-corrected for alignment and not designed for marker detection as it may suppress biological variation in favor of alignment.
    # only.pos = TRUE: Returns only positive markers (regions more accessible in the cluster).
    # min.pct: Minimum fraction of cells expressing the feature in the cluster.
    # logfc.threshold: Minimum log fold-change required to call a feature significant.
    markers <- list()
    for (i in cluster_range) {
        cluster <- paste0('Olig_cluster', suffix, "_", i)
        Seurat::Idents(integrated_rna) <- integrated_rna[[cluster]][,1]
        markers[[cluster]] <- Seurat::FindAllMarkers(
            integrated_rna,
            assay = 'integrated',
            slot = 'data',
            only.pos = TRUE,
            min.pct = 0.25,
            logfc.threshold = 0.25)
    }

    # Save object
    if (save) {qs::qsave(integrated_rna, file = file.path(qsave_dir, paste0("08_clustered", suffix, ".qs")))
        message('Saved processed integraed rna objec to ', file.path(qsave_dir, paste0("08_clustered", suffix, ".qs")))

        for (name in names(markers)) {
            write.csv(markers[[name]], file = file.path(markers_dir, paste0(name, "_markers", suffix, ".csv")))
            message('Saved processed markers to ', file.path(markers_dir, paste0(name, "_markers", suffix, ".csv")))
        }
    }
    return(integrated_rna)
}
