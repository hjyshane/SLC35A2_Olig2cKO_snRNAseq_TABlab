# read this function with
# source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/06__process_integrated.R')

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
process_oligsub <- function(
    integrated_rna,
    save = TRUE,
    suffix = NULL,
    qsave_dir = NULL) {
  # Save check
  if (save) {if (is.null(qsave_dir)) {stop("You must provide 'qsave_dir' when save = TRUE.")}
    if (!dir.exists(qsave_dir)) {dir.create(qsave_dir, recursive = TRUE)}
  }

  # Process
  Seurat::DefaultAssay(integrated_rna) <- 'integrated'

  integrated_rna <- integrated_rna %>%
    # Scale data
    Seurat::ScaleData(verbose = TRUE) %>%

    # PCA
    Seurat::RunPCA(
      assay = 'integrated',
      verbose = T) %>%
    #	Identify cell neighborhoods and perform clustering:
    #	Compute the k-nearest neighbor graph based on the PCA results.
    #	kNN graph (RNA_nn): Represents direct neighbors of each cell.	A k-nearest neighbor graph
    #	SNN graph (RNA_snn): Shared nearest neighbor graph, used for clustering. derived from the kNN graph
    Seurat::FindNeighbors(
      dims = 1:30,
      reduction = 'pca',
      graph.name = c('nn', 'snn'),
      verbose = T) %>%
    # Dimensional recudtion for visual
    Seurat::RunUMAP(
      reduction = 'pca',
      reduction.name = 'umap',
      dims = 1:30,
      verbose = T) %>%
    Seurat::RunTSNE(
      reduction = 'pca',
      dims = 1:30,
      verbose = T)

  # Save object
  if (save) {
    qs::qsave(integrated_rna, file = file.path(qsave_dir, paste0("12_olig_processed", suffix, ".qs")))
    message('Saved processed integraed rna objec to ', file.path(qsave_dir, paste0("12_olig_processed", suffix, ".qs")))
    }
  return(integrated_rna)
}

