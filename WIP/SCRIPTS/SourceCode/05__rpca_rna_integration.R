# read this function with
# source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/05__rpca_rna_integration.R')

#' Integrate with RPCA Seurat objects
#'
#' This function performs RPCA-based integration of a list of Seurat objects.
#' RunPCA: npcs = 48 for RunPCA bc smallest cell count is 49
#' IntegrateDate: lowered k.weight to 32.
#' @param sobj_list A named list of Seurat objects.
#' @param save Logical. Whether to save the integrated object (default: FALSE).
#' @param save_dir Directory to save the integrated object if `save = TRUE`.
#'
#' @return Integrated Seurat object.
#' @export
rpca_integration <- function(
        sobj_list,
        npc = 50,
        save = TRUE,
        suffix = NULL,
        save_dir = NULL) {
    # Save check
    if (save) {if (is.null(save_dir)) {stop("You must provide save_dir when save = TRUE.")}
        if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
    }
    
    # Normalize and Findvariablefeatures
    sobj_list <- lapply(sobj_list, function(x) {
        x <- Seurat::NormalizeData(x)
        x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    
    # Select features for integration
    # Identifies common features (genes) across datasets that are informative for integration.
    # These features are used for anchor identification and subsequent data integration.
    features <- Seurat::SelectIntegrationFeatures(object.list = sobj_list)
    
    
    # Prepare for integration
    # Prepares the datasets for integration by ensuring that the same SCT model is applied to the selected features across datasets.
    sobj_list <- lapply(sobj_list, function(x) {
        x <- Seurat::ScaleData(x, features = features, verbose = TRUE)
        x <- Seurat::RunPCA(x, features = features, verbose = TRUE, npcs = npc)
    })
    
    # Find integration anchors
    # Identifies anchors (shared cells or cell populations) across datasets based on the selected features.
    anchors <- Seurat::FindIntegrationAnchors(sobj_list, anchor.features = features, reduction = 'rpca')
    
    # Find integration anchors
    # Identifies anchors (shared cells or cell populations) across datasets based on the selected features.
    integrated_rna <- IntegrateData(anchorset = anchors, k.weight = 32)
    
    # Save object
    if (save) {qs::qsave(integrated_rna, file = file.path(save_dir, paste0("05__integrated_rna", suffix, ".qs")))
        message('Saved integrated rna object to ', file.path(save_dir, paste0("05__integrated_rna", suffix, ".qs")))
    }
    return(integrated_rna)
}