# read this function with
# source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/03_qc_threshold_test.R')

#' Function for calculating how many cells pass the qc thresholds
#' @param metadata
#' @param min_count
#' @param max_count
#' @param min_feature
#' @param max_feature
#' @param max_mt
#'
#' @return
#'
#'@examples
#'
#' @export
#'
qc_threshold <-
    function(
        seurat_list,
        min_count = 0,
        max_count = 100000,
        min_feature = 0,
        max_mt = 5) {
        # Create empty table
        table <- data.frame(
            MouseID = character(),
            TotalCells = numeric(),
            PassedCells = numeric(),
            PassedRate = numeric())

        # Loop through the list
        for (sample in  names(seurat_list)){
            # MouseID
            ID <- sample

            # Get metadata
            metadata <- seurat_list[[sample]]@meta.data

            # Get Total number of cells
            TotalCells <- nrow(metadata)

            # Proportion that pass the thresholds
            PassedCells <-
                sum(
                    metadata$nCount_RNA >= min_count &
                    metadata$nCount_RNA <= max_count &
                    metadata$nFeature_RNA >= min_feature &
                    metadata$percent_mt_rna <= max_mt)

            # Ratio
            PassedRate <- PassedCells/TotalCells

            # Put in the table
            dt <- data.frame(
                MouseID = ID,
                TotalCells = TotalCells,
                PassedCells = PassedCells,
                PassedRate = PassedRate)

            table <- bind_rows(table, dt)
        }
        return(table)
    }


