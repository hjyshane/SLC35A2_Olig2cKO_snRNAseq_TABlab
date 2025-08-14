# read this function with
# source('~/SLC35A2_Olig2cKO_snRNA_JY/WIP/SCRIPTS/SourceCode/99_doublet_process_rna.R')
#' Run DoubletFinder analysis on a list of Seurat objects
#'
#' This function performs doublet detection using DoubletFinder on each Seurat object
#' in a list. It optimizes the pK parameter through parameter sweeping, calculates
#' expected doublet rates based on cell count, and identifies potential doublets
#' using SCTransform-compatible settings.
#'
#' @param processed_list A list of Seurat objects to process for doublet detection.
#'   Each object should have SCTransform normalization already applied.
#'
#' @return A list of Seurat objects with doublet classifications added to metadata.
#'   New columns will include pANN scores and DF.classifications indicating
#'   "Singlet" or "Doublet" status for each cell.
#' @param save Logical; whether to save results (default = TRUE).
#' @param save_dir Directory to save volcano plots.
#' @details
#' The function performs the following steps for each Seurat object:
#' \itemize{
#'   \item Parameter sweep to find optimal pK value using BCmvn optimization
#'   \item Calculates expected doublet rate based on cell count (0.8% per 1000 cells)
#'   \item Runs DoubletFinder with optimized parameters
#'   \item Uses PCs 1:30, pN=0.25, and SCTransform-compatible settings
#' }
#'
#' @examples
#' \dontrun{
#' # Process a list of Seurat objects for doublet detection
#' processed_list <- find_dbl(processed_list)
#' 
#' # Check doublet classification results
#' lapply(processed_list, function(x) {
#'   doublet_col <- grep("DF.classifications", colnames(x@meta.data), value = TRUE)
#'   if(length(doublet_col) > 0) table(x@meta.data[[doublet_col[1]]])
#' })
#' }
#'
#' @seealso \code{\link[DoubletFinder]{doubletFinder}}, \code{\link[DoubletFinder]{paramSweep}}
#'
#' @export

find_dbl <- function(processed_list,
                     save = FALSE,
                     save_dir = NULL){
  # Save check
  if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
    if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
  }
  for (i in seq_along(processed_list)) {
    
    # Set sequential processing to avoid cluster issues
    library(future)
    plan("sequential")
    
    sample_name <- names(processed_list)[i]  # Get sample name for messaging
    message("Processing sample: ", ifelse(is.null(sample_name), paste0("Sample_", i), sample_name))
    
    test_sample <- processed_list[[i]]
    
    
    # Set the correct assay for DoubletFinder
    DefaultAssay(test_sample) <- "RNA"  # DoubletFinder works on raw counts
    
    # Parameter sweep across pK values
    # Run parameter sweep
    # pN = 0.25 means we create artificial doublets equal to 25% of real cells
    # pK range: tests from 0.005 to 0.3 in increments
    sweep_res <- paramSweep(
      seu = test_sample,           # Seurat object
      PCs = 1:30,                  # Principal components to use
      sct = TRUE,                  # SCTransform was used
      num.cores = 1             
    )
    
    # Summarize sweep results
    sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
    
    # Find optimal pK value
    bcmvn <- find.pK(sweep_stats)
    optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    
    message("Optimal pK found: ", optimal_pK)
    
    # Calculate expected doublet rate (auto-estimation)
    cell_count <- ncol(test_sample)
    expected_doublet_rate <- (cell_count / 1000) * 0.008
    nExp_poi <- round(expected_doublet_rate * cell_count)
    
    # Run DoubletFinder
    test_sample <- doubletFinder(
      seu = test_sample,
      PCs = 1:30,
      pN = 0.25,
      pK = optimal_pK,
      nExp = nExp_poi,
      reuse.pANN = FALSE,
      sct = TRUE,
      annotations = 
    )
    
    # rename column
    df_col <- grep("DF.classifications", colnames(test_sample@meta.data), value = TRUE)
    if (length(df_col) > 0) {
      # Rename to standardized name
      colnames(test_sample@meta.data)[colnames(test_sample@meta.data) == df_col[1]] <- "Doublet_Status"
      message("Renamed ", df_col[1], " to Doublet_Status")
    }
    
    pann_cols <- grep("pANN", colnames(test_sample@meta.data), value = TRUE)
    if (length(pann_cols) > 0) {
      colnames(test_sample@meta.data)[colnames(test_sample@meta.data) == pann_cols[1]] <- "pANN"
      message("Renamed ", pann_cols[1], " to pANN")
    }
    
    processed_list[[i]] <- test_sample
    
  }
  
  # Save if save is enabled
  if (save){
    qs::qsave(processed_list, file = file.path(save_dir, "04_processed_rna_list_dbl.qs"))
  }
  return(processed_list)
}
