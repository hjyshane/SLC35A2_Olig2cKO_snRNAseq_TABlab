# read this function with
# source("~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/01_load_rna_from_h5.R")
#' Load multiple RNA expression matrices from 10X .h5 files
#'
#' This function loads RNA count matrices (e.g., "Gene Expression") from multiple
#' 10X Genomics HDF5 files. Optionally, it can save each matrix using `qs::qsave()`.
#'
#' @param input_dir Character. Directory containing sample folders.
#' @param samples Named list of character vectors. Names are index identifiers,
#'        and each vector contains sample names for that index. Each sample must have
#'        the file structure: `/input_dir/index/outs/per_sample_outs/sample/count/sample_filtered_feature_bc_matrix.h5`.
#' @param assay Character. Assay name to extract from each H5 file (default: "Gene Expression").
#' @param save Logical. If TRUE, saves each matrix using `qs::qsave()` (default: FALSE).
#' @param save_dir Character or NULL. Output directory for saved `.qs` files (one per sample). Required if `save = TRUE`.
#'
#' @return Named list of sparse matrices (class `dgCMatrix`), one per sample.
#'
#' @examples
#' # Define sample structure
#' samples <- list("FF061" = c("SOM333", "SOM318", "SOM303", "SOM301"),
#'                 "FF062" = c("SOM326", "SOM325", "SOM323", "SOM322"))
#'
#' # Load RNA data
#' rna_list <- load_rna_from_h5(input_dir = "~/project/data",
#'                              samples = samples,
#'                              save_dir = "~/project/saved_data")
#'
#' @export

load_rna_from_h5 <-
  function(
    input_dir,
    samples,
    save = FALSE,
    save_dir = NULL){
    # Save check
    if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
      if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
      }

    # Create empty list to store returned object
    rna_list <- list()
    total_samples <- sum(lengths(samples))

    # Read each sample's h5 files and store into rna_list[[sample]]
    message("Loading ", total_samples, " samples.")
    for (index in names(samples)) {
      for (sample in samples[[index]]) {

        # Set file path
        h5_path <- file.path(input_dir, index, 'outs', 'per_sample_outs', sample, 'count', 'sample_filtered_feature_bc_matrix.h5')

        if (file.exists(h5_path)) {
          message("Loading sample: ", index, "-", sample, " from ", h5_path)

          # Read h5 and Store in rna_list
          h5 <- Seurat::Read10X_h5(h5_path)
          rna_list[[sample]] <- h5
          # Save if save is enabled
          if (save){
            save_name <- file.path(save_dir, "INDIVIDUAL_h5", paste0(sample, "_h5.qs"))
            dir.create(file.path(save_dir, "INDIVIDUAL_h5"), recursive = T)

            qs::qsave(h5, file = save_name)
            message("Saved object: ", sample, " to ", save_name, ".")
            }
          } else {
            message("File not found for sample: ", sample, "\nMissing in ", h5_path)
            next
          }
        }
      message("Finished loading ", length(rna_list), " of ", total_samples, " samples.")
    }
    return(rna_list)
  }



