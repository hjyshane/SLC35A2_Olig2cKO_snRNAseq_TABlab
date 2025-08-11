# read this function with
# source('~/SLC35A2_Olig2cKO_snRNA/WIP/SCRIPTS/SourceCode/02_qc_Initial_VlnPlot.R')

#' Function to create VlnPlots.
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
# Normalize for visualization
qc_VlnPlot <-
    function(
        seurat_list,
        feature_list = c("nFeature_RNA", "nCount_RNA", "percent_mt_rna"),
        ylim = c(20000, 200000, 20),
        prefix = "Vln_initial",
        save = TRUE,
        save_dir = NULL) {
        # Save check
        if (save) {if (is.null(save_dir)) {stop("You must provide 'save_dir' when save = TRUE.")}
            if (!dir.exists(save_dir)) {dir.create(save_dir, recursive = TRUE)}
        }

        # Initiate plot list
        plots <- list()

        # Normalize for visualization
        for (sample in names(seurat_list)) {
            seurat_list[[sample]] <- NormalizeData(seurat_list[[sample]])
        }

        # create VlnPlot for nCount, nFeature, and mt_rna
        for (i in seq_along(feature_list)) {
            for (sample in names(seurat_list)) {
                plots[[feature_list[[i]]]][[sample]] <-
                    VlnPlot(seurat_list[[sample]],
                            features = feature_list[i],
                            alpha = 0) +
                    ylim(0, ylim[i]) +
                    NoLegend() +
                    theme(plot.title = element_blank(),
                          axis.title.x = element_blank(),
                          axis.text.x = element_text(angle = 0,
                                                     hjust = 0.5))
            }
            # Combine them
            feature_combined <-
                patchwork::wrap_plots(plots[[i]],
                                      col = 4) +
                patchwork::plot_annotation(feature_list[i],
                                           theme = theme(plot.title = element_text(size = 24,
                                                                                   hjust = 0.5)))
            # Save
            ggsave(feature_combined,
                   file = file.path(save_dir, paste0(prefix, "_", feature_list[i], ".png")),
                   width = 12,
                   height = 12,
                   dpi = 300,
                   bg = "white")
        }
        return(plots)
    }
