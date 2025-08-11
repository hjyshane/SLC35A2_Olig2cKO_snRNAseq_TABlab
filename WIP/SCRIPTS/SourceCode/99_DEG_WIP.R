# Cell type and Timpoints
clusters_to_keep <- cell_counts %>%
    filter(P00_cKO >= 10 & P00_Control >= 10 |
               P08_cKO >= 10 & P08_Control >= 10 |
               P23_cKO >= 10 & P23_Control >= 10 |
               P56_cKO >= 10 & P56_Control >= 10) %>%
    pull(Cell_Type)

timepoints <- unique(oligsub_man2$TimePoint)

# Set DEG save directory
deg_dir <- file.path(save_dir, "DEG")
bg_dir <- file.path(save_dir, "DEG", "bg_csv")
sig_dir <- file.path(save_dir, "DEG", "sig_csv")
plot_dir <- file.path(save_dir, "DEG", "Plots", "VolcanoPlot")

dir.create(deg_dir, recursive = TRUE)
dir.create(bg_dir, recursive = TRUE)
dir.create(sig_dir, recursive = TRUE)
dir.create(plot_dir, recursive = TRUE)

Seurat::DefaultAssay(oligsub_man2) <- "RNA"
oligsub_man2 <-Seurat:: NormalizeData(oligsub_man2)
oligsub_man2 <- Seurat::ScaleData(oligsub_man2)
oligsub_man2 <- Seurat::SCTransform(oligsub_man2)
oligsub_man2 <- Seurat::PrepSCTFindMarkers(oligsub_man2, assay = "SCT", verbose = T)

# Cell count data frame
cell_counts <-
    table(oligsub_man2[["ref_short"]][,1], oligsub_man2[["TimePoint"]][,1]) %>%
    as.data.frame.matrix() %>%
    tibble::rownames_to_column(var = "Cell_Type")

# Loop through cell type
for (time in timepoints){
    print(time)

    time_data <- subset(oligsub_man2,
                        subset = TimePoint %in% time)

    # For each cell type, loop through each comparison group
    for (ct in clusters_to_keep)  {
        # print(ct)
        # print(time_data[[]] %>% filter(ref_short == ct) %>% group_by(GenoType) %>% summarize(n = n()) %>% dplyr::select(n))
        # print(time_data[[]] %>% filter(ref_short == ct) %>% group_by(GenoType) %>% summarize(n = n()) %>% dplyr::select(n) < 10)

        if (!(ct %in% time_data[[]]$ref_short) | any(time_data[[]] %>% filter(ref_short == ct) %>% group_by(GenoType) %>% summarize(n = n()) %>% dplyr::select(n) < 10)) {
            message(glue::glue("{ct} @ {time}: fewer than 10 cells, skipped."))
            next
        }

        dg_data <- subset(time_data,
                          subset = ref_short == ct)
        dg_data[[]] %>% group_by(GenoType) %>% summarize(n = n()) %>% print

        # if (ncol(dg_data) < 10) {
        #
        #   next
        # }

        # subsetting cell type in each condition
        dg_data <- Seurat::NormalizeData(dg_data)


        if (length(unique(dg_data$GenoType)) == 2) {
            # time_data <- PrepSCTFindMarkers(time_data)

            # Set identities
            Seurat::Idents(dg_data) <- "GenoType"
            Seurat::DefaultAssay(dg_data) <- "RNA"
            # DE analysis
            bg_results <- Seurat::FindMarkers(
                object = dg_data,
                assay = "RNA",
                group.by = "GenoType",
                ident.1 = "cKO",
                ident.2 = "Control",
                logfc.threshold = 0.0,
                min.pct = 0.25,
                test.use = "wilcox"
            )

            # Add metadata info
            bg_results$gene <- rownames(bg_results)
            bg_results$cell_type <- ct
            bg_results$TimePoint <- time

            # Add regulation info
            bg_results <- bg_results %>%
                mutate(Regulation = case_when(
                    avg_log2FC > 0 ~ "Up",
                    avg_log2FC < 0 ~ "Down",
                    TRUE ~ "NoChange"))

            # Get significant genes - DEGs
            sig_results <- bg_results %>%
                filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.2)

            # Create VolcanoPlot
            if (nrow(bg_results) == 0) {
                message("There are no genes/peaks to plot")
                next
            } else {
                vplot <- EnhancedVolcano::EnhancedVolcano(
                    bg_results,
                    lab = bg_results$gene,
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = paste0('VolcanoPlot for ', ct),
                    subtitle = time,
                    pCutoff = 0.05,
                    FCcutoff = 0.2,
                    pointSize = 1.0,
                    labSize = 3.0)
            }

            # Save
            if (TRUE) {
                if (nrow(bg_results) == 0) {
                    message("There are no genes/peaks to found")
                    next
                } else {
                    write.csv(bg_results, file = file.path(bg_dir, paste0(ct, "_", time, ".csv")))
                    ggplot2::ggsave(plot = vplot, file = file.path(plot_dir, paste0(ct, "_", time, ".png")), width = 10, height = 12, dpi = 300)
                }
                if (nrow(sig_results) == 0) {
                    message("There are no deg genes/peaks to found")
                    next
                } else {
                    write.csv(sig_results, file = file.path(sig_dir, paste0(ct, "_", time, ".csv")))
                }

            }


        } else {next}
    }
}

# Run DEG
run_dg(
    oligsub_man2,
    cell_types,
    timepoints,
    minimum_cnt = 3, # 10 or higher for MAST
    mode = "wilcox", # or "MAST"
    p_value = 0.05,
    fc_value = 0.2,
    save = TRUE,
    bg_dir = bg_dir,
    sig_dir = sig_dir,
    plot_dir = plot_dir)