#### directory ####
qc_dir <- "./QC"
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
  dir.create(file.path(qc_dir, "Plots"), recursive = TRUE)
  }

# read object
qcobject <- qs::qread("./QSAVE/10_metadata_edited_obj.qs")

#### QC fueature plots ####
mt_rna <- FeaturePlot(object = qcobject,
                      features = "percent_mt_rna",
                      reduction = "umap",  # Use UMAP or tSNE for dimensionality reduction
                      label = F, # Label clusters
                      cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                      pt.size = 0.2               # Adjust point size for better visibility
                      )

nCount <- FeaturePlot(object = qcobject,
                      features = "nCount_RNA",
                      reduction = "umap", # Use UMAP or tSNE for dimensionality reduction
                      label = F, # Label clusters
                      cols = c("lightgrey", "blue"),# Adjust colors to highlight high/low values
                      pt.size = 0.2               # Adjust point size for better visibility
                      )

nFeature <- FeaturePlot(object = qcobject,
                        features = "nFeature_RNA",
                        reduction = "umap", # Use UMAP or tSNE for dimensionality reduction
                        label = F, # Label clusters
                        cols = c("lightgrey", "blue"), # Adjust colors to highlight high/low values
                        pt.size = 0.2               # Adjust point size for better visibility
                        )

# combined
RNAqcplot <- mt_rna | nCount | nFeature

# save
ggsave(mt_rna, filename = file.path(qc_dir, "Plots", "mt_rna.png"), width = 8, height = 8, dpi = 300)
ggsave(nCount, filename = file.path(qc_dir, "Plots", "nCount.png"), width = 8, height = 8, dpi = 300)
ggsave(nFeature, filename = file.path(qc_dir, "Plots", "nFeature.png"), width = 8, height = 8, dpi = 300)

#### matrics ####
qc_metrics <- c( "nCount_RNA", "nFeature_RNA", "percent_mt_rna")
qc_thresholds <- list(c(1000, 25000), c(400, NA), c(NA, 5))

# get cell counts for each matrics
qc_pass <- data.frame(
  RNA_count_pass = sum(qcobject$nCount_RNA >= 1000 & qcobject$nCount_RNA <= 25000),
  RNA_feature_pass = sum(qcobject$nFeature_RNA >= 400),
  MT_percent_pass = sum(qcobject$percent_mt_rna <= 5),
  )

# Calculate percentages
qc_pass_percent <- round(qc_pass / ncol(qcobject) * 100, 2)
qc_pass_table <- rbind(qc_pass, qc_pass_percent)

# Calculate percentages
qc_pass_percent <- round(qc_pass / ncol(qcobject) * 100, 2)
qc_pass_table <- rbind(qc_pass, qc_pass_percent)

summary <-
  data.frame(summary(qcobject@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent_mt_rna")]))

write.csv(summary, file = paste0(qc_dir, "QC_matrics.csv"))

#### Distribution ####
p1 <- ggplot(qcobject@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mt_rna)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "RNA metrics",
       x = "UMI count",
       y = "Gene count",
       color = "Percent mitochondrial")

p2 <- ggplot(qcobject@meta.data, aes(x = orig.ident, y = nCount_RNA)) +
  geom_violin(aes(fill = orig.ident), show.legend = FALSE) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  theme_minimal() +
  labs(title = "RNA count distribution",
       x = "Sample",
       y = "UMI count")

p3 <- ggplot(qcobject@meta.data, aes(x = orig.ident, y = percent_mt_rna)) +
  geom_violin(aes(fill = orig.ident), show.legend = FALSE) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  theme_minimal() +
  labs(title = "Mitochondrial percentage",
       x = "Sample",
       y = "Percent mitochondrial")


# Combine plots
combined_plot <- p1 + p2 + p3
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

# distribution plots
p <- ggplot(qcobject@meta.data, aes(x = .data[[metric]])) +
  geom_histogram(bins = 100, fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = paste(metric, "distribution"),
       x = metric,
       y = "Count")

p <- p + geom_vline(xintercept = threshold_low, color = "red", linetype = "dashed")
p <- p + geom_vline(xintercept = threshold_high, color = "red", linetype = "dashed")