library(slingshot)
library(ggplot2)
library(tidyverse)
library(viridis)

# Read previous object if needed.
integrated_man_npc <- qs::qread(file.path(qsave_dir, "OligoSubset", "10__man_npc.qs"))

#Pseudotime prep subset object
OPC_subset <- subset(integrated_man_npc, subset = CellTypes2 %in% c("OPC_1", "OPC_2", "OPC_3"))
umap <- Embeddings(OPC_subset, "umap")
celltype <- OPC_subset$CellTypes2

# Pseudotime Sligshot
lineages <- slingshot(umap, clusterLabels = celltype, reducedDim = "UMAP")

# Convert the UMAP embeddings into a data frame for plotting
umap_df <- as.data.frame(umap)
colnames(umap_df) <- c("umap_1", "umap_2")
umap_df$cell_type <- celltype
umap_df <- cbind(umap_df, OPC_subset@meta.data)

# Extract pseudotime from the first lineage (if multiple lineages exist, you can choose the one of interest)
umap_df$pseudotime <-  lineages@assays@data@listData[["pseudotime"]]

# Plot UMAP with pseudotime as color
pseudo <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point(size = 1) +
  scale_color_viridis_c(na.value = "grey50") +
  labs(title = "Pseudotime Trajectory (Slingshot)",
       color = "Pseudotime") +
  theme_minimal()

condition <-ggplot(umap_df, aes(x = umap_1, y = umap_2, color = pseudotime)) +
  geom_point(size = 1) +
  scale_color_viridis_c(na.value = "grey50") +
  facet_wrap(~ TimeGeno) +  # adjust "condition" as needed
  labs(title = "Pseudotime by Condition",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()

# save
ggsave(filename = file.path(plot_dir, "Trajectory", "pseudotime.png"), plot = pseudo, width = 12, height = 12, dpi = 300)
ggsave(filename = file.path(plot_dir, "Trajectory", "pseudotime2.png"), plot = condition, width = 12, height = 12, dpi = 300)

# density plot
density <- ggplot(umap_df, aes(x = pseudotime, color = TimeGeno, fill = TimeGeno)) +
  geom_density(alpha = 0.3) +
  labs(title = "Pseudotime Density by Time Point", x = "Pseudotime", y = "Density") +
  theme_minimal() +
  scale_color_viridis_d() +
  scale_fill_viridis_d()

# bar plot
bar <- ggplot(umap_df, aes(x = TimeGeno, y = pseudotime, fill = TimeGeno)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Pseudotime Distribution by Time Point", x = "Time", y = "Pseudotime") +
  theme_minimal() +
  scale_fill_viridis_d()

# save
ggsave(filename =  file.path(plot_dir, "Trajectory", "density.png"), plot = density, width = 12, height = 12, dpi = 300)
ggsave(filename = file.path(plot_dir, "Trajectory", "bar.png"), plot = bar, width = 12, height = 12, dpi = 300)
