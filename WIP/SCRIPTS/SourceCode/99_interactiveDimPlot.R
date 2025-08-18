library(Seurat)
library(plotly)
library(ggplot2)

# Read previous object
integrated <- qs::qread(file.path(qsave_dir, "10_metadata_edited_obj.qs"))

# Extract UMAP coordinates and metadata
umap_coords <- Embeddings(integrated, reduction = "umap")
metadata <- integrated@meta.data

# Combine coordinates with metadata
plot_data <- data.frame(
  UMAP_1 = umap_coords[,1],
  UMAP_2 = umap_coords[,2],
  cluster = Idents(integrated),
  cell_barcode = rownames(metadata),
  metadata
)

# Create interactive plot with custom hover
p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, 
                           color = cluster,
                           text = paste("Cluster:", cluster,
                                        "<br>TimeGeno:", TimeGeno,
                                        "<br>nFeature_RNA:", nFeature_RNA,
                                        "<br>nCount_RNA:", nCount_RNA,
                                        "<br>mt_ran:", percent_mt_rna,
                                        "<br>CellTypes:", CellTypes,
                                        "<br>ref_short:", ref_short
                                        ))) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "Interactive UMAP")

ggplotly(p, tooltip = "text")
