library(tidyverse)
library(ggplot2)
library(speckle)

set.seed(0827)

# Read object
integrated <- qs::qread(file.path(qsave_dir, "09_ref_annotated_obj.qs"))

# Directory
save_dir <-
    '~/SLC35A2_Olig2cKO_snRNA/WIP'
qsave_dir <-
    file.path(save_dir, "QSAVE")
csv_dir <-
    file.path(save_dir, 'CSV')
plot_dir <-
    file.path(save_dir, 'PLOTS')

#### Cell count Visualization based on ref_short####
# Metadata
meta <- integrated@meta.data %>%
    mutate(
        lineage_group = case_when(
            ref_short %in% c("Ex_CA1", "Ex_CA3", "Ex_Cortex", "Granule_DG", "Granule_Blast_DG") ~ "Excitatory Neuron",
            ref_short %in% c("Interneuron", "CGE", "MGE") ~ "Inhibitory Neuron",
            ref_short %in% c("Astrocytes", "Microglia", "NPC", "RGL_DG") ~ "Glial",
            ref_short %in% c("OPC", "COP", "MFOL", "NFOL", "Mature_OL") ~ "Oligodendrocytes",
            ref_short %in% c("Macrophages") ~ "Immune",
            ref_short %in% c("Vascular", "Pericytes") ~ "Vascular",
            TRUE ~ "Other"
        )
    )

lineage_colors <- c(
    "Excitatory Neuron" = "#1f78b4",
    "Inhibitory Neuron" = "#e31a1c",
    "Glial" = "#33a02c",
    "Oligodendrocytes" = "orange",
    "Immune" = "yellow",
    "Vascular" = "#6a3d9a",
    "Other" = "gray70"
)

# Cell counts per type
celltype_GenoType_counts <- meta %>%
    count(Cell_type = ref_short, GenoType, TimePoint, MouseID, TimeGeno, name = "n_cells") %>%
    arrange(desc(n_cells))

all_cell <- ggplot(celltype_GenoType_counts, aes(x = MouseID, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
        title = "Stacked Cell Type Counts per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_cell_count <- ggplot(celltype_GenoType_counts, aes(x = MouseID, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Cell Type Counts per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_cell_tg <- ggplot(celltype_GenoType_counts, aes(x = TimeGeno, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
        title = "Stacked Cell Type Counts GenoXTime",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_cell_tg_count <- ggplot(celltype_GenoType_counts, aes(x = TimeGeno, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Cell Type Counts GenoXTime",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

mouse_lineage_counts <- meta %>%
    count(GenoType, TimePoint, MouseID, TimeGeno, lineage_group)

cell_gruop <- ggplot(mouse_lineage_counts, aes(x = MouseID, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Stacked Cell Type Counts Mouse (Lineage grouped)",
        y = "Proportion of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

cell_gruop_count <- ggplot(mouse_lineage_counts, aes(x = MouseID, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Cell Type Counts Mouse (Lineage grouped)",
        y = "Proportion of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

cell_gruop_geno <- ggplot(mouse_lineage_counts, aes(x = TimeGeno, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Stacked Cell Type Counts TimeXGeno (Lineage grouped)",
        y = "Proportion of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

cell_gruop_geno_numb <- ggplot(mouse_lineage_counts, aes(x = TimeGeno, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Cell Type Counts TimeXGeno (Lineage grouped)",
        y = "Number of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

Olig <- c("OPC", "COP", "NFOL", "MFOL", "Mature_OL")

mouse_celltype_counts <- celltype_GenoType_counts %>%
    filter(Cell_type %in% Olig)

fill_colors <- c(
    "OPC" = "#1b9e77",
    "COP" = "#d95f02",
    "NFOL" = "#7570b3",
    "MFOL" = "#e7298a",
    "Mature_OL" = "#66a61e")

geno_cell_mouse <- ggplot(mouse_celltype_counts, aes(x = MouseID, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    # scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_ratio <- ggplot(mouse_celltype_counts, aes(x = MouseID, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_lineage <- ggplot(mouse_celltype_counts, aes(x = TimeGeno, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Oligodendrocytes Composition TimeXGeno",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell <- ggplot(mouse_celltype_counts, aes(x = TimeGeno, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Cell Type Composition TimeXGeno",
        x = "Mouse ID",
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

Olig_NPC <- c("OPC", "COP", "NFOL", "MFOL", "Mature_OL", "NPC")

mouse_celltype_counts_NPC <- celltype_GenoType_counts %>%
    filter(Cell_type %in% Olig_NPC)

fill_colors_NPC <- c(
    "OPC" = "#1b9e77",
    "COP" = "#d95f02",
    "NFOL" = "#7570b3",
    "MFOL" = "#e7298a",
    "Mature_OL" = "#66a61e",
    "NPC" = "yellow")

geno_cell_mouse_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = MouseID, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    # scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_ratio_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = MouseID, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_lineage_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = TimeGeno, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Oligodendrocytes Composition TimeXGeno",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = TimeGeno, y = n_cells, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Cell Type Composition TimeXGeno",
        x = "Mouse ID",
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Cell count Visualization based on CellTypes####
# Metadata
meta <- integrated@meta.data %>%
    mutate(
        lineage_group_2 = case_when(
            CellTypes %in% c("ExNeurons") ~ "Excitatory Neuron",
            CellTypes %in% c("InNeurons") ~ "Inhibitory Neuron",
            CellTypes %in% c("Astrocytes", "Microglia", "Ependymal") ~ "Glial",
            CellTypes %in% c("OPC", "MFOL", "NFOL") ~ "Oligodendrocytes",
            CellTypes %in% c("NPC") ~ "NPC",
            CellTypes %in% c("Pericytes") ~ "Pericytes",
            TRUE ~ "Other"
        )
    )

lineage_colors <- c(
    "Excitatory Neuron" = "#1f78b4",
    "Inhibitory Neuron" = "#e31a1c",
    "Glial" = "#33a02c",
    "Oligodendrocytes" = "orange",
    "NPC" = "yellow",
    "Pericytes" = "#6a3d9a",
    "Other" = "gray70"
)

# Cell counts per type
celltype_GenoType_counts <- meta %>%
    dplyr::count(CellTypes, GenoType, TimePoint, MouseID, TimeGeno, name = "n_cells") %>%
    arrange(desc(n_cells))

all_cell <- ggplot(celltype_GenoType_counts, aes(x = MouseID, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
        title = "Stacked Cell Type Counts per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_cell_count <- ggplot(celltype_GenoType_counts, aes(x = MouseID, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Cell Type Counts per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_cell_tg <- ggplot(celltype_GenoType_counts, aes(x = TimeGeno, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
        title = "Stacked Cell Type Counts GenoXTime",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_cell_tg_count <- ggplot(celltype_GenoType_counts, aes(x = TimeGeno, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Cell Type Counts GenoXTime",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

mouse_lineage_counts <- meta %>%
    dplyr::count(GenoType, TimePoint, MouseID, TimeGeno, lineage_group)

cell_gruop <- ggplot(mouse_lineage_counts, aes(x = MouseID, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Stacked Cell Type Counts Mouse (Lineage grouped)",
        y = "Proportion of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

cell_gruop_count <- ggplot(mouse_lineage_counts, aes(x = MouseID, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Cell Type Counts Mouse (Lineage grouped)",
        y = "Proportion of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

cell_gruop_geno <- ggplot(mouse_lineage_counts, aes(x = TimeGeno, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Stacked Cell Type Counts TimeXGeno (Lineage grouped)",
        y = "Proportion of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

cell_gruop_geno_numb <- ggplot(mouse_lineage_counts, aes(x = TimeGeno, y = n, fill = lineage_group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = lineage_colors) +
    labs(
        title = "Cell Type Counts TimeXGeno (Lineage grouped)",
        y = "Number of Cells",
        x = NULL,
        fill = "Lineage Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))

Olig <- c("OPC", "NFOL", "MFOL")

mouse_celltype_counts <- celltype_GenoType_counts %>%
    filter(CellTypes %in% Olig)

fill_colors <- c(
    "OPC" = "#1b9e77",
    "COP" = "#d95f02",
    "NFOL" = "#7570b3",
    "MFOL" = "#e7298a",
    "Mature_OL" = "#66a61e")

geno_cell_mouse <- ggplot(mouse_celltype_counts, aes(x = MouseID, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    # scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_ratio <- ggplot(mouse_celltype_counts, aes(x = MouseID, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_lineage <- ggplot(mouse_celltype_counts, aes(x = TimeGeno, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Oligodendrocytes Composition TimeXGeno",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell <- ggplot(mouse_celltype_counts, aes(x = TimeGeno, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors) +
    labs(
        title = "Cell Type Composition TimeXGeno",
        x = "Mouse ID",
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

Olig_NPC <- c("OPC", "NFOL", "MFOL", "NPC")

mouse_celltype_counts_NPC <- celltype_GenoType_counts %>%
    filter(CellTypes %in% Olig_NPC)

fill_colors_NPC <- c(
    "OPC" = "#1b9e77",
    "NFOL" = "#7570b3",
    "MFOL" = "#e7298a",
    "NPC" = "yellow")

geno_cell_mouse_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = MouseID, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    # scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_ratio_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = MouseID, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Oligodendrocytes Composition per Mouse",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_lineage_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = TimeGeno, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Oligodendrocytes Composition TimeXGeno",
        x = NULL,
        y = "Proportion of Cells",
        fill = "Cell Type"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

geno_cell_NPC <- ggplot(mouse_celltype_counts_NPC, aes(x = TimeGeno, y = n_cells, fill = CellTypes)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_colors_NPC) +
    labs(
        title = "Cell Type Composition TimeXGeno",
        x = "Mouse ID",
        y = "Number of Cells",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Cell count propeller stats ref_short ####
# Extract the vectors we need from your Seurat object
# for Control only
cell_clusters <- integrated[, integrated$GenoType == "Control"]$ref_short        # Individual test batch
cell_samples <- integrated[, integrated$GenoType == "Control"]$MouseID           # Sample ID for each cell
cell_groups <- integrated[, integrated$GenoType == "Control"]$TimePoint           # Comparison group

# cKO only
# for Control only
cell_clusters <- integrated[, integrated$GenoType == "cKO"]$ref_short       # Individual test batch
cell_samples <- integrated[, integrated$GenoType == "cKO"]$MouseID           # Sample ID for each cell
cell_groups <- integrated[, integrated$GenoType == "cKO"]$TimePoint           # Comparison group

# control and cKO
cell_clusters <- integrated$ref_short        # Individual test batch
cell_samples <- integrated$MouseID           # Sample ID for each cell
cell_groups <- integrated$TimePoint            # Comparison group


# Overall analysis across all TimePoints
# From data(cluster and sample IDs), build a table of counts per sample per cluster type.
# Compute proportions (count / total cells per sample).
propeller_overall <- propeller(
    clusters = cell_clusters,
    sample = cell_samples,
    group = cell_groups
)

# Normalize proportion values before modeling.
# For two groups: moderated t-test
# For >2 groups: moderated ANOVA
# Adjust p-values for multiple testing using Benjamini–Hochberg FDR

props_for_plot <- getTransformedProps(
    clusters = cell_clusters,
    sample = cell_samples,
    transform = "logit" # asin is more robust when proportions include 0 or 1; logit has more power when proportions are moderate
)

# asin test
props_for_plot <- getTransformedProps(
    clusters = cell_clusters,
    sample = cell_samples,
    transform = "asin" # asin is more robust when proportions include 0 or 1; logit has more power when proportions are moderate
)

sample_order <- colnames(props_for_plot$Proportions)
sample_info <- integrated@meta.data %>%
    select(MouseID, TimePoint, GenoType, TimeGeno) %>%
    distinct() %>%
    arrange(MouseID)
plot_groups <- sample_info$TimePoint[match(sample_order, sample_info$MouseID)]

# Get significant cell types (FDR < 0.05)
sig_celltypes <- rownames(propeller_overall)[propeller_overall$FDR < 0.05]

# If no significant results, show top 5 by p-value for demonstration
if(length(sig_celltypes) == 0) {
    top_celltypes <- rownames(propeller_overall)[order(propeller_overall$P.Value)][1:5]
    print("Top 5 cell types by p-value:")
    print(propeller_overall[top_celltypes, ])
}

# Define colors for your 8 developmental groups
group_colors <- c("P00_cKO" = "red", "P00_Control" = "lightcoral",
                  "P08_cKO" = "blue", "P08_Control" = "lightblue",
                  "P23_cKO" = "darkgreen", "P23_Control" = "lightgreen",
                  "P56_cKO" = "purple", "P56_Control" = "plum")

group_colors <- c("P00" = "lightcoral", "P08" = "lightblue", "P23" = "lightgreen", "P56" = "plum")

png(file.path(plot_dir, "Cell_Counts", "ProportionTest", "Propeller_oligodendrocyte_lineage_asin.png"), width = 1200, height = 800, res = 100)
par(mfrow = c(2, 3))
par(mar = c(8, 5, 3, 2))
for(celltype in ol_types) {
   stripchart(as.numeric(props_for_plot$Proportions[celltype, ]) ~ plot_groups,
               vertical = TRUE,
               pch = 16,
               method = "jitter",
               col = group_colors[col_groups],
               cex = 2,
               ylab = "Proportions",
               cex.axis = 1.2,
               cex.lab = 1.4,
               las = 2)

    title(paste(celltype, "Development"), cex.main = 1.4, adj = 0)

    # Add FDR annotation
    fdr_val <- propeller_overall[celltype, "FDR"]
    text(4, max(props_for_plot$Proportions[celltype, ]) * 0.85,
         labels = paste("FDR =", round(fdr_val, 4)),
         cex = 1.2, col = "black", bg = "white")

}

dev.off()


png(file.path(plot_dir, "Cell_Counts", "ProportionTest", "Propeller_all_Control_asin.png"), width = 1200, height = 800, res = 100)
par(mfrow = c(7, 3))
par(mar=c(1,1,1,1))
for(celltype in sig_celltypes) {
    stripchart(as.numeric(props_for_plot$Proportions[celltype, ]) ~ plot_groups,
               vertical = TRUE,
               pch = 16,
               method = "jitter",
               col = group_colors[plot_groups],
               cex = 2,
               ylab = "Proportions",
               cex.axis = 1.2,
               cex.lab = 1.4,
               las = 2)

    title(celltype, cex.main = 1.4, adj = 0)

    # Add FDR annotation
    fdr_val <- propeller_overall[celltype, "FDR"]
    text(4, max(props_for_plot$Proportions[celltype, ]) * 0.85,
         labels = paste("FDR =", formatC(fdr_val, format = "e", digits = 2)),
         cex = 1.2)
}
dev.off()

#### Cell count propeller stats CellTypes ####
# Extract the vectors we need from your Seurat object
# for Control only
cell_clusters <- integrated[, integrated$GenoType == "Control"]$CellTypes        # Individual test batch
cell_samples <- integrated[, integrated$GenoType == "Control"]$MouseID           # Sample ID for each cell
cell_groups <- integrated[, integrated$GenoType == "Control"]$TimePoint           # Comparison group

# cKO only
# for Control only
cell_clusters <- integrated[, integrated$GenoType == "cKO"]$CellTypes       # Individual test batch
cell_samples <- integrated[, integrated$GenoType == "cKO"]$MouseID           # Sample ID for each cell
cell_groups <- integrated[, integrated$GenoType == "cKO"]$TimePoint           # Comparison group

# control and cKO
cell_clusters <- integrated$CellTypes        # Individual test batch
cell_samples <- integrated$MouseID           # Sample ID for each cell
cell_groups <- integrated$TimePoint            # Comparison group


# Overall analysis across all TimePoints
# From data(cluster and sample IDs), build a table of counts per sample per cluster type.
# Compute proportions (count / total cells per sample).
propeller_overall <- propeller(
    clusters = cell_clusters,
    sample = cell_samples,
    group = cell_groups
)

# Normalize proportion values before modeling.
# For two groups: moderated t-test
# For >2 groups: moderated ANOVA
# Adjust p-values for multiple testing using Benjamini–Hochberg FDR

props_for_plot <- getTransformedProps(
    clusters = cell_clusters,
    sample = cell_samples,
    transform = "logit" # asin is more robust when proportions include 0 or 1; logit has more power when proportions are moderate
)

# asin test
props_for_plot <- getTransformedProps(
    clusters = cell_clusters,
    sample = cell_samples,
    transform = "asin" # asin is more robust when proportions include 0 or 1; logit has more power when proportions are moderate
)

sample_order <- colnames(props_for_plot$Proportions)
sample_info <- integrated@meta.data %>%
    select(MouseID, TimePoint, GenoType, TimeGeno) %>%
    distinct() %>%
    arrange(MouseID)
plot_groups <- sample_info$TimePoint[match(sample_order, sample_info$MouseID)]

# Get significant cell types (FDR < 0.05)
sig_celltypes <- rownames(propeller_overall)[propeller_overall$FDR < 0.05]

# If no significant results, show top 5 by p-value for demonstration
if(length(sig_celltypes) == 0) {
    top_celltypes <- rownames(propeller_overall)[order(propeller_overall$P.Value)][1:5]
    print("Top 5 cell types by p-value:")
    print(propeller_overall[top_celltypes, ])
}

# Define colors for your 8 developmental groups
group_colors <- c("P00_cKO" = "red", "P00_Control" = "lightcoral",
                  "P08_cKO" = "blue", "P08_Control" = "lightblue",
                  "P23_cKO" = "darkgreen", "P23_Control" = "lightgreen",
                  "P56_cKO" = "purple", "P56_Control" = "plum")

group_colors <- c("P00" = "lightcoral", "P08" = "lightblue", "P23" = "lightgreen", "P56" = "plum")


png(file.path(plot_dir, "Cell_Counts", "ProportionTest", "Propeller_oligodendrocyte_lineage_asin_cKO.png"), width = 1200, height = 800, res = 100)
par(mfrow = c(2, 3))
par(mar = c(8, 5, 3, 2))
for(celltype in Olig) {
    stripchart(as.numeric(props_for_plot$Proportions[celltype, ]) ~ plot_groups,
               vertical = TRUE,
               pch = 16,
               method = "jitter",
               col = group_colors[plot_groups],
               cex = 2,
               ylab = "Proportions",
               cex.axis = 1.2,
               cex.lab = 1.4,
               las = 2)

    title(paste(celltype, "Development"), cex.main = 1.4, adj = 0)

    # Add FDR annotation
    fdr_val <- propeller_overall[celltype, "FDR"]
    text(4, max(props_for_plot$Proportions[celltype, ]) * 0.85,
         labels = paste("FDR =", round(fdr_val, 4)),
         cex = 1.2, col = "black", bg = "white")

}

dev.off()


png(file.path(plot_dir, "Cell_Counts", "ProportionTest", "Propeller_all_cKO_asin.png"), width = 1200, height = 800, res = 100)
par(mfrow = c(2, 3))
par(mar=c(1,1,1,1))
for(celltype in sig_celltypes) {
    stripchart(as.numeric(props_for_plot$Proportions[celltype, ]) ~ plot_groups,
               vertical = TRUE,
               pch = 16,
               method = "jitter",
               col = group_colors[plot_groups],
               cex = 2,
               ylab = "Proportions",
               cex.axis = 1.2,
               cex.lab = 1.4,
               las = 2)

    title(celltype, cex.main = 1.4, adj = 0)

    # Add FDR annotation
    fdr_val <- propeller_overall[celltype, "FDR"]
    text(4, max(props_for_plot$Proportions[celltype, ]) * 0.85,
         labels = paste("FDR =", formatC(fdr_val, format = "e", digits = 2)),
         cex = 1.2)
}
dev.off()

#### Cell count scProportionTest ####
prop_test <- sc_utils(integrated)
## todo ## # sanity check among time point
# n.s
prop_test <- permutation_test(
    prop_test, cluster_identity = "ref_short",
    sample_1 = "cKO", sample_2 = "Control",
    sample_identity = "TimeGeno"
)

ckocontrol <-  permutation_plot(prop_test)


prop_test_P00 <- permutation_test(
    prop_test, cluster_identity = "ref_short",
    sample_1 = "P00_cKO", sample_2 = "P00_Control",
    sample_identity = "TimeGeno"
)

p00 <-  permutation_plot(prop_test)

prop_test_P08 <- permutation_test(
    prop_test, cluster_identity = "ref_short",
    sample_1 = "P08_cKO", sample_2 = "P08_Control",
    sample_identity = "TimeGeno"
)

p08 <-  permutation_plot(prop_test)

prop_test_P23 <- permutation_test(
    prop_test, cluster_identity = "ref_short",
    sample_1 = "P23_cKO", sample_2 = "P23_Control",
    sample_identity = "TimeGeno"
)

p23 <-  permutation_plot(prop_test)

prop_test_P56 <- permutation_test(
    prop_test, cluster_identity = "ref_short",
    sample_1 = "P56_cKO", sample_2 = "P56_Control",
    sample_identity = "TimeGeno"
)

p56 <-  permutation_plot(prop_test)

p00 / p08 | p23 / p56

#### Normalized counts ####
# Find which row Slc35a2 is in
gene_row <- which(rownames(integrated[["RNA"]]) == "Slc35a2")

# Get counts data of Slc35a2
Slc35a2_data <- integrated[["RNA"]]@layers$data[gene_row, ]

# see the distribution
length(Slc35a2_data) # Total cells : 86899
hist(Slc35a2_data)

# get logical vector for Slc35a2 expressing cells
Slc35a2_expressed_data <- Slc35a2_data > 0

# Basic stats
sum(Slc35a2_expressed_data) # Slc35a2 expressing cell : 34837
mean(Slc35a2_expressed_data) # 0.4008907 -> ~40% cell express Slc35a2
mean(Slc35a2_data) # Average expression level of Slc35a2 : 0.3903355
mean(Slc35a2_data[Slc35a2_data > 0]) # Slc35a2 expressing cell average expression level : 0.9736706

# Add meta data
meta$Slc35a2_expressed_data <- Slc35a2_expressed_data
meta$Slc35a2_data <- as.numeric(Slc35a2_data)

# Cell tyep X GenoType
# ref_short
expression_summary_data <- meta %>%
    group_by(ref_short, GenoType, TimePoint, TimeGeno, lineage_group) %>%
    summarise(
        n_cells = n(),
        n_expressed = sum(Slc35a2_expressed_data),
        percent_expressed = round(100 * mean(Slc35a2_expressed_data), 2),
        expression_level = round(mean(Slc35a2_data[Slc35a2_expressed_data]), 2),
        expression_all_cells = round(mean(Slc35a2_data), 2)
    ) %>%
    arrange(desc(percent_expressed))

expression_summary_data_subset <- expression_summary_data %>%
    filter(ref_short %in% c("OPC", "COP", "NFOL", "MFOL", "Mature_OL", "NPC"))

# CellTypes
expression_summary_data_2 <- meta %>%
    group_by(CellTypes, GenoType, TimePoint, TimeGeno, lineage_group) %>%
    summarise(
        n_cells = n(),
        n_expressed = sum(Slc35a2_expressed_data),
        percent_expressed = round(100 * mean(Slc35a2_expressed_data), 2),
        expression_level = round(mean(Slc35a2_data[Slc35a2_expressed_data]), 2),
        expression_all_cells = round(mean(Slc35a2_data), 2)
    ) %>%
    arrange(desc(percent_expressed))

expression_summary_data_subset_2 <- expression_summary_data %>%
    filter(ref_short %in% c("OPC", "NFOL", "MFOL", "NPC"))


# simple viz
slcexpcell <- ggplot(expression_summary_data_subset, aes(x = ref_short, y = percent_expressed, fill = GenoType)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(~TimePoint)+
    labs(
        title = "Slc35a2 Expression Rate per Cell Type × GenoType",
        x = "Cell Type",
        y = "Percent of Expressing Cells"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

slcexc <- ggplot(expression_summary_data_subset, aes(x = ref_short, y = expression_level, fill = GenoType)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(~TimePoint) +
    labs(
        title = "Mean Slc35a2 Expression in Expressing Cells",
        x = "Cell Type",
        y = "Normalized counts"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

slcexpcell <- ggplot(expression_summary_data_subset_2, aes(x = ref_short, y = percent_expressed, fill = GenoType)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(~TimePoint)+
    labs(
        title = "Slc35a2 Expression Rate per Cell Type × GenoType",
        x = "Cell Type",
        y = "Percent of Expressing Cells"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

slcexc <- ggplot(expression_summary_data_subset_2, aes(x = ref_short, y = expression_level, fill = GenoType)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(~TimePoint) +
    labs(
        title = "Mean Slc35a2 Expression in Expressing Cells",
        x = "Cell Type",
        y = "Normalized counts"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Raw counts ####
# Find which row Slc35a2 is in
gene_row <- which(rownames(integrated[["RNA"]]) == "Slc35a2")

# Get counts data of Slc35a2
Slc35a2_counts <- integrated[["RNA"]]@layers$counts[gene_row, ]

# see the distribution
length(Slc35a2_counts) # Total cells : 86899
table(Slc35a2_counts) # summary of each counts
hist(Slc35a2_counts)

# get logical vector for Slc35a2 expressing cells
Slc35a2_expressed <- Slc35a2_counts > 0

# Basic stats
sum(Slc35a2_expressed) # Slc35a2 expressing cell : 34837
mean(Slc35a2_expressed) # 0.4008907 -> ~40% cell express Slc35a2
mean(Slc35a2_counts) # Average expression level of Slc35a2 : 0.5972911
mean(Slc35a2_counts[Slc35a2_expressed]) # Slc35a2 expressing cell average expression level : 1.48991

# Add meta data
meta$Slc35a2_expressed <- Slc35a2_expressed
meta$Slc35a2_counts <- as.numeric(Slc35a2_counts)

# Cell tyep X GenoType
expression_summary_counts <- meta %>%
    group_by(ref_short, GenoType, TimePoint, TimeGeno, lineage_group) %>%
    summarise(
        n_cells = n(),
        n_expressed = sum(Slc35a2_expressed),
        percent_expressed = round(100 * mean(Slc35a2_expressed),2),
        expression_level = round(mean(Slc35a2_counts[Slc35a2_expressed]),2),
        expression_all_cells = round(mean(Slc35a2_counts), 2)
    ) %>%
    arrange(desc(percent_expressed))

# simple viz
ggplot(expression_summary_counts, aes(x = ref_short, y = percent_expressed, fill = GenoType)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Slc35a2 Expression Rate per Cell Type × GenoType",
        x = "Cell Type",
        y = "Percent of Expressing Cells"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(expression_summary_counts, aes(x = ref_short, y = expression_level, fill = GenoType)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
        title = "Mean Slc35a2 Expression in Expressing Cells",
        x = "Cell Type",
        y = "Mean UMI Count (Expressing Cells)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### save ####
integrated@meta.data <- meta

# Save object
qs::qsave(integrated, file = file.path(qsave_dir, "10_metadata_edited_obj.qs"))

# Save plots
ggsave(all_cell, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "Ratio", "PerMouse", "CellCountsPerMouseAllTypes.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(all_cell_count, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts", "PerMouse", "CellCountsPerMouseAllTypesCounts.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(all_cell_tg, filename = file.path(plot_dir, "Cell_Counts","BarPlots", "Ratio", "CellCountsTimeGenoAllTypes.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(all_cell_tg_count, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts", "CellCountsTimeGenoAllTypesCounts.png"), width = 7, height = 8, dpi = 300, bg = "white")

ggsave(cell_gruop, filename = file.path(plot_dir, "Cell_Counts","BarPlots", "Ratio", "PerMouse", "CellCountsPerMouseAllTypes_LineageGrouped.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(cell_gruop_count, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts", "PerMouse", "CellCountsPerMouseCountsOlig.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(cell_gruop_geno_numb, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts", "CellCountsTimeGenoCountsOlig.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(cell_gruop_geno, filename = file.path(plot_dir, "Cell_Counts","BarPlots", "Ratio", "CellCountsTimeGenoAllTypes_LineageGrouped.png"), width = 7, height = 8, dpi = 300, bg = "white")

ggsave(geno_cell_mouse, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts",  "PerMouse", "CellCountsPerMouseOlig.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(geno_cell_ratio, filename = file.path(plot_dir, "Cell_Counts","BarPlots", "Ratio", "PerMouse", "CellCountsPerMouseOligRatio.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(geno_cell_lineage, filename = file.path(plot_dir, "Cell_Counts","BarPlots", "Ratio", "CellCountsTimeGenoOligRatio.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(geno_cell, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts", "CellCountsTimeGenoOlig.png"), width = 7, height = 8, dpi = 300, bg = "white")

ggsave(geno_cell_mouse_NPC, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts", "PerMouse", "CellCountsPerMouseOligNPC.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(geno_cell_ratio_NPC, filename = file.path(plot_dir, "Cell_Counts","BarPlots", "Ratio", "PerMouse", "CellCountsPerMouseOligNPCRatio.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(geno_cell_lineage_NPC, filename = file.path(plot_dir, "Cell_Counts","BarPlots", "Ratio", "CellCountsTimeGenoOligNPCRatio.png"), width = 7, height = 8, dpi = 300, bg = "white")
ggsave(geno_cell_NPC, filename = file.path(plot_dir, "Cell_Counts", "BarPlots", "RawCounts", "CellCountsTimeGenoNPCOlig.png"), width = 7, height = 8, dpi = 300, bg = "white")

# Save summary
write.csv(expression_summary_data, file = file.path(csv_dir, "Slc35a2_expression_summary.csv"))
