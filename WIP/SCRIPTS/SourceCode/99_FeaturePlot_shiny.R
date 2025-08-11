# Load required libraries
library(shiny)
library(plotly)
library(dplyr)
library(tibble)
library(Seurat)

# Function to extract plot data from Seurat object
extract_plot_data <- function(seurat_obj, gene, resolution_col = "seurat_clusters") {
    # Extract UMAP coordinates
    umap_coords <- Embeddings(seurat_obj, reduction = "umap")

    # Extract gene expression
    gene_expr <- GetAssayData(seurat_obj, slot = "data")[gene, ]

    # Extract cluster information
    clusters <- seurat_obj@meta.data[[resolution_col]]

    # Combine into a data frame
    plot_data <- data.frame(
        cell_barcode = rownames(umap_coords),
        UMAP_1 = umap_coords[, 1],
        UMAP_2 = umap_coords[, 2],
        gene_expression = gene_expr,
        cluster = clusters,
        stringsAsFactors = FALSE
    )

    return(plot_data)
}

# Function to create interactive plot with highlighting
create_interactive_plot_with_highlighting <- function(plot_data, gene_name, selected_clusters = character(0)) {
    library(plotly)

    # Create highlighting column
    plot_data$highlight_status <- ifelse(
        as.character(plot_data$cluster) %in% selected_clusters,
        "Selected",
        "Unselected"
    )

    # Create the plot with conditional styling
    p <- plot_ly(
        data = plot_data,
        x = ~UMAP_1,
        y = ~UMAP_2,
        color = ~gene_expression,
        colors = c("lightgrey", "blue"),
        type = "scatter",
        mode = "markers",
        # Dynamic marker size based on selection
        marker = list(
            size = ~ifelse(highlight_status == "Selected", 6, 3),
            line = list(
                width = ~ifelse(highlight_status == "Selected", 2, 0),
                color = "red"
            )
        ),
        text = ~paste("Cell:", cell_barcode,
                      "<br>Cluster:", cluster,
                      "<br>Expression:", round(gene_expression, 3),
                      "<br>Status:", highlight_status),
        hovertemplate = "%{text}<extra></extra>",
        source = "main_plot"
    ) %>%
        layout(
            title = paste("Feature Plot for", gene_name,
                          ifelse(length(selected_clusters) > 0,
                                 paste("| Selected clusters:", paste(selected_clusters, collapse = ", ")),
                                 "")),
            xaxis = list(title = "UMAP_1"),
            yaxis = list(title = "UMAP_2")
        )

    return(p)
}

# Function to get cluster information for modal
get_cluster_information <- function(plot_data, selected_clusters, resolution_col) {
    # Get cells in selected clusters
    selected_cells <- plot_data[plot_data$cluster %in% selected_clusters, ]

    # Summary statistics
    summary_table <- data.frame(
        Cluster = selected_clusters,
        Number.of.Cells = sapply(selected_clusters, function(x) {
            sum(plot_data$cluster == x)
        }),
        Avg.Gene.Expression = sapply(selected_clusters, function(x) {
            cluster_cells <- plot_data[plot_data$cluster == x, ]
            round(mean(cluster_cells$gene_expression), 3)
        }),
        stringsAsFactors = FALSE
    )

    # Current annotations (using your predicted.id and ref_short columns)
    # Get metadata for selected cells
    metadata_subset <- integrated@meta.data[selected_cells$cell_barcode,
                                            c("predicted.id", "ref_short", "lineage_group")]
    metadata_subset$cell_barcode <- rownames(metadata_subset)

    # Combine with cluster info
    annotation_data <- merge(selected_cells[, c("cell_barcode", "cluster")],
                             metadata_subset,
                             by = "cell_barcode")

    # Summarize annotations by cluster
    annotation_summary <- annotation_data %>%
        group_by(cluster, predicted.id, ref_short, lineage_group) %>%
        summarise(Count = n(), .groups = "drop") %>%
        arrange(cluster, desc(Count))

    return(list(
        summary_table = summary_table,
        annotation_table = annotation_summary,
        selected_cells = selected_cells
    ))
}

# Function to update Seurat metadata
update_seurat_metadata <- function(seurat_obj, selected_clusters, new_identity, column_name, resolution_col) {
    tryCatch({
        # Get cells in selected clusters
        current_clusters <- seurat_obj@meta.data[[resolution_col]]
        cells_to_update <- which(current_clusters %in% selected_clusters)

        # Create column if it doesn't exist
        if (!column_name %in% colnames(seurat_obj@meta.data)) {
            seurat_obj@meta.data[[column_name]] <- NA
        }

        # Update the metadata
        seurat_obj@meta.data[cells_to_update, column_name] <- new_identity

        # Update the global object (important!)
        assign("integrated", seurat_obj, envir = .GlobalEnv)

        return(TRUE)
    }, error = function(e) {
        print(paste("Error updating metadata:", e$message))
        return(FALSE)
    })
}

# UI
ui <- fluidPage(
    titlePanel("Interactive Feature Plot with Multi-Selection"),
    sidebarLayout(
        sidebarPanel(
            selectInput("resolution", "Choose Clustering Resolution:",
                        choices = c("ref_short",
                                    paste0("RNA_cluster_", c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2))),
                        selected = "ref_short"),

            textInput("gene", "Enter Gene Name:", value = "Bcl6"),
            actionButton("plotButton", "Generate Interactive Plot"),

            # Selection management
            h4("Cluster Selection:"),
            verbatimTextOutput("clickInfo"),

            # Button to clear selections
            actionButton("clearSelection", "Clear All Selections"),

            # Button to proceed to assignment (Phase 4)
            actionButton("assignIdentity", "Assign Identity to Selected Clusters",
                         class = "btn-primary"),

            br(), br(),
            downloadButton("savePlot", "Save Feature Plot")
        ),
        mainPanel(
            plotlyOutput("interactivePlot", width = "100%", height = "800px")
        )
    )
)

# Server
server <- function(input, output, session) {
    # Reactive values for state management
    plot_data <- reactiveVal()
    selected_clusters <- reactiveVal(character(0))

    # Generate plot when button is clicked
    observeEvent(input$plotButton, {
        req(input$gene, input$resolution)

        # Reset selections when new plot is generated
        selected_clusters(character(0))

        # Extract data
        data <- extract_plot_data(integrated, input$gene, input$resolution)
        plot_data(data)

        # Create interactive plot
        output$interactivePlot <- renderPlotly({
            create_interactive_plot_with_highlighting(data, input$gene, selected_clusters())
        })
    })

    # Enhanced click handler with multi-selection
    observeEvent(event_data("plotly_click", source = "main_plot"), {
        click_data <- event_data("plotly_click", source = "main_plot")

        if (!is.null(click_data)) {
            # Get clicked point and its cluster
            point_index <- click_data$pointNumber + 1
            current_data <- plot_data()
            clicked_cluster <- as.character(current_data[point_index, "cluster"])

            # Multi-selection logic (toggle behavior)
            current_selection <- selected_clusters()

            if (clicked_cluster %in% current_selection) {
                # Remove cluster from selection (toggle off)
                new_selection <- current_selection[current_selection != clicked_cluster]
            } else {
                # Add cluster to selection
                new_selection <- c(current_selection, clicked_cluster)
            }

            selected_clusters(new_selection)

            # Update plot with new highlighting
            output$interactivePlot <- renderPlotly({
                create_interactive_plot_with_highlighting(plot_data(), input$gene, new_selection)
            })

            # Display selection info
            output$clickInfo <- renderText({
                if (length(new_selection) == 0) {
                    "No clusters selected. Click clusters to select them."
                } else {
                    paste(
                        "Selected Clusters:", paste(new_selection, collapse = ", "), "\n",
                        "Number of selected clusters:", length(new_selection), "\n",
                        "Click selected clusters again to deselect them."
                    )
                }
            })
        }
    })

    # Clear selection handler
    observeEvent(input$clearSelection, {
        selected_clusters(character(0))

        # Update plot
        if (!is.null(plot_data())) {
            output$interactivePlot <- renderPlotly({
                create_interactive_plot_with_highlighting(plot_data(), input$gene, character(0))
            })

            output$clickInfo <- renderText("All selections cleared.")
        }
    })

    # Modal dialog handler
    observeEvent(input$assignIdentity, {
        req(length(selected_clusters()) > 0)

        # Get information about selected clusters
        current_data <- plot_data()
        current_resolution <- input$resolution

        # Calculate cluster information
        cluster_info <- get_cluster_information(current_data, selected_clusters(), current_resolution)

        showModal(modalDialog(
            title = "Assign Cell Type Identity",
            size = "l",

            # Cluster information display
            h4("Selected Cluster Information:"),
            tableOutput("clusterInfoTable"),

            br(),

            # Current annotations in these clusters
            h4("Current Cell Type Predictions:"),
            tableOutput("currentAnnotations"),

            br(),

            # Input fields
            fluidRow(
                column(6,
                       textInput("newIdentity", "New Cell Type Identity:",
                                 placeholder = "e.g., Oligodendrocytes, Astrocytes")
                ),
                column(6,
                       textInput("columnName", "Metadata Column Name:",
                                 value = paste0("manual_annotation_", Sys.Date()),
                                 placeholder = "e.g., manual_annotation_v1")
                )
            ),

            # Action buttons
            footer = tagList(
                actionButton("cancelAssignment", "Cancel", class = "btn-secondary"),
                actionButton("confirmAssignment", "Save Assignment", class = "btn-primary")
            )
        ))

        # Render cluster information table
        output$clusterInfoTable <- renderTable({
            cluster_info$summary_table
        })

        # Render current annotations table
        output$currentAnnotations <- renderTable({
            cluster_info$annotation_table
        })
    })

    # Handle assignment confirmation
    observeEvent(input$confirmAssignment, {
        req(input$newIdentity, input$columnName)

        # Update Seurat object metadata
        success <- update_seurat_metadata(
            seurat_obj = integrated,
            selected_clusters = selected_clusters(),
            new_identity = input$newIdentity,
            column_name = input$columnName,
            resolution_col = input$resolution
        )

        if (success) {
            # Show success message
            showNotification(
                paste("Successfully assigned", input$newIdentity,
                      "to", length(selected_clusters()), "clusters"),
                type = "message",
                duration = 5
            )

            # Clear selections
            selected_clusters(character(0))

            # Update plot
            output$interactivePlot <- renderPlotly({
                create_interactive_plot_with_highlighting(plot_data(), input$gene, character(0))
            })

            output$clickInfo <- renderText("Assignment completed. Selections cleared.")
        } else {
            showNotification("Error updating metadata", type = "error")
        }

        removeModal()
    })

    # Handle assignment cancellation
    observeEvent(input$cancelAssignment, {
        removeModal()
    })
}

# Run the application
shinyApp(ui = ui, server = server)