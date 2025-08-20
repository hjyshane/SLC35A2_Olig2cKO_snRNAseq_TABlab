library(Seurat)
library(shiny)
library(tidyverse)

# Define UI for the app
ui <- fluidPage(titlePanel("Interactive Feature Plot"),
                sidebarLayout(
                    sidebarPanel(
                        textInput("gene", "Enter Gene Name:", value = "Bcl6"),
                        actionButton("plotButton", "Generate Feature Plot"),
                        downloadButton("savePlot", "Save Feature Plot")
                    ),
                    mainPanel(plotOutput(
                        "featurePlot", width = "auto", height = "auto"
                    ))
                ))
# Define server logic
server <- function(input, output, session) {
    observeEvent(input$plotButton, {
        output$featurePlot <- renderPlot({
            req(input$gene)
            FeaturePlot(
                object = integrated_man_npc, # Change object here
                features = input$gene,  # Fixed: use the actual gene input
                reduction = "umap",
                label = TRUE,
                repel = TRUE,
                cols = c("lightgrey", "blue")
            ) + ggtitle(paste("Feature Plot for", input$gene)) +
                theme(plot.title = element_text(size = 16, face = "bold"))
        }, height = 800, width = 800)
    })
    output$savePlot <- downloadHandler(
        filename = function() {
            paste(input$gene, "_featureplot.png", sep = "")
        },
        content = function(file) {
            ggsave(
                file,
                plot = last_plot(),
                width = 12,
                height = 12,
                dpi = 300
            )
        }
    )
}

# Run the application
shinyApp(ui = ui, server = server)
