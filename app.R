
library(Seurat)
library(cytokit)
library(shiny)
library(DT)
library(tidyverse)

load("data/seurat_genes.Rda")

# Define UI for cytoscope
ui <- fluidPage(
  
  # Application title
  titlePanel("My dataset"),
  
  # Sidebar with input
  sidebarLayout(
    sidebarPanel(width = 3,
                 h3("Data"),
                 selectInput("sample", "Dataset", multiple = FALSE, selected = "ct_p3",
                             choices = list(
                               "Sample collection 1" = c("Sample 1 name" = "s1_id",
                                                         "Sample 2 name" = "s2_id"),
                               "Sample collection 2" = c("Sample 3 name" = "s3_id",
                                                         "Sample 4 name" = "s4_id",
                                                         "Sample 5 name" = "s5_id")
                             )),
      selectInput("gene", "Genes (max 3)", choices = character(0), multiple = TRUE),
      selectInput("dr", "Dimensionality reduction", multiple = FALSE, choices = c("tsne", "pca"), selected = "tsne"),
      h3("Feature plot"),
      selectInput("label", "Label clusters", multiple = FALSE,
                  choices = list("Long", "Short", "None")),
      selectInput("palette", "Colour palette", choices = list("Grey-red" = "redgrey",
                                                              "Blues" = "blues",
                                                              "Viridis" = "viridis"), selected = "redgrey"),
      selectInput("stat", "Colour cells by", selected = "mean",
                  choices = list("Mean expression" = "mean",
                                 "Percentile bin of expression" = "percentiles")),
      h3("Violin plot"),
      numericInput("point_size", "Point size", value = 0.1, min = 0, max = 1, step = 0.1),
      selectInput("sort", "Sort clusters by expression", choices = c(TRUE, FALSE), selected = FALSE)
    ),
    
    # Output plots
    mainPanel(tabsetPanel(

                tabPanel("Visualize gene expression",
                         plotOutput("dr_plot", width = "5in", height = "5in"),
                         plotOutput("feature"),
                         downloadLink("download_feature", "Download PDF"),
                         plotOutput("vln"),
                         downloadLink("download_vln", "Download PDF")
                ),
                
                tabPanel("Sample info",
                         tableOutput("ncell_table")),
                
                tabPanel("Cluster markers",
                         p("Use the navigation tools to search, filter, and order the table of gene markers."),
                         DT::dataTableOutput("markers"))

              )
    )
  )
)

# Define server logic for cytoscope
server <- function(input, output, session) {

  observe({

    sample <- input$sample

    if (is.null(sample)) sample <- character(0)

    updateSelectInput(session, "gene",
                      choices = genes[[sample]])

  })
  
  # Reactive expressions
  gene   <- reactive({head(input$gene, 3)})
  seurat <- reactive({  get(load(paste0("data/seurat/",  input$sample, ".seurat_small.Rda")))})
  mk     <- reactive({read.delim(paste0("data/markers/", input$sample, ".markers.tsv"),
                                 stringsAsFactors = FALSE)})

  # tSNE
  output$dr_plot <- renderPlot({cytokit::plotDR(seurat(), reduction = input$dr,
                                                colours = seurat()@misc$colours,
                                                title = seurat()@project.name,
                                                point_size = 1.1)})

  # Number of cells
  output$ncell_table <- renderTable({
    n_cells <- as.data.frame(table(seurat()@ident))
    colnames(n_cells) <- c("Cluster", "Number of cells")
    return(n_cells)
  })
  
  # Markers
  output$markers <- DT::renderDataTable({
    mk() %>%
      dplyr::select(cluster, external_gene_name, avg_logFC, p_val_adj, pct.1, pct.2, ensembl_gene_id, description) %>% 
      DT::datatable(filter = "top",
                  rownames = FALSE,
                  selection = "none") %>% 
      formatStyle("cluster",  fontWeight = "bold") %>% 
      formatStyle("external_gene_name",  fontWeight = "bold")
    })

  # Feature plot
  feature_plot <- reactive({cytokit::feature(seurat(), genes = gene(),
                                             palette = input$palette,
                                             reduction = input$dr,
                                             legend = TRUE,
                                             point_size = 1.3,
                                             ncol = length(gene()),
                                             statistic = input$stat,
                                             label = ifelse(input$label == "None", FALSE, TRUE),
                                             label_short = ifelse(input$label == "Short", TRUE, FALSE))})

  output$feature <- renderPlot({feature_plot()})

  # Feature plot download
  output$download_feature <- downloadHandler(filename = "feature.pdf",
                                             content = function(file) {
                                               ggsave(filename = file, plot = feature_plot(), width = 5*length(gene()), height = 4)
                                             })

  # Violin plot
  vln_plot <- reactive({Seurat::VlnPlot(seurat(), features.plot = gene(),
                                        do.sort = input$sort,
                                        point.size.use = input$point_size,
                                        nCol = length(gene()),
                                        cols.use = seurat()@misc$colours, x.lab.rot = TRUE)})

  output$vln <- renderPlot({vln_plot()})

  # Violin plot download
  output$download_vln <- downloadHandler(filename = "violin.pdf",
                                         content = function(file) {
                                           ggsave(filename = file, plot = vln_plot(), width = 7*length(gene()), height = 5)
                                         })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

