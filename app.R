
library(Seurat)
library(cytokit)
library(shiny)

# Define UI for cytoscope
ui <- fluidPage(
  
  # Application title
  titlePanel("SC dev/PBT dataset"),
  
  # Sidebar with input
  sidebarLayout(
    sidebarPanel(
      h3("Data"),
      selectInput("sample", "Dataset", multiple = FALSE,
                  choices = list(
                    "Mouse forebrain samples" = c("Forebrain E12.5" = "data/ct_e12.seurat_small.Rda",
                                                  "Forebrain E15.5" = "data/ct_e15.seurat_small.Rda",
                                                  "Forebrain P0" = "data/ct_p0.seurat_small.Rda"),
                    "Mouse pons/hindbrain samples" = c("Hindbrain E12.5" = "data/po_e12.seurat_small.Rda",
                                                       "Pons E15.5" = "data/po_e15.seurat_small.Rda",
                                                       "Pons P0" = "data/po_p0.seurat_small.Rda"),
                    "Human fetal brain samples" = c("Human 19pcw brainstem" = "data/hg19w.seurat_small.Rda"),
                    "Re-embedded clusters" = c("Embryonic pons progenitors" = "data/pons_prog.seurat_small.Rda",
                                               "Human brainstem astrocytes" = "data/hg19w_astro.seurat_small.Rda"),
                    "Patient tumor samples" = c("ETMR1" = "data/etmr1.seurat_small.Rda",
                                                "WNT-MB-1" = "data/wnt1.seurat_small.Rda",
                                                "ATRT1" = "data/atrt1.seurat_small.Rda")
                  )),
      textInput("gene", "Gene", placeholder = "Stmn2"),
      selectInput("dr", "Dimensionality reduction", multiple = FALSE, choices = c("tsne", "pca"), selected = "tsne"),
      h3("Feature plot"),
      selectInput("label", "Label clusters", multiple = FALSE,
                  choices = list("Long", "Short", "None")),
      h3("Violin plot"),
      numericInput("point_size", "Point size", value = 0.1, min = 0, max = 1, step = 0.1)
    ),
    
    # Output plots
    mainPanel(width = 5,
              tabsetPanel(
                
                tabPanel("Gene expression",
                         plotOutput("dr_plot"),
                         plotOutput("feature"),
                         downloadLink("download_feature", "Download PDF"),
                         plotOutput("vln"),
                         downloadLink("download_vln", "Download PDF")
                ),
                tabPanel("Sample info",
                         tableOutput("ncell_table"))
                
              )
    )
  )
)

# Define server logic for cytoscope
server <- function(input, output) {
  
  # Reactive expressions
  gene <- reactive({input$gene})
  seurat <- reactive({get(load(input$sample))})
  
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
  
  # Feature plot
  feature_plot <- reactive({cytokit::feature(seurat(), genes = gene(),
                                             reduction = input$dr,
                                             legend = TRUE,
                                             point_size = 1.3,
                                             label = ifelse(input$label == "None", FALSE, TRUE),
                                             label_short = ifelse(input$label == "Short", TRUE, FALSE))})
  
  output$feature <- renderPlot({feature_plot()})
  
  # Feature plot download
  output$download_feature <- downloadHandler(filename = "feature.pdf",
                                             content = function(file) {
                                               ggsave(filename = file, plot = feature_plot(), width = 8, height = 7)
                                             })
  
  # Violin plot
  vln_plot <- reactive({Seurat::VlnPlot(seurat(), features.plot = gene(),
                                        point.size.use = input$point_size,
                                        size.x.use = 12,
                                        cols.use = seurat()@misc$colours, x.lab.rot = TRUE)})
  
  output$vln <- renderPlot({vln_plot()})
  
  # Violin plot download
  output$download_vln <- downloadHandler(filename = "violin.pdf",
                                         content = function(file) {
                                           ggsave(filename = file, plot = vln_plot(), width = 7, height = 5)
                                         })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

