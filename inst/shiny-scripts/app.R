# The code is adapted from
# RStudio Inc. (2013). Tabsets. Shiny Gallery. URL:https://shiny.rstudio.com/gallery/tabsets.html

# I have example rds files for DESeq2 and edgeR result in inst/shiny-scipts

library(shiny)

# Define UI for random distribution app ----
ui <- fluidPage(

  # App title ----
  titlePanel("DEVisComp compMA and compVolcano"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("Plot Venn diagram, MA plots and Volcano plots comparing results
              from DESeq2 and edgeR."),
      tags$p("Upload RDS files containing the result from DESeq2, a DESeqResults
             object, and a data frame taken from the table component of the
             returned value of edgeR topTags function (edgeR::topTags(..)$table
             )"),

      br(),

      # Input files
      fileInput("deseq2File", "Choose RDS File for DESeq2 Result",
                accept = ".rds"
      ),

      fileInput("edgerFile", "Choose RDS File for edgeR Result",
                accept = ".rds"
      ),

      # Input padj cutoff
      textInput(inputId = "padj",
                label = "Enter padj cutoff (0 - 1)", "0.05"),

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output:
      tabsetPanel(type = "tabs",
                  tabPanel("MA plot", plotOutput("ma")),
                  tabPanel("Volcano plot", plotOutput("volcano"))
      )

    )
  )
)

# Define server
server <- function(input, output) {

  # Plotting MA plots
  output$ma <- renderPlot({

    deseq2File <- input$deseq2File
    edgerFile <- input$edgerFile

    # wait until we got all files
    if (is.null(deseq2File) | is.null(edgerFile)) {
      return(NULL)
    }

    deseq2Result <- readRDS(deseq2File$datapath)
    edgerResult <- readRDS(edgerFile$datapath)
    cutoff <- as.numeric(input$padj)

    if (is.null(deseq2Result) | is.null(edgerResult)) {
      return(NULL)
    }

    DEVisComp::compMA(deseq2Result, edgerResult, cutoff)
  })

  # Plotting Volcano plots
  output$volcano <- renderPlot({

    deseq2File <- input$deseq2File
    edgerFile <- input$edgerFile

    # wait until we got all files
    if (is.null(deseq2File) | is.null(edgerFile)) {
      return(NULL)
    }

    deseq2Result <- readRDS(deseq2File$datapath)
    edgerResult <- readRDS(edgerFile$datapath)
    cutoff <- as.numeric(input$padj)

    if (is.null(deseq2Result) | is.null(edgerResult)) {
      return(NULL)
    }

    DEVisComp::compVolcano(deseq2Result, edgerResult, cutoff)
  })
}

# Create Shiny app ----
shinyApp(ui, server)

# [END]
