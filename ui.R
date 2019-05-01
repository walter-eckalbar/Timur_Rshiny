
library(shinythemes)
library(shiny)

options(shiny.maxRequestSize=1000*1024^2)

ui <- fluidPage(
  titlePanel("RNAseqPlotting", windowTitle = "RNA-seq plotting"),
  
  navbarPage("", theme = shinytheme("spacelab"), id = "plotType",
             tabPanel("Gene Plots",
                      titlePanel(h3("Individual Gene Plots")),
                      ##### LOOK UP GENES SIDEBAR ####
                      sidebarPanel(h4("Plot Options"),
                                   hr(),
                                   fluidRow(
                                     column(width = 4, offset = 0,
                                            radioButtons("nameType", label = h5("Name Type"), 
                                                         choices = list("Ensembl ID" = "Ensembl_ID", "Gene Name" = "Gene"), selected = "Gene")
                                     )
                                   ),
                                   fluidRow(
                                     selectInput("geneLookUpX", label = "Enter condition or gene for the X-axis", c("this with update")),
                                                    br(),
                                                    helpText("This input is case sensitive")
                                   ),
                                   fluidRow(
                                     selectInput("geneLookUpY", label = "Enter gene for the Y-axis", c("this with update")),
                                                    br(),
                                                    helpText("This input is case sensitive")
                                   )
                      ),
                      mainPanel(
                                  fluidRow(
                                    column(4, downloadButton('downloadBoxPlot', 'Download Box Plot'),
                                                        p("Download your box plot image file"))
                                  ),
                                  fluidRow(
                                    conditionalPanel(condition = "length(output.boxText > '0'",
                                                          textOutput("boxText")
                                    ),
                                    uiOutput("plotBox.ui")
                                  ),
                                  fluidRow(
                                    dataTableOutput('tableSingleBox')
                                  )
                      )
             )
  )
)
                                   
