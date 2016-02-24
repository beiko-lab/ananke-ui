library(shinyFiles)

shinyUI(
  navbarPage("timeclust-ui",
    tabPanel("Setup/Load Data",
      fixedPage(
        h3("File Selection"),
        p("Select the appropriate files as indicated"),
        fixedRow(
          column(3,
            shinyFilesButton("timeseries_database", title="Select time series database file", label = "Time-series database file", multiple=FALSE)
          ),
          column(5,
            verbatimTextOutput("tsdb_filename")
          )
        ),
        fixedRow(
          column(3,
            actionButton("loadFiles","Load data!")
          )
        )
      )
    ),
    tabPanel("Data Summary",
      fluidPage(
        titlePanel("Dataset Summary and Statistics"),
        actionButton("plotnclust","Plot Number of Clusters vs Clustering Parameter"),
        plotOutput("clusterVsEps", width= "50%"),
        uiOutput("plotConsistButton"),
        plotOutput("temptaxconsist",width="75%")
      )
    ),
    tabPanel("Explore Clusters",
      tabsetPanel(id="onionTabs",
        tabPanel("Sequences by Time Cluster",value="panel1",
          fluidPage(
            sidebarPanel(
              titlePanel("Explore Time Series Clustering"),
              fluidRow(
                uiOutput("clusterParamSelector")
              ),
              fluidRow(
                 uiOutput("clusterSelector")
              ),
              br(),
              downloadLink("saveTable", "Save current cluster (.csv)"),
              br(),
              downloadLink("saveMainPlot", "Save current plot (.svg)"),
            width=2),                                
            mainPanel(
              fluidRow(
                column(8,
                  plotOutput("main_plot", width = "150%", height="900px")
                )
              ),
              fluidRow(
                column(6,
                  dataTableOutput("infotable")
                )
              )
            )
          )                                        
        ),
        tabPanel("Sequences by OTU", value="panel2",
          fluidPage(
            sidebarPanel(
              titlePanel("Explore OTU Clustering"),
              fluidRow(
                uiOutput("OTUSelector")
              ),
              br(),
              downloadLink("saveOTUTable", "Save current OTU (.csv)"),
              br(),
              downloadLink("saveOTUPlot", "Save current plot (.svg)"),
            width=2),
            mainPanel(
              fluidRow(
                column(8,
                  plotOutput("otu_plot", width="150%", height="900px")
                )
              ),
              fluidRow(
                column(6,
                  dataTableOutput("otutable")
                )
              )
            )
          )
        )
      )
    )
  )
)
