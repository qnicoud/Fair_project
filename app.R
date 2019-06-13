########################################################
# Fair_bioinfo - Shiny
# Quentin Nicoud
# 06/2019
# GitHub :
########################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
########################################################
# Library
########################################################

rm(list = ls())

if( !require(anyLlib) )
{
  install.packages("anyLib")

  library(anyLib)
}

anyLib::anyLib(c("shiny", "shinydashboard", "shinyWidgets", "DT", "plotly", "ggplot2", "googleVis", "colourpicker"))

# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinyWidgets")
# install.packages("DT")
# install.packages("plotly")
# install.packages("ggplot2")
# install.packages("googleVis")
# install.packages("colourpicker")
# 
# library(shiny)
# library(shinydashboard)
# library(shinyWidgets)
# library(googleVis)
# library(colourpicker)


########################################################
# Shiny app
########################################################


ui <- dashboardPage(
	dashboardHeader(title = "Fair_bioinfo App"),
	dashboardSidebar(
	  sidebarMenu(
	    menuItem("Import data", tabName = "readData", icon = icon("readme")),
	    menuItem("Data visualization", tabName = "visualization", icon = icon("poll"))
	  )
	),
	dashboardBody(
	  tabItems(
	    #Import data
	    tabItem(tabName = "readData",
	            h1("Import data"),
	            fileInput("dataFile", label = NULL,
	                      buttonLabel = "Browse ...",
	                      placeholder = "No File Selected"),
      fluidRow(           
	            column(3,
    	            h3("Parameters"),
    	            
    	            #Input: Checkbox if header
    	            radioButtons(inputId = "header",
    	                         label = "Header",
    	                         choices = c("Yes" = TRUE,
    	                                     "No" = FALSE),
    	                         selected = TRUE, 
    	                         inline = T),
    	            
    	            #Input: select separator
    	            radioButtons(inputId = "sep",
    	                         label = "Separator",
    	                         choices = c(Comma = ",",
    	                                     Semicolon = ".",
    	                                     Tab = "\t"),
    	                         selected = "\t",
    	                         inline = T),
    	            
    	            #Input: Select quotes
    	            radioButtons(inputId = "quotes",
    	                         label = "Quotes",
    	                         choices = c(None = "",
    	                                     "Single quote" = "'",
    	                                     "Double quote" = '"'),
    	                         selected = "",
    	                         inline = T)
	            ),
	            
	            column(9,
	                   h3("File preview"),
      	             dataTableOutput(outputId = "preview")
	            )
	    ),
      
      tags$br(),
      
      div(actionButton(inputId = "actBtnVisu", label = "Visualization",icon = icon("play") ), align = "center")
      
      ),

	    #Visualization
	    tabItem(tabName = "visualization",
	            h1("Data visualization"),
	            h2("Table exploration"),
	            dataTableOutput('dataTable'),
	            h2("Graphics"),
	            fluidRow(
	              column(8,plotOutput("plotAvecR")),
	              column(4,
	                     colourpicker::colourInput("colR", "Couleur graphique R", "black",allowTransparent = T),
	                     sliderInput("cex", "Taille",
	                                 min = 0.5, max = 3,
	                                 value = 1,step = 0.2)
	                     ),
	              column(4, 
	                     selectInput(inputId = "pch", choices = 1:20, label = "Type de points",selected = 1),
	                     textInput("title", "Titre", "Sepal length vs Petal length (R)")
	                     )
	            )
	    )	  
	  )
	)
)

################################################################################
# Server
################################################################################

server <- function(input, output, session) 
{
  data = reactiveValues()
  
  #=============================================================================
  # Preview
  #=============================================================================
  
  output$preview <- renderDataTable({
    req(input$dataFile)
    
    df <- read.csv(input$dataFile$datapath,
                   header = as.logical(input$header),
                   sep = input$sep,
                   quote = input$quote,
                   nrows = 10)
  }, options = list(scrollX = TRUE, dom = 't'))
  
  #=============================================================================
  # Lecture
  #=============================================================================
  
  observeEvent(input$actBtnVisu, {
    
    if(!is.null(input$dataFile$datapath)) {
      data$table = read.csv(input$dataFile$datapath,
                            header = as.logical(input$header),
                            sep = input$sep,
                            quote = input$quote)
      sendSweetAlert(
        session = session,
        title = "Done!",
        text = "File was readen succesfully.",
        type = "success"
      )
      
      updateTabItems(session, "tabs", selected = "visualization")
    }
    })
  
  #=============================================================================
  # Exploration du tableau
  #=============================================================================
  
  output$dataTable = DT::renderDataTable({
    datatable(data$table, filter = 'top') %>%
      formatStyle('Sepal.Length',
                  background = styleColorBar(data$table$Sepal.Length, 'lightcoral'),
                  backgroudSize = '100% 90%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center'
      ) %>%
      formatStyle('Sepal.Width',
                  backgroundColor = styleInterval(c(3,4), c('white', 'red', "firebrick")),
                  color = styleInterval(c(3,4), c('black', 'white', "white"))
      ) %>%
      formatStyle('Petal.Length',
                  background = styleColorBar(data$table$Petal.Length, 'lightcoral'),
                  backgroundSize = '100% 90%',
                  backgroundRepeat = 'no-repeat',
                  backgroundPosition = 'center'
      ) %>%
      formatStyle('Petal.Width',
                  backgroundColor = styleInterval(c(1,2), c('white', 'red', "firebrick")),
                  color = styleInterval(c(1,2), c('black', 'white', "white"))
      ) %>%
      formatStyle('Species',
                  backgroundColor = styleEqual(
                      levels = unique(data$table$Species), values = c('lightblue', 'lightgreen', 'lavender')
                  )
      )
      
  })
  
  #=============================================================================
  # Graphiques
  #=============================================================================
  # dataTable_tab_rows_selected
  output$plotAvecR <- renderPlot({
    if(!is.null(data$table)) {
      plot(data$table$Petal.Length[input$dataTable_rows_all],
           data$table$Sepal.Length[input$dataTable_rows_all],
           main = input$title,
           ylab = "Sepal length",
           xlab = "Petal length",
           pch = as.numeric(input$pch),
           col = input$colR, 
           cex = input$cex)
    } else {
      NULL
    }
  })

}

shinyApp(ui, server)
