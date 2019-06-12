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

# if( !require(anyLlib) )
# {
#   install.packages("anyLib")
#   
#   library(anyLib)
# }

# anyLib::anyLib(c("shiny", "shinydashboard", "shinyWidgets", "DT", "plotly", "ggplot2", "googleVis", "colourpicker"))

# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinyWidgets")
# install.packages("DT")
# install.packages("plotly")
# install.packages("ggplot2")
# install.packages("googleVis")
# install.packages("colourpicker")

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(googleVis)
library(colourpicker)


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
    	                         inline = T),
    	            
    	            actionButton(inputId = "actBtnVisu", label = "Visualization", icon = icon("play"))
	            ),
	            
	            column(9,
	                   h3("File preview"),
      	             dataTableOutput(outputId = "preview")
	            )
	    )),

	    #Visualization
	    tabItem(tabName = "visualization",
	            h1("Data visualization")
	    )	  
	  )
	)
)

server <- function(input, output, session) 
{
  output$preview <- renderDataTable({
    req(input$dataFile)
    
    df <- read.csv(input$dataFile$datapath,
                   header = as.logical(input$header),
                   sep = input$sep,
                   quote = input$quote,
                   nrows = 10)
  }, options = list(scrollX = TRUE, dom = 't'))
  
  data = reactiveValues()
  
  observeEvent(input$actBtnVisu, {
    data$table = read.csv(input$dataFile$datapath,
                          header = as.logical(input$header),
                          sep = input$sep,
                          quote = input$quote,
                          nrows=10)
    sendSweetAlert(
      session = session,
      title = "Done!",
      text = "File was readen succesfully.",
      type = "success"
      )
    
    updateTabItems(session, "tabs", selected = "visualization")
    })

}

shinyApp(ui, server)
