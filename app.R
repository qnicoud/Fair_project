########################################################
# Fair_bioinfo - Shiny
# Quentin Nicoud
# 06/2019
# GitHub :
########################################################

########################################################
# Library
########################################################

rm(list = ls())

# install.packages("shiny")
install.packages("shinydashboard")
install.packages("shinyWidgets")
# install.packages("DT")
# install.packages("plotly")
# install.packages("ggplot2")
# install.packages("googleVis")
# install.packages("colourpicker")

library(shiny)
library(shinydashboard)


########################################################
# Shiny app
########################################################


ui <- dashboardPage(
	dashboardHeader(title = "Fair_bioinfo App"),
	dashboardSidebar(),
	dashboardBody()
)

server <- function(input, output) {}

shinyApp(ui, server)
