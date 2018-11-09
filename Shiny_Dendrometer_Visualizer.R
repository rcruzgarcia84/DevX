### Shiny Dendrometer curves app!

#install.packages("shiny")
library(shiny)
library(tidyverse)


# Define UI ----
ui <- fluidPage(
  
  titlePanel("Dendrometer Visualizer"),
  
  fluidRow(
    
    column(1, wellPanel(
      h5("Select a tree"), 
      checkboxGroupInput(inputId = "tree", 
                         label = NULL,
                         choices = list("FSY A" = "FSYA", 
                                        "FSY B" = "FSYB", 
                                        "FSY C" = "FSYC", 
                                        "QRO A" = "QROA", 
                                        "QRO B" = "QROB", 
                                        "QRO C" = "QROC", 
                                        "APS A" = "APSA", 
                                        "APS B" = "APSB"), 
                         selected = "FSYA"))),
    column(2, wellPanel(dateRangeInput("date", 
                                       label = "Choose a date range within 2016", 
                                       start = "2016-03-01", end = "2016-09-30", 
                                       min = "2016-03-01", max = "2016-09-30" )), 
           wellPanel(radioButtons("radio",
                                  label = "Show derived phenology?", 
                                  choices = list("Yes" = TRUE, "No" = FALSE)))),
    
    column(5, imageOutput("image"), 
           verbatimTextOutput("name_image"))
  ),
  
  hr(),
  
  plotOutput(outputId = "dendro_curve", click = clickOpts(id = "plot_click"))
)


# Define server logic ----
server <- function(input, output) {
  dendro_curves <- read_csv("/home/mochomo/Doktorarbeit/Microcore vs Dendrometer/Dendro_Visualizer/dendrom_curves_tidy.csv")
  
  days_wimage <- read.csv("/home/mochomo/Doktorarbeit/Microcore vs Dendrometer/days_wimage.csv", 
                          header = T, stringsAsFactors = F, colClasses = c("character", "POSIXct", "character", "POSIXct", 
                                                                           "POSIXct", "numeric", "numeric", "character"))
  
  growth_phenology <- read_csv("/home/mochomo/Doktorarbeit/Microcore vs Dendrometer/growth_phenology.csv")
  
  
  output$dendro_curve <- renderPlot({
    
    tree_of_choice <- input$tree
    
    date_begin <- input$date[1]
    date_end <- input$date[2]
    
    if(input$radio == T){
      dendro_curves %>% filter(tree %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
        ggplot(aes(color = tree)) + geom_line(aes(x = tiempo, y = Dendrometer)) +
        geom_point(data = days_wimage %>% filter(tree %in% tree_of_choice, date >= date_begin & date <= date_end), aes(x = date, y = ymax-50, colour = tree), 
                   size = 4, alpha = 0.75, inherit.aes = F)+
        geom_vline(data = growth_phenology %>% filter(tree %in% tree_of_choice), aes(xintercept = date_onset, color = tree)) +
        geom_vline(data = growth_phenology %>% filter(tree %in% tree_of_choice), aes(xintercept = date_cessation, color = tree))+
        scale_alpha(guide = F) + scale_x_datetime(breaks = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d"), labels = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d")) +
        theme(panel.background = element_rect(fill = "white", color = "black"), 
              axis.text.x = element_text(angle = 90, hjust = 1)) 
    } else {
      dendro_curves %>% filter(tree %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
        ggplot(aes(color = tree)) + geom_line(aes(x = tiempo, y = Dendrometer)) +
        geom_point(data = days_wimage %>% filter(tree %in% tree_of_choice, date >= date_begin & date <= date_end), aes(x = date, y = ymax-50, colour = tree), 
                   size = 4, alpha = 0.75, inherit.aes = F)+
        scale_alpha(guide = F) + scale_x_datetime(breaks = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d"), labels = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d")) +
        theme(panel.background = element_rect(fill = "white", color = "black"), 
              axis.text.x = element_text(angle = 90, hjust = 1))}
    
    
  })
  
  output$image <- renderImage({
    tree_of_choice <- input$tree
    # With base graphics, need to tell it what the x and y variables are.
    filename <- normalizePath(file.path('/home/mochomo/Doktorarbeit/Microcore vs Dendrometer/images',
                                        as.character(nearPoints(dendro_curves %>% filter(tree %in% tree_of_choice),
                                                                input$plot_click, xvar = "tiempo", yvar = "Dendrometer",  threshold = 100000, maxpoints = 1, addDist = T)[,5])))
    filename <- na.omit(filename)[[1]]
    # nearPoints() also works with hover and dblclick events
    list(src = filename,
         contentType = "image/jpg",
         width = 500,
         height = 400,
         alt = "Image file")
    
  }, deleteFile = F)
  
  output$name_image <- renderText({
    tree_of_choice <- input$tree
    # With base graphics, need to tell it what the x and y variables are.
    filename <- normalizePath(file.path('/home/mochomo/Doktorarbeit/Microcore vs Dendrometer/images',
                                        as.character(nearPoints(dendro_curves %>% filter(tree %in% tree_of_choice),
                                                                input$plot_click, xvar = "tiempo", yvar = "Dendrometer", threshold = 100000, maxpoints = 1, addDist = T)[,5])))
    substr(filename, 60, nchar(filename))
                              })
 
}

# Run the app ----
shinyApp(ui = ui, server = server)

