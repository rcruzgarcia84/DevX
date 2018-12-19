### Shiny Dendrometer curves app!

#install.packages("shiny")
library(shiny)
library(tidyverse)
library(broom)
library(lubridate)

### Load data
dendrom_curves <- read_csv("dendrom_curves_tidy.csv")


### ### Using Weibull function
### Fit a curve to each Dendrometer curve
### and get growth phenology (5 and 95% amount of change for growth begin and cessation)

dendrom_curves_models <- dendrom_curves %>% select(tiempo, Dendrometer, tree) %>% group_by(tree) %>% mutate(steps = seq_along(tiempo)) %>% ungroup() ## select columns of interest
### and add a "steps" column for each tree, to model in base of that instead of the POSIXct data

dendrom_curves_models_1 <- dendrom_curves_models %>% na.omit() %>% group_by(tree) %>% nest() %>% 
  mutate(weibull = map(data, ~ nls(Dendrometer ~ SSweibull(steps, Asym, Drop, lrc, pwr), data = .)), tidied_weibull = map(weibull, tidy), weibull_fitted  = map(weibull, augment)) %>%
  unnest(weibull_fitted) ###  apply the weibull models to each courve and get the fitted values with augment()

pos_steps <- match(dendrom_curves_models_1$steps,dendrom_curves_models$steps) ## getting back the tiempo column, get matching positions with orignal table

dendrom_curves_models_1[,"tiempo"] <- as.vector(dendrom_curves_models[pos_steps, "tiempo"]) ### put them back carefully (tibbles are messy)


### get the growth phases (onset, duration and cessation)

growth_phenology <- dendrom_curves_models_1 %>% group_by(tree) %>% 
  summarize(max(.fitted), min(.fitted), total_amplitude = max(.fitted) - min(.fitted), five_percent = total_amplitude * 0.05, ninety_five = total_amplitude *0.95, begin = max(.fitted) - ninety_five, 
            cessation = max(.fitted) - five_percent, date_onset = tiempo[which.min(abs(.fitted - begin))], date_cessation = tiempo[which.min(abs(.fitted - cessation))], 
            growth_duration = date_cessation - date_onset, doy_begin = yday(date_onset), doy_cessation = yday(date_cessation)) ##using which.min(abs("values" - "value to query")) you get the position of the closes value!

rm(dendrom_curves_models)


# Define UI ----
ui <- fluidPage(
  
  titlePanel(h1(strong("Dendrometer Visualizer"))),
  
  hr(),
  
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
                                  choices = list("Yes" = TRUE, "No" = FALSE))),
           wellPanel(radioButtons("image_anot",
                                  label = "Show annotated image?", 
                                  choices = list("Yes" = TRUE, "No" = FALSE), selected = FALSE))
           ),
    column(2, wellPanel(radioButtons("weibull", label = "Show Weibull fit?", 
                        choices = list("Yes" = TRUE, "No" = FALSE), selected = FALSE))),
    
    column(3, imageOutput("image", height = "500px"), 
           textOutput("name_image")), 
    column(2, p("This app helps visualize Dendrometer Data and Microcoring from the same tree/time period. On the left you can control the options to Display and below you can see the Dendro
                meter graph. If you click on the points on the graph, you can see a thin-section of the wood on that sampling date and compare it to derived phenology."),
           br(), 
           div(h3(strong("Click on the dots to see the corresponding thin-section!"), style = "color:blue")),
            offset = 1)
  ),
  
  hr(),
  
  plotOutput(outputId = "dendro_curve", click = clickOpts(id = "plot_click"))
)


# Define server logic ----
server <- function(input, output) {
  
   output$dendro_curve <- renderPlot({
    
    tree_of_choice <- input$tree
    
    date_begin <- input$date[1]
    date_end <- input$date[2]
    
    if(input$radio == T){
      plot_true <- dendrom_curves %>% filter(tree %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
        ggplot(aes(color = tree)) + geom_line(aes(x = tiempo, y = Dendrometer)) +
        geom_point(aes(x = x_images, y = y_value_points + 5, colour = tree), 
                   size = 4, alpha = 0.75, inherit.aes = T)+
        geom_vline(data = growth_phenology %>% filter(tree %in% tree_of_choice), aes(xintercept = date_onset, color = tree)) +
        geom_vline(data = growth_phenology %>% filter(tree %in% tree_of_choice), aes(xintercept = date_cessation, color = tree))+
        scale_alpha(guide = F) + #scale_x_datetime(breaks = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d"), labels = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d")) +
        theme(panel.background = element_rect(fill = "white", color = "black"), 
              axis.text.x = element_text(angle = 90, hjust = 1)) 
      if(input$weibull == TRUE){
        plot_true <- plot_true + geom_line(aes(x = tiempo, y = .fitted, group = tree), color = "black", alpha = 0.7, data = dendrom_curves_models_1 %>%
                    filter(tree %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
        } else {
        plot_true
      }
      plot_true
      
    } else {
      plot_false <- dendrom_curves %>% filter(tree %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
        ggplot(aes(color = tree)) + geom_line(aes(x = tiempo, y = Dendrometer)) + 
        geom_point(aes(x = x_images, y = y_value_points + 5, colour = tree), 
                   size = 4, alpha = 0.75, inherit.aes = T)+
        scale_alpha(guide = F) + #scale_x_datetime(breaks = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d"), labels = as.POSIXct(unique(days_wimage$date), "%Y-%m-%d")) +
        theme(panel.background = element_rect(fill = "white", color = "black"), 
              axis.text.x = element_text(angle = 90, hjust = 1))
      if (input$weibull == TRUE){
        plot_false <- plot_false + geom_line(aes(x = tiempo, y = .fitted, group = tree), color = "black", alpha = 0.7, data = dendrom_curves_models_1 %>%
                                             filter(tree %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
        } else {
          plot_false
        }
      
      plot_false}
    
    
  })
  
  output$image <- renderImage({
    tree_of_choice <- input$tree
    # With base graphics, need to tell it what the x and y variables are.
    if(input$image_anot != TRUE){
      filename <- normalizePath(file.path('./images',
                                        as.character(unique(na.omit(nearPoints(dendrom_curves %>% filter(tree %in% tree_of_choice),
                                                                               input$plot_click, xvar = "tiempo", yvar = "Dendrometer", 
                                                                               threshold = 3000000, maxpoints = 5, addDist = T)[,5])[1])[1])))
    filename <- na.omit(filename)[[1]]
    # nearPoints() also works with hover and dblclick events
    list(src = filename,
         contentType = "image/jpg",
         width = 500,
         height = 400,
         alt = "Image file")
    } else {
      filename <- normalizePath(file.path( './Annot_Images',
                                          as.character(unique(na.omit(nearPoints(dendrom_curves %>% filter(tree %in% tree_of_choice),
                                                                                 input$plot_click, xvar = "tiempo", yvar = "Dendrometer", 
                                                                                 threshold = 3000000, maxpoints = 5, addDist = T)[,5])[1])[1])))
      filename <- na.omit(filename)[[1]]
      # nearPoints() also works with hover and dblclick events
      list(src = filename,
           contentType = "image/jpg",
           width = 500,
           height = 400,
           alt = "Image file")
    }
    
  }, deleteFile = F)
  
  output$name_image <- renderText({
    tree_of_choice <- input$tree
    # With base graphics, need to tell it what the x and y variables are.
    filename <- normalizePath(file.path('./images',
                                        as.character(unique(na.omit(nearPoints(dendrom_curves %>% filter(tree %in% tree_of_choice),
                                                                               input$plot_click, xvar = "tiempo", yvar = "Dendrometer", 
                                                                               threshold = 3000000, maxpoints = 5, addDist = T)[,5])[1])[1])))
    substr(filename, nchar(filename)-21, nchar(filename))
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
