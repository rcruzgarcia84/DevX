### Shiny Dendrometer curves app!

#install.packages("shiny")
library(shiny)
library(tidyverse)
library(broom)
library(lubridate)
source("auxiliary_functions.R", local = T)
###

### Load data

dendrom_curves <- read_csv("dendro_curves.csv")

dendrom_curves <- dendrom_curves %>% mutate(step_locator = seq_along(tiempo))

dendrom_curves <- dendrom_curves %>% select(tiempo, Dendrometer, tree, step_locator) %>% na.omit() %>% 
  group_by(tree) %>% mutate(min_max_norm = min_max_norm(Dendrometer), 
                            sri = Dendrometer - min(Dendrometer)) %>% 
  ungroup() 


### Bind Image Data and Dendrometer Data

image_dates <- read_csv("Microcore_sampling_dates.csv")

image_dates <- image_dates %>% gather(key = "date", value = "image",  -Baum, -ID_RCG) %>% unite("tree", Baum, ID_RCG, sep = "")

image_dates$date <- as.POSIXct(image_dates$date, format = "%d.%m.%y", tz = "UTC")


### Bind with DendrometerData

### Create the day_wimage object

days_wimage <- image_dates %>% filter(!is.na(image)) 

rm(image_dates)

### Add a column to the dataframe with the path or name of the image file for the date!

### Load image files, or names at least 

bilder <- list.files("./images")


### you create a list of the files

##then you compare which trees match (substr) and which dates (complex posixct and substr)
### wrap around which and you have the rows where the files are supposed to go!

#step_1 <- substr(bilder, 1, 5) %>% str_replace_all("_", "") %>% 
# str_replace("1", "A") %>% str_replace("2", "B") %>% str_replace("3", "C")
step_1 <- substr(bilder, 1, 4)

step_2 <- as.POSIXct(substr(bilder, 6, nchar(bilder)-8), format = "%d%m%Y", tz = "UTC")

pos_bilder <- match(paste(step_1, step_2), paste(days_wimage$tree, days_wimage$date))

days_wimage[pos_bilder, "file"] <- bilder ### place the filenames in the right position

rm(bilder, pos_bilder, step_1, step_2)

file_wimage <- days_wimage %>% select(tree, date, file) %>% mutate(dates = as.Date(date))### extract only important columns

colnames(file_wimage) <- c("tree", "date_posixct", "file", "date")

dendrom_curves$date <- as.Date(dendrom_curves$tiempo)

dendrom_curves <- left_join(dendrom_curves, file_wimage, by = c("tree", "date")) ### left joining the dataframes

### Add point dates for images!

point_dates <- dendrom_curves$date[!is.na(dendrom_curves$file)]### get the dates for 
### points where we have images

row_with_date <- which(!is.na(dendrom_curves$file))## get the row positions

length(point_dates);length(row_with_date)### same length

dendrom_curves[row_with_date, "x_images"] <- as.POSIXct(point_dates, format = "%Y-%m-%d", tz = "UTC")### save it as character

dendrom_curves$x_images <- as.POSIXct(dendrom_curves$x_images, origin = "1970-01-01")


### Add a clear Y value for points with image

y_values <- dendrom_curves %>% group_by(tree, file) %>% summarize(y_value_points_raw  = mean(Dendrometer, na.rm = T), 
                                                                  y_value_points_norm  = mean(min_max_norm, na.rm = T), 
                                                                  y_value_points_sri  = mean(sri, na.rm = T)) %>%  filter(file != "NA")


dendrom_curves <- left_join(dendrom_curves, y_values, by = "file")
### after here is tree.x

rm(point_dates, row_with_date, file_wimage, y_values)



### Prepare data, modelling

dendrom_curves_models <- dendrom_curves %>% select(tiempo, Dendrometer, tree.x, 
                                                   step_locator, min_max_norm, sri) %>% na.omit() %>% 
  group_by(tree.x) %>% mutate(steps = seq_along(tiempo)) %>% 
  ungroup() ## select columns of interest
### and add a "steps" column for each tree, to model in base of that instead of the POSIXct data

### Gompertz Models

dendrom_curves_models_g <- dendrom_curves_models %>% group_by(tree.x) %>% nest() %>%
  mutate(gompertz = map(data, ~ nls_rcg(.)), tidied_gompertz = map(gompertz, tidy), gompertz_fitted  = map(gompertz, augment)) %>%
  unnest(gompertz_fitted)

### Weibull


### Primero get initials pa todos

initials_weibull <- vector(mode = "list", length = length(unique(dendrom_curves_models$tree.x)))
names(initials_weibull) <- unique(dendrom_curves_models$tree.x)

for(i in unique(dendrom_curves_models$tree.x)){
  pos <- which(dendrom_curves_models$tree.x == i)
  initials_weibull[[i]]  <- tryCatch(getInitial(min_max_norm ~ SSweibull(steps, Asym, Drop, lrc, pwr),
                                                data = dendrom_curves_models[pos,]), 
                                     error = function(e) paste("error"))
}


errors <- which(initials_weibull == "error")
replace_par <- which(initials_weibull != "error")[1]

initials_weibull <- replace(initials_weibull, errors, initials_weibull[replace_par])


dendrom_curves_models_w <- dendrom_curves_models %>% group_by(tree.x) %>% nest()  %>%
  mutate(weibull = map2(.y = .$data, .x = initials_weibull,  ~ nls(weibull_formula, data = .y, start = unlist(.x))), 
         tidied_weibull = map(weibull, tidy), weibull_fitted  = map(weibull, augment)) %>%
  unnest(weibull_fitted)

### Putting both model and Raw Measurements together

dendrom_curves_models_fertig <- bind_cols(dendrom_curves_models_w, dendrom_curves_models_g, dendrom_curves_models[, c("Dendrometer", "sri", "step_locator")]) %>%
  select(tree.x, min_max_norm, steps, "weibull_fit" = .fitted, 
         "weibull_res" = .resid, "gompertz_fit"= .fitted1, 
         "gompertz_res" = .resid1, "raw_data" = Dendrometer, 
         "SRI" = sri, step_locator)

pos_steps <- match(dendrom_curves_models_fertig$steps,dendrom_curves_models$steps) ## getting back the tiempo column, get matching positions with orignal table

dendrom_curves_models_fertig[,"tiempo"] <- as.vector(dendrom_curves_models[pos_steps, "tiempo"]) ### put them back carefully (tibbles are messy)

### Joining Photograph info and data for Points in graph

pos <- match(dendrom_curves_models_fertig$step_locator,dendrom_curves$step_locator)

dendrom_curves_models_fertig <- bind_cols(dendrom_curves_models_fertig, dendrom_curves[pos, c("date", "file", "x_images", "y_value_points_norm", "y_value_points_sri", "y_value_points_sri")])

rm(nls_rcg, dendrom_curves, pos_steps, gompertz_initials, min_max_norm, dendrom_curves_models, dendrom_curves_models_g, dendrom_curves_models_w, 
   initials_weibull, errors, gompertz_formula, i, pos, replace_par, weibull_formula)

### Getting Growth phenology


growth_phen_gompertz <- growth_phenology_calc(dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = gompertz_fit)

growth_phen_weibull <- growth_phenology_calc(dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = weibull_fit)

growth_phen_change <- growth_phenology_calc(dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = SRI)

rm(growth_phenology_calc)



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
      plot_true <- dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
        ggplot(aes(color = tree.x)) + geom_line(aes(x = tiempo, y = min_max_norm)) +
        geom_point(aes(x = x_images, y = y_value_points_norm, colour = tree.x), 
                   size = 4, alpha = 0.75)+
        geom_vline(data = growth_phen_weibull %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_onset, color = tree.x)) +
        geom_vline(data = growth_phen_weibull %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_cessation, color = tree.x))+
        # scale_x_datetime(breaks = date, labels = date) +
        theme(panel.background = element_rect(fill = "white", color = "black"), 
              axis.text.x = element_text(angle = 90, hjust = 1)) 
      if(input$weibull == TRUE){
        plot_true <- plot_true + geom_line(aes(x = tiempo, y = weibull_fit, group = tree.x), color = "black", alpha = 0.7, data = dendrom_curves_models_fertig %>%
                                             filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
      } else {
        plot_true
      }
      plot_true
      
    } else {
      plot_false <- dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
        ggplot(aes(color = tree.x)) + geom_line(aes(x = tiempo, y = min_max_norm)) + 
        geom_point(aes(x = x_images, y = y_value_points_norm, colour = tree.x), 
                   size = 4, alpha = 0.75) + #scale_x_datetime(breaks = date, labels = date) +
        theme(panel.background = element_rect(fill = "white", color = "black"), 
              axis.text.x = element_text(angle = 90, hjust = 1))
      if (input$weibull == TRUE){
        plot_false <- plot_false + geom_line(aes(x = tiempo, y = weibull_fit, group = tree.x), color = "black", alpha = 0.7, data = dendrom_curves_models_fertig %>%
                                               filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
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
                                          as.character(unique(na.omit(nearPoints(dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
                                                                                 input$plot_click, xvar = "x_images", yvar = "y_value_points_norm", 
                                                                                 threshold = 3000000, maxpoints = 2, addDist = T)[,13])[1])[1])))
      filename <- na.omit(filename)[[1]]
      # nearPoints() also works with hover and dblclick events
      list(src = filename,
           contentType = "image/jpg",
           width = 500,
           height = 400,
           alt = "Image file")
    } else {
      filename <- normalizePath(file.path( './Annot_Images',
                                           as.character(unique(na.omit(nearPoints(dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
                                                                                  input$plot_click, xvar = "x_images", yvar = "y_value_points_norm", 
                                                                                  threshold = 3000000, maxpoints = 2, addDist = T)[,13])[1])[1])))
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
                                        as.character(unique(na.omit(nearPoints(dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
                                                                               input$plot_click, xvar = "x_images", yvar = "y_value_points_norm", 
                                                                               threshold = 3000000, maxpoints = 2, addDist = T)[,13])[1])[1])))
    #substr(filename, nchar(filename)-21, nchar(filename))
    filename
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)

