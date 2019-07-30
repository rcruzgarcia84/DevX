### Shiny Dendrometer curves app!


#install.packages("shinydashboard")
library(shinydashboard)
library(shiny)
library(tidyverse)
library(broom)
library(lubridate)
source("auxiliary_functions.R", local = T)

###



#### Graphical Interface

ui <- dashboardPage(
  dashboardHeader(title = "DevX"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dendro/Thinsection\nViewer", tabName = "dendro", icon = icon("area-chart")), 
      menuItem("Phenology Comparison", tabName = "pheno", icon = icon("envira")), 
      menuItem("Source Code", icon = icon("code"), href = "https://github.com/rcruzgarcia84/DevX/blob/master/DevX_v0.1.R"),
      menuItem("Upload User Data", icon = icon("upload"), tabName = "upload")
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem("dendro",
              fluidRow(
                box(solidHeader = T, width = 12, plotOutput(outputId = "dendro_curve", click = clickOpts(id = "plot_click")), 
                    status = "primary")
              ),
              fluidRow(
                box(solidHeader = T, width = 4, title = "Image", column(width = 12, imageOutput("image", height = "100%", width = "100%")), 
                    textOutput("name_image"), 
                    status = "primary"), 
                #style = "height:50vh"), 
                tabBox(width = 4, 
                       title = "Controls",
                       side = "right",
                       height = "450px",
                       tabPanel("Tree/Date", wellPanel(
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
                                            selected = "FSYA")), 
                         wellPanel(dateRangeInput("date", 
                                                  label = "Choose a date range within 2016", 
                                                  start = "2016-03-01", end = "2016-09-30", 
                                                  min = "2016-03-01", max = "2016-09-30" ))), 
                       tabPanel("Phenology/Models/Anot", wellPanel(checkboxGroupInput("pheno",
                                                                                      label = "Show derived phenology?", 
                                                                                      choices = list("Weibull" = "weibull_pheno", "Gompertz" = "gompertz_pheno", "Raw" = "raw_pheno"), 
                                                                                      selected = "Weibull")),
                                wellPanel(checkboxGroupInput("models", "Select a model fit:", 
                                                             choices = list("Weibull" = "weibull", "Gompertz" = "gompertz"), selected = NULL)),
                                wellPanel(radioButtons("image_anot",
                                                       label = "Show annotated image?", 
                                                       choices = list("Yes" = TRUE, "No" = FALSE), selected = FALSE)))), 
                box(width = 4, title = "Info", solidHeader = T,
                    status = "warning", 
                    p("This app helps visualize Dendrometer Data and Microcoring from the same tree/time period. On the middle panel you can control the options to display
                      (by choosing between the tabs) and below you can see the Dendrometer graph. If you click on the points on the graph,
                      you can see a thin-section of the wood on that sampling date and compare it to derived phenology."),
                    br(), 
                    div(h3(strong("Click on the dots to see the corresponding thin-section!"), style = "color:blue")))
                )
              ), 
      tabItem("pheno", 
              fluidRow(
                box(width = 2, 
                    title = "Derived Phenology", 
                    solidHeader = T, 
                    status = "primary",
                    wellPanel(radioButtons("pheno_comparison", 
                                           label = "Choose Phenology\nestimation method", 
                                           choices = list("Weibull" = "weibull", "Gompertz" = "gompertz", "Raw" = "raw"), 
                                           selected = "weibull")
                    )
                ), 
                box(width = 8, 
                    height = 500,
                    title = "Compared Phenology", 
                    solidHeader = T, 
                    status = "primary", 
                    plotOutput(outputId = "derived_phenology"))
              )
      ), 
      tabItem("upload", 
              fluidRow(
                box(width = 4, 
                    title = "Upload User Data",
                    solidHeader = T, 
                    status = "primary", 
                    wellPanel(fileInput("dendro_csv", 
                                        label = "Upload Dendrometer Data", 
                                        multiple = F,
                                        buttonLabel = "Select CSV file"
                                        )
                              )
                    ), 
                box(width = 4, 
                    title = "Upload Xylogenesis Images",
                    solidHeader = T,
                    status = "primary", 
                    wellPanel(fileInput("xylo_images", 
                                        label = "Upload Images", 
                                        multiple = T, 
                                        buttonLabel = "Select Images (.jpeg)"
                      
                    ))
                  
                ), 
                box(width = 4, 
                    title = "Upload Sampling Table",
                    solidHeader = T,
                    status = "primary", 
                    wellPanel(fileInput("sampling_csv", 
                                        label = "Upload CSV", 
                                        multiple = F, 
                                        buttonLabel = "Select CSV"
                                        
                    ))
                    
                )
                      )
              )
      )
    )
  )



# Define server logic
server <- function(input, output) {
  
  ### Re-location of code to wrangle raw data
  
  ### Load data
  
  dendrom_curves <- read_csv("dendro_curves.csv")
  
  dendrom_curves <- dendrom_curves %>% mutate(step_locator = seq_along(tiempo))
  
  dendrom_curves <- dendrom_curves %>% select(tiempo, Dendrometer, tree, step_locator) %>% na.omit() %>% 
    group_by(tree) %>% mutate(min_max_norm = min_max_norm(Dendrometer), 
                              sri = Dendrometer - min(Dendrometer)) %>% ungroup() 
  
  ### Bind Image Data and Dendrometer Data
  
  image_dates <- read_csv("Microcore_sampling_dates.csv")
  
  image_dates <- image_dates %>% gather(key = "date", value = "image",  -Baum, -ID_RCG) %>% unite("tree", Baum, ID_RCG, sep = "")
  
  image_dates$date <- as.POSIXct(image_dates$date, format = "%d.%m.%y", tz = "UTC")
  
  ### Bind with DendrometerData
  
  ### Create the day_wimage object
  
  days_wimage <- image_dates %>% filter(!is.na(image)) 
  
  rm(image_dates)
  
  ### Add a column to the dataframe with the path or name of the image file for the date!
  ###
  
  ### Load image files, or names at least 
  
  bilder <- list.files("./images")
  
  ### you create a list of the files
  ### then you compare which trees match (substr) and which dates (complex posixct and substr)
  ### wrap around which and you have the rows where the files are supposed to go!
  
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
  
  dendrom_curves[row_with_date, "x_images"] <- as.POSIXct(point_dates, format = "%Y-%m-%d", tz = "UTC")### save it as character
  
  dendrom_curves$x_images <- as.POSIXct(dendrom_curves$x_images, origin = "1970-01-01")
  
  ### Add a clear Y value for points with image
  
  y_values <- dendrom_curves %>% filter(file != "NA") %>% group_by(tree, file) %>% summarize(y_value_points_raw  = mean(Dendrometer, na.rm = T), 
                                                                                             y_value_points_norm  = mean(min_max_norm, na.rm = T), 
                                                                                             y_value_points_sri  = mean(sri, na.rm = T))  
  
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
  ### nls_rcg and gompertz_formula are coded in "auxiliary_functions.R" file
  
  dendrom_curves_models_g <- dendrom_curves_models %>% group_by(tree.x) %>% nest() %>%
    mutate(gompertz = map(data, ~ nls_rcg(.)), tidied_gompertz = map(gompertz, tidy), gompertz_fitted  = map(gompertz, augment)) %>%
    unnest(gompertz_fitted)
  
  ### Weibull
  
  ### Getting Weibull initials
  
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
  
  ### Putting both model and raw measurements together
  
  dendrom_curves_models_fertig <- bind_cols(dendrom_curves_models_w, dendrom_curves_models_g, dendrom_curves_models[, c("Dendrometer", "sri", "step_locator")]) %>%
    select(tree.x, min_max_norm, steps, "weibull_fit" = .fitted, 
           "weibull_res" = .resid, "gompertz_fit"= .fitted1, 
           "gompertz_res" = .resid1, "raw_data" = Dendrometer, 
           "SRI" = sri, step_locator)
  
  pos_steps <- match(dendrom_curves_models_fertig$step_locator,dendrom_curves_models$step_locator) ## getting back the tiempo column, get matching positions with orignal table
  
  dendrom_curves_models_fertig[,"tiempo"] <- as.vector(dendrom_curves_models[pos_steps, "tiempo"]) ### put them back carefully (tibbles are messy)
  
  ### Joining Photograph info and data for Points in graph
  
  dendrom_curves_models_fertig <- left_join(dendrom_curves_models_fertig, dendrom_curves[, c("step_locator", "date", "file", "x_images", "y_value_points_norm", "y_value_points_sri")], by = "step_locator")
  
  rm(nls_rcg, dendrom_curves, pos_steps, gompertz_initials, min_max_norm, dendrom_curves_models, dendrom_curves_models_g, dendrom_curves_models_w, 
     initials_weibull, errors, gompertz_formula, i, pos, replace_par, weibull_formula)
  
  ###
  ### Getting Growth phenology
  
  
  growth_phen_gompertz <- growth_phenology_calc(dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = gompertz_fit)
  
  growth_phen_weibull <- growth_phenology_calc(dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = weibull_fit)
  
  growth_phen_change <- growth_phenology_calc(dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = SRI)
  
  
  
  rm(growth_phenology_calc)
  
  ### Server own functions
  ### 
  
  output$dendro_curve <- renderPlot({
    
    tree_of_choice <- input$tree
    
    date_begin <-input$date[1]
    date_end <- input$date[2]
    
    
    
    plot_plain <- dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
      ggplot(aes(color = tree.x, group = tree.x)) + theme(panel.background = element_rect(fill = "white", color = "black"), 
                                                          axis.text.x = element_text(angle = 90, hjust = 1), 
                                                          text = element_text(size = 20)) + labs(color = "Tree")+
      geom_line(aes(x = tiempo, y = min_max_norm)) + geom_point(aes(x = x_images, y = y_value_points_norm), size = 4, alpha = 0.75) + xlab("Time") + 
      ylab("Min-Max Normalized\nStem Radial Increment") + scale_x_datetime(date_breaks = "2 weeks", date_labels = "%b/%d")
    
    
    weibull_pheno_1 <- geom_vline(data = growth_phen_weibull %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_onset, color = tree.x), 
                                  linetype = 1) 
    weibull_pheno_2 <- geom_vline(data = growth_phen_weibull %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_cessation, color = tree.x), 
                                  linetype = 1)  
    
    gompertz_pheno_1 <- geom_vline(data = growth_phen_gompertz %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_onset, color = tree.x), 
                                   linetype = 2) 
    gompertz_pheno_2 <- geom_vline(data = growth_phen_gompertz %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_cessation, color = tree.x), 
                                   linetype = 2) 
    
    raw_pheno_1 <- geom_vline(data = growth_phen_change %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_onset, color = tree.x), 
                              linetype = 3) 
    raw_pheno_2 <- geom_vline(data = growth_phen_change %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_cessation, color = tree.x), 
                              linetype = 3) 
    
    pheno_layers <- list("weibull_pheno" = c(weibull_pheno_1, weibull_pheno_2), "gompertz_pheno" = c(gompertz_pheno_1, gompertz_pheno_2), "raw_pheno" = c(raw_pheno_1, raw_pheno_2))
    
    
    if(is.null(input$pheno)){
      plot_plain
    } else {
      plot_plain <- plot_plain + pheno_layers[input$pheno]
      plot_plain 
    }
    
    
    if(is.null(input$models)){
      plot_plain
    } else {
      weibull <- geom_line(aes(x = tiempo, y = weibull_fit, group = tree.x), color = "black", alpha = 0.7, data = dendrom_curves_models_fertig %>%
                             filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
      gompertz <- geom_line(aes(x = tiempo, y = gompertz_fit, group = tree.x), color = "red", alpha = 0.7, data = dendrom_curves_models_fertig %>%
                              filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
      models <- c("weibull" = weibull, "gompertz" = gompertz)
      
      plot_plain + models[input$models]
    }
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
           contentType = "image/png",
           width = 350,
           height = 300,
           alt = "Image file")
    } else {
      filename <- normalizePath(file.path( './Annot_Images',
                                           as.character(unique(na.omit(nearPoints(dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
                                                                                  input$plot_click, xvar = "x_images", yvar = "y_value_points_norm", 
                                                                                  threshold = 3000000, maxpoints = 2, addDist = T)[,13])[1])[1])))
      filename <- na.omit(filename)[[1]]
      # nearPoints() also works with hover and dblclick events
      list(src = filename,
           contentType = "image/png",
           width = 350,
           height = 300,
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
  
  output$derived_phenology <- renderPlot({
    
    derived_pheno <- input$pheno_comparison 
    
    phenos <- list("gompertz" = growth_phen_gompertz, "weibull" = growth_phen_weibull, "raw" = growth_phen_change)
    
    phenos[[derived_pheno]] %>% select(tree.x, doy_begin, doy_cessation, growth_duration) %>% 
      mutate(species = as.factor(substr(tree.x, 1,3)))%>% gather(key = "parameter", value = "DOY", -tree.x, -species) %>% 
      ggplot(aes(x = species, y = DOY, group = parameter)) + geom_point(aes(color = species)) + facet_wrap(~ parameter) + 
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            text = element_text(size = 20), strip.background = element_rect(color = "black", fill = "white")) + 
      labs(color = "Species")
    
    
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
