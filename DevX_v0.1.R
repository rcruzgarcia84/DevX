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
      menuItem("Data Selection", icon = icon("upload"), tabName = "upload"),
      menuItem("Dendro/Thinsection\nViewer", tabName = "dendro", icon = icon("area-chart")), 
      menuItem("Phenology Comparison", tabName = "pheno", icon = icon("envira")), 
      menuItem("Source Code", icon = icon("code"), href = "https://github.com/rcruzgarcia84/DevX/blob/master/DevX_v0.1.R")
    )
  ),
  
  dashboardBody(
    tabItems(
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
              ), 
              fluidRow(
                box(width = 4, 
                    title = "Load User Data", 
                    solidHeader = T, 
                    status = "primary", 
                    wellPanel(actionButton(inputId = "user_data",  
                                           label = "Use User Uploaded Data?", 
                                           icon = icon("hat-wizard")
                    ), 
                    actionButton(inputId = "pre_loaded", 
                                 label = "Use pre-loaded Data?")
                    
                    ))
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
                          # checkboxGroupInput(inputId = "tree", 
                          #                    label = NULL,
                          #                    choices = names_trees, 
                          #                    )
                          uiOutput("tree")
                        ), 
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
      ))
    )
  )



# Define server logic
server <- function(input, output) {
  
 ### Create reactive object to allow to choose between pre-loaded data and user uploaded Data!
  
  v <- reactiveValues(dendrom_curves_models_fertig = NULL, 
                      names_trees = NULL,
                      growth_phen_gompertz = NULL, 
                      growth_phen_weibull = NULL, 
                      growth_phen_change = NULL)
  
  ### In case user selects pre_loaded button, do pre-loaded data vis.
  observeEvent(input$pre_loaded, {
    
    dendrom_curves <- read_csv("dendro_curves.csv")
    
    v$names_trees <- as.list(unique(dendrom_curves$tree)) %>% setNames(., unique(dendrom_curves$tree))
    
    image_dates <- read_csv("Microcore_sampling_dates.csv")
    
    ### Load image files, or names at least
    bilder <- list.files("./images")
    
    v$dendrom_curves_models_fertig <- wrangling_dev_x(dendrom_curves, image_dates, bilder)  ### Wrangling converted to function and pushed to "auxiliary functions.R"
    
    ### Getting Growth phenology
    
    v$growth_phen_gompertz <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = gompertz_fit)
    
    v$growth_phen_weibull <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = weibull_fit)
    
    v$growth_phen_change <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = SRI)
    
  })
  
  ### In case user uploads data, use his! (if there is no uploaded data it crashes!)
  
  observeEvent(input$user_data, {
    
    dendrom_curves <- read_csv(input$dendro_csv$datapath)
    
    v$names_trees <- as.list(unique(dendrom_curves$tree)) %>% setNames(., unique(dendrom_curves$tree))
    
    image_dates <- read_csv(input$sampling_csv$datapath)
    
    ### Load image files, or names at least
    
    
    bilder <- input$xylo_images[,1] ### "xylo_images"
    file_path <- input$xylo_images[,4]
    
    v$dendrom_curves_models_fertig <- wrangling_dev_x(dendrom_curves, image_dates, bilder, user_data = T, file_path = file_path)
    
    ### Getting Growth phenology
    
    v$growth_phen_gompertz <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = gompertz_fit)
    
    v$growth_phen_weibull <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = weibull_fit)
    
    v$growth_phen_change <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = SRI)
    
  })
  
  
  

  

  ### Server own functions
  ### 
  
  ### Make it reactively show on UI
  
  output$tree <- renderUI({checkboxGroupInput(inputId = "tree", 
                                              label = NULL,
                                              choices = v$names_trees, 
                                              selected = v$names_trees[1]
  )})
  
  output$dendro_curve <- renderPlot({
    

    tree_of_choice <- input$tree
    
    date_begin <-input$date[1]
    date_end <- input$date[2]
    
    
    
    plot_plain <- v$dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end) %>%
      ggplot(aes(color = tree.x, group = tree.x)) + theme(panel.background = element_rect(fill = "white", color = "black"), 
                                                          axis.text.x = element_text(angle = 90, hjust = 1), 
                                                          text = element_text(size = 20)) + labs(color = "Tree")+
      geom_line(aes(x = tiempo, y = min_max_norm)) + geom_point(aes(x = x_images, y = y_value_points_norm), size = 4, alpha = 0.75) + xlab("Time") + 
      ylab("Min-Max Normalized\nStem Radial Increment") + scale_x_datetime(date_breaks = "2 weeks", date_labels = "%b/%d")
    
    
    weibull_pheno_1 <- geom_vline(data = v$growth_phen_weibull %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_onset, color = tree.x), 
                                  linetype = 1) 
    weibull_pheno_2 <- geom_vline(data = v$growth_phen_weibull %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_cessation, color = tree.x), 
                                  linetype = 1)  
    
    gompertz_pheno_1 <- geom_vline(data = v$growth_phen_gompertz %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_onset, color = tree.x), 
                                   linetype = 2) 
    gompertz_pheno_2 <- geom_vline(data = v$growth_phen_gompertz %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_cessation, color = tree.x), 
                                   linetype = 2) 
    
    raw_pheno_1 <- geom_vline(data = v$growth_phen_change %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_onset, color = tree.x), 
                              linetype = 3) 
    raw_pheno_2 <- geom_vline(data = v$growth_phen_change %>% filter(tree.x %in% tree_of_choice), aes(xintercept = date_cessation, color = tree.x), 
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
      weibull <- geom_line(aes(x = tiempo, y = weibull_fit, group = tree.x), color = "black", alpha = 0.7, data = v$dendrom_curves_models_fertig %>%
                             filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
      gompertz <- geom_line(aes(x = tiempo, y = gompertz_fit, group = tree.x), color = "red", alpha = 0.7, data = v$dendrom_curves_models_fertig %>%
                              filter(tree.x %in% tree_of_choice, tiempo >= date_begin & tiempo <= date_end)) 
      models <- c("weibull" = weibull, "gompertz" = gompertz)
      
      plot_plain + models[input$models]
    }
  })
  
  
  output$image <- renderImage({
    
    
    ### Separate method for Pre-loaded data and for User uploaded data
    
    if(input$pre_loaded){
    
    tree_of_choice <- input$tree
    # With base graphics, need to tell it what the x and y variables are.
    if(input$image_anot != TRUE){
      filename <- normalizePath(file.path('./images',
                                          as.character(unique(na.omit(nearPoints(v$dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
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
                                           as.character(unique(na.omit(nearPoints(v$dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
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
    } else {
      
       tree_of_choice <- input$tree
       # With base graphics, need to tell it what the x and y variables are.
       if(input$image_anot != TRUE){
         filename <- normalizePath(file.path(as.character(unique(na.omit(nearPoints(v$dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
                                                                                    input$plot_click, xvar = "x_images", yvar = "y_value_points_norm",
                                                                                    threshold = 3000000, maxpoints = 2, addDist = T)[,"datapath"])[1])[1])))
         filename <- na.omit(filename)[[1]]
         # nearPoints() also works with hover and dblclick events
         list(src = filename,
              contentType = "image/png",
              width = 350,
              height = 300,
              alt = "Image file")
       } else {
         filename <- normalizePath(file.path(as.character(unique(na.omit(nearPoints(v$dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
                                                                                     input$plot_click, xvar = "x_images", yvar = "y_value_points_norm",
                                                                                     threshold = 3000000, maxpoints = 2, addDist = T)[,"datapath"])[1])[1])))
         filename <- na.omit(filename)[[1]]
         # nearPoints() also works with hover and dblclick events
         list(src = filename,
              contentType = "image/png",
              width = 350,
              height = 300,
              alt = "Image file")
       }
    }
  }, deleteFile = F)
  
  output$name_image <- renderText({
    tree_of_choice <- input$tree
    # With base graphics, need to tell it what the x and y variables are.
    filename <- normalizePath(file.path('./images',
                                        as.character(unique(na.omit(nearPoints(v$dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice),
                                                                               input$plot_click, xvar = "x_images", yvar = "y_value_points_norm", 
                                                                               threshold = 3000000, maxpoints = 2, addDist = T)[,13])[1])[1])))
    #substr(filename, nchar(filename)-21, nchar(filename))
    filename
  })
  
  output$derived_phenology <- renderPlot({
    
    derived_pheno <- input$pheno_comparison 
    
    phenos <- list("gompertz" = v$growth_phen_gompertz, "weibull" = v$growth_phen_weibull, "raw" = v$growth_phen_change)
    
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
