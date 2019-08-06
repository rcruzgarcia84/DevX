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
    sidebarMenu(id = "menu",
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
                box(width = 12, title = "Instructions", solidHeader = T,
                    status = "primary", 
                    h3(strong("Instructions")),
                    p("To upload user data and combine dendrometer measurements and xylogenesis records 
                      you will need the following items:"),
                    p(strong("- Dendrometer records"), " in a .csv file with three columns ('tree', 'dendrometer' and 'time', all lower case)"),
                    img(src = "dendro_csv.png", align = "center", width = "35%", heigth = "35%"),
                    br(),
                    p(strong("- Xylogenesis images"), " in .png format (Tree ID and date of sampling have to be in the file name, as in 'TREEID_DDMMYYYY.png'.
                      On the upload box for Xylogenesis images, you need to indicate DevX on which character the TREEID ends
                      (4 in the example image below) and 
                      from which to which character can it find the date 
                      (6 to 13 in the example) - TREEID has to match the IDs in the .csv file!! "),
                    img(src = "xylo_files.png", align = "center", width = "35%", heigth = "35%"),
                                        br(), 
                    div(h3(strong("Click on the button below to start the visualization!"), style = "color:blue"), 
                        br(),
                        h4("If you uploaded data, click on the left button after finishing uploading"), 
                        br(),
                        h4("If you want to use pre-loaded data from northeastern Germany, click the right button!")
                        
                        )
                    )
              ),
              fluidRow(  box(width = 4, 
                             title = "Load User Data", 
                             solidHeader = T, 
                             status = "primary", 
                             wellPanel(actionButton(inputId = "user_data",  
                                                    label = "Use User Uploaded Data?", 
                                                    icon = icon("hat-wizard")
                             ), 
                             actionButton(inputId = "pre_loaded", 
                                          label = "Use pre-loaded Data?")
                             
                             )),
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
                    title = "Upload Annotated Xylogenesis Images",
                    solidHeader = T,
                    status = "primary", 
                    wellPanel(fileInput("xylo_annot_images", 
                                        label = "Upload Images", 
                                        multiple = T, 
                                        buttonLabel = "Select Images (.jpeg)"
                                        
                    ))
                    
                )
              ), 
              fluidRow(box(width = 4, 
                           title = "Upload Xylogenesis Images",
                           solidHeader = T,
                           status = "primary", 
                           wellPanel(fileInput("xylo_images", 
                                               label = "Upload Images", 
                                               multiple = T, 
                                               buttonLabel = "Select Images (.png)"
                                               
                           ), 
                           numericInput(inputId = "tree_id_pos", 
                                        value = 1,
                                        label = "Tree ID Characters End",
                                        min = 1,width = "40%"
                                        
                           ), 
                           numericInput(inputId = "date_start_pos", 
                                        value = 2,
                                        label = "Position of Date within file name",
                                        min = 2,width = "40%"
                                        
                           ), 
                           numericInput(inputId = "date_end_pos", 
                                        value = 3,
                                        label = "Position of last Date character within file name",
                                        min = 2,width = "40%"
                                        
                           )
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
                        tabPanel("Tree/Date", 
                                 wellPanel(
                          h5("Select a tree"), 
                          # checkboxGroupInput(inputId = "tree", 
                          #                    label = NULL,
                          #                    choices = names_trees, 
                          #                    )
                          uiOutput("tree")
                        ), 
                        wellPanel(
                          uiOutput("date")
                          )
                        ),
                        # wellPanel(dateRangeInput("date", 
                        #                          label = "Choose a date range within 2016", 
                        #                          start = "2016-03-01", end = "2016-09-30", 
                        #                          min = "2016-03-01", max = "2016-09-30" ))), 
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
server <- function(input, output, session) {
  
 ### Create reactive object to allow to choose between pre-loaded data and user uploaded Data!
  
  v <- reactiveValues(dendrom_curves_models_fertig = NULL, 
                      names_trees = NULL,
                      growth_phen_gompertz = NULL, 
                      growth_phen_weibull = NULL, 
                      growth_phen_change = NULL, 
                      dendrom_curves = NULL)
  
  ### In case user selects pre_loaded button, do pre-loaded data vis.
  observeEvent(input$pre_loaded, {
    
    v$dendrom_curves <- read_csv("dendro_curves.csv")
    
    v$names_trees <- as.list(unique(v$dendrom_curves$tree)) %>% setNames(., unique(v$dendrom_curves$tree))
    
    #image_dates <- read_csv("Microcore_sampling_dates.csv")
    
    ### Load image files, or names at least
    bilder <- list.files("./images")
    
    v$dendrom_curves_models_fertig <- wrangling_dev_x(v$dendrom_curves, bilder)  ### Wrangling converted to function and pushed to "auxiliary functions.R"
    
    ### Getting Growth phenology
    
    v$growth_phen_gompertz <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = gompertz_fit)
    
    v$growth_phen_weibull <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = weibull_fit)
    
    v$growth_phen_change <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = SRI)
    
    
    updateTabItems(session, inputId = "menu", selected = "dendro")
    
  })
  
  ### In case user uploads data, use his! (if there is no uploaded data it crashes!)
  
  observeEvent(input$user_data, {
    
    v$dendrom_curves <- read_csv(input$dendro_csv$datapath)
    
    v$names_trees <- as.list(unique(v$dendrom_curves$tree)) %>% setNames(., unique(v$dendrom_curves$tree))
    
    #image_dates <- read_csv(input$sampling_csv$datapath)
    
    ### Load image files, or names at least
    
    
    bilder <- input$xylo_images[,"name"] ### "xylo_images"
    file_path <- input$xylo_images[,"datapath"]
    annot_file_path <- input$xylo_annot_images[, "datapath"]
    
    
    
    v$dendrom_curves_models_fertig <- wrangling_dev_x(v$dendrom_curves, bilder = input$xylo_images[,"name"],
                                                      user_data = T, file_path = input$xylo_images[,"datapath"], 
                                                      annot_file_path = input$xylo_annot_images[, "datapath"], 
                                                      tree_id_pos_last = input$tree_id_pos, date_start_pos = input$date_start_pos, 
                                                      date_end_pos = input$date_end_pos)
    
    ### Getting Growth phenology
    
    v$growth_phen_gompertz <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = gompertz_fit)
    
    v$growth_phen_weibull <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = weibull_fit)
    
    v$growth_phen_change <- growth_phenology_calc(v$dendrom_curves_models_fertig, group_by_var = tree.x, fitted_col = SRI)
    
    updateTabItems(session, inputId = "menu", selected = "dendro")
  })
  
  
  

  

  ### Server own functions
  ### 
  
  ### Make it reactively show on UI
  
  output$tree <- renderUI({checkboxGroupInput(inputId = "tree", 
                                              label = NULL,
                                              choices = v$names_trees, 
                                              selected = v$names_trees[1]
  )})
  
  output$date <- renderUI({
    dateRangeInput(inputId = "date", 
  label = "Choose a date range", 
   start=  min(as.Date(v$dendrom_curves$time)), end = max(as.Date(v$dendrom_curves$time)), 
  min = min(as.Date(v$dendrom_curves$time)), max = max(as.Date(v$dendrom_curves$time)))
})
  
  output$dendro_curve <- renderPlot({
    

    tree_of_choice <- input$tree
    
    date_begin <-input$date[1]
    date_end <- input$date[2]
    
    
    plot_plain <- v$dendrom_curves_models_fertig %>% filter(tree.x %in% tree_of_choice, time >= date_begin & time <= date_end) %>%
      ggplot(aes(color = tree.x, group = tree.x)) + theme(panel.background = element_rect(fill = "white", color = "black"), 
                                                          axis.text.x = element_text(angle = 90, hjust = 1), 
                                                          text = element_text(size = 20)) + labs(color = "Tree")+
      geom_line(aes(x = time, y = min_max_norm)) + geom_point(aes(x = x_images, y = y_value_points_norm), size = 4, alpha = 0.75) + xlab("Time") + 
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
      weibull <- geom_line(aes(x = time, y = weibull_fit, group = tree.x), color = "black", alpha = 0.7, data = v$dendrom_curves_models_fertig %>%
                             filter(tree.x %in% tree_of_choice, time >= date_begin & time <= date_end)) 
      gompertz <- geom_line(aes(x = time, y = gompertz_fit, group = tree.x), color = "red", alpha = 0.7, data = v$dendrom_curves_models_fertig %>%
                              filter(tree.x %in% tree_of_choice, time >= date_begin & time <= date_end)) 
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
                                                                                 threshold = 3000000, maxpoints = 2, addDist = T)[,"file"])[1])[1])))
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
                                                                                  threshold = 3000000, maxpoints = 2, addDist = T)[,"file"])[1])[1])))
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
                                                                                    threshold = 3000000, maxpoints = 2, addDist = T)[,"data_path"])[1])[1])))
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
                                                                                     threshold = 3000000, maxpoints = 2, addDist = T)[,"annot_datapath"])[1])[1])))
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
                                                                               threshold = 3000000, maxpoints = 2, addDist = T)[,"file"])[1])[1])))
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
