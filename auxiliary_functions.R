#### Secondary functions for Dendro Viewer

### Wrangling of data starts

wrangling_dev_x <- function(dendrom_curves, image_dates,  bilder){
  
  dendrom_curves <- dendrom_curves %>% mutate(step_locator = seq_along(tiempo)) %>% select(tiempo, Dendrometer, tree, step_locator) %>% na.omit() %>% 
    group_by(tree) %>% mutate(min_max_norm = min_max_norm(Dendrometer), 
                              sri = Dendrometer - min(Dendrometer)) %>% ungroup() 
  
  ### Bind Image Data and Dendrometer Data
  
  image_dates <- image_dates %>% gather(key = "date", value = "image",  -Baum, -ID_RCG) %>% unite("tree", Baum, ID_RCG, sep = "")
  
  image_dates$date <- as.POSIXct(image_dates$date, format = "%d.%m.%y", tz = "UTC")
  
  ### Bind with DendrometerData
  
  ### Create the day_wimage object
  
  days_wimage <- image_dates %>% filter(!is.na(image)) 
  
  ### Add a column to the dataframe with the path or name of the image file for the date!
  ###
  
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
  
  dendrom_curves_models_fertig <- bind_cols(dendrom_curves_models_w, dendrom_curves_models_g,
                                            dendrom_curves_models[, c("Dendrometer", "sri", "step_locator")]) %>%
    select(tree.x, min_max_norm, steps, "weibull_fit" = .fitted, 
           "weibull_res" = .resid, "gompertz_fit"= .fitted1, 
           "gompertz_res" = .resid1, "raw_data" = Dendrometer, 
           "SRI" = sri, step_locator)
  
  pos_steps <- match(dendrom_curves_models_fertig$step_locator,dendrom_curves_models$step_locator) ## getting back the tiempo column, get matching positions with orignal table
  
  dendrom_curves_models_fertig[,"tiempo"] <- as.vector(dendrom_curves_models[pos_steps, "tiempo"]) ### put them back carefully (tibbles are messy)
  
  ### Joining Photograph info and data for Points in graph
  
  dendrom_curves_models_fertig <- left_join(dendrom_curves_models_fertig, dendrom_curves[, c("step_locator", "date", "file", "x_images", "y_value_points_norm", "y_value_points_sri")], by = "step_locator")
  
  dendrom_curves_models_fertig
  
  
}


### Range normalization function

min_max_norm <- function(x){
  (x - min(x))/(max(x) - min(x))
}

### Gompertz Requirements
### getting initials for gompertz

gompertz_initials <- function(x){
  i <- which.max(diff(x$min_max_norm))
  starting.values <- c(A = max(x$min_max_norm), 
                       mu = max(diff(x$min_max_norm)/(x[i+1,"steps"]-x[i, "steps"])), 
                       lambda = i)
  starting.values
}

gompertz_formula <- formula(min_max_norm ~ A * exp(-exp(mu * exp(1) / A * (lambda - steps) + 1)))

nls_rcg <- function(x){
  nls(gompertz_formula, data = as.data.frame(x), 
      start = gompertz_initials(as.data.frame(x)))
}

### WEibull functions


weibull_formula <- formula(min_max_norm ~ Asym-Drop*exp(-exp(lrc)* steps^pwr))

### Growth Phenology



### Easy function for growth phenology (speficfy the "fitted_col" function to get 
### the estimates for agiven model)
growth_phenology_calc <- function(x, group_by_var, fitted_col){
  group_by_var <- enquo(group_by_var)
  fitted_col <- enquo(fitted_col)
  
  x %>% group_by(!!group_by_var) %>%
    summarise(max(!!fitted_col), min(!!fitted_col), total_amplitude = max(!!fitted_col) - min(!!fitted_col), five_percent = total_amplitude * 0.05, ninety_five = total_amplitude *0.95, begin = max(!!fitted_col) - ninety_five, 
              cessation = max(!!fitted_col) - five_percent, date_onset = tiempo[which.min(abs((!!fitted_col) - begin))], date_cessation = tiempo[which.min(abs((!!fitted_col) - cessation))], 
              growth_duration = date_cessation - date_onset, doy_begin = yday(date_onset), doy_cessation = yday(date_cessation)) ##using which.min(abs("values" - "value to query")) you get the position of the closest value!
}
