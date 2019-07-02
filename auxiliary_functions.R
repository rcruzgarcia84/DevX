#### Secondary functions for Dendro Viewer

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