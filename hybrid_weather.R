library(ncdf4)
library(ncdf4.helpers)


extract_variable_from_year <- function(nc, variable, year) {
  origin <- ncatt_get(nc, "time")$time_origin
  # var <- ncvar_get(nc, variable)

  dates <- as.Date(nc$dim$time$vals, origin = origin)
  oldest_date <- min(dates[as.integer(format(dates, "%Y")) == year])
  recent_date <- max(dates[as.integer(format(dates, "%Y")) == year])
  start_idx_time <- as.integer(oldest_date - as.Date(origin)) + 1
  end_idx_time <- as.integer(recent_date - as.Date(origin)) + 1

  # var_year <- var[,,start_idx_time:end_idx_time]
  var_year <- nc.get.var.subset.by.axes(nc, variable, list(T = start_idx_time:end_idx_time))
  
  return(var_year)
}

combine_year_var <- function(base_var, fill_var, cutoff_doy) {
  base <- base_var[,,1:cutoff_doy]
  fill <- fill_var[,,(cutoff_doy + 1):dim(fill_var)[3]]
  
  hybrid <- abind::abind(base, fill, along = 3)
  
  return(hybrid)
}

make_hybrid_weather <- function(nc, years, cutoff_date) {
  cutoff_doy <- as.integer(format(cutoff_date, "%j"))
  base_year <- as.integer(format(as.Date(cutoff_date), "%Y"))
  
  for (variable in names(nc$var)) {
    base_var <- extract_variable_from_year(nc, variable, base_year)
    
    for (year in years) {
      fill_var <- extract_variable_from_year(nc, variable, year)
      hybrid <- combine_year_var(base_var, fill_var, cutoff_doy)
      
      origin <- ncatt_get(nc, "time")$time_origin
      start_idx_time <- as.integer(as.Date(paste0(year, "-01-01")) - as.Date(origin)) + 1
      end_idx_time <- as.integer(as.Date(paste0(year, "-12-31")) - as.Date(origin)) + 1
      
      nc.put.var.subset.by.axes(nc, variable, hybrid, list(T = start_idx_time:end_idx_time))
    }
  }
}

run_hybrid_weather <- function(weather_path) {
  nc <- nc_open(weather_path, write = TRUE)
  oldest_date <- as.Date(ncatt_get(nc, "time")$time_origin)
  recent_date <- oldest_date + max(nc$dim$time$vals)
  oldest_year <- as.integer(format(oldest_date, "%Y"))
  recent_year <- as.integer(format(recent_date, "%Y"))
  years <- oldest_year:(recent_year-1)
  make_hybrid_weather(nc, years, cutoff_date = recent_date)
  nc_close(nc)
}


# weather_path <- "merged_ERA5_2000-01-01–2025-04-17.nc"
# run_hybrid_weather(weather_path)



