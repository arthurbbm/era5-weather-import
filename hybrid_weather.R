library(ncdf4)
library(ncdf4.helpers)
library(lubridate)


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

map_fill_to_base <- function(doy_fill, n_fill, n_base) {
  if (n_fill == n_base) return(doy_fill)
  # fill is leap, base not → remove Feb 29 (DOY 60)
  if (n_fill == 366 && n_base == 365 && doy_fill >= 60) return(doy_fill - 1)
  # fill is 365, base is leap → insert Feb 29 after DOY 59
  if (n_fill == 365 && n_base == 366 && doy_fill >= 59) return(doy_fill + 1)
  return(doy_fill)
}

combine_year_var <- function(base_var, fill_var, start_fill_doy, end_fill_doy) {
  n_base <- dim(base_var)[3]
  n_fill <- dim(fill_var)[3]
  
  # start_base <- map_fill_to_base(start_fill_doy, n_fill, n_base)
  # end_base   <- map_fill_to_base(end_fill_doy,   n_fill, n_base)
  # slice three pieces:
  # seg1 <- if (start_base >  1)            base_var[,,      1:(start_base-1)] else NULL
  # seg2 <-                                   fill_var[,, start_fill_doy:end_fill_doy]
  # seg3 <- if (end_base   < n_base)       base_var[,,    (end_base+1):n_base]   else NULL
  
  start_adj_fill <- map_fill_to_base(start_fill_doy, n_fill, n_base)
  end_adj_fill   <- map_fill_to_base(end_fill_doy,   n_fill, n_base)
  
  seg1 <- if (start_adj_fill >  1)            base_var[,,      1:(start_fill_doy-1)] else NULL
  seg2 <-                                   fill_var[,, start_fill_doy:end_fill_doy]
  seg3 <- if (end_adj_fill   < n_base)       base_var[,,    (end_fill_doy+1):n_base]   else NULL
  
  # bind along time
  parts <- Filter(Negate(is.null), list(seg1, seg2, seg3))
  hybrid <- do.call(abind::abind, c(parts, along = 3))
  
  return(hybrid)
}


make_hybrid_weather <- function(nc, years, start_cutoff_date, end_cutoff_date) {
  start_fill_doy <- as.integer(format(start_cutoff_date, "%j"))
  end_fill_doy   <- as.integer(format(end_cutoff_date,   "%j"))
  base_year      <- as.integer(format(start_cutoff_date, "%Y"))
  
  for (variable in names(nc$var)) {
    base_var <- extract_variable_from_year(nc, variable, base_year)
    
    for (year in years) {
      fill_var <- extract_variable_from_year(nc, variable, year)
      hybrid  <- combine_year_var(base_var, fill_var,
                                  start_fill_doy, end_fill_doy)
      
      origin     <- as.Date(ncatt_get(nc, "time")$time_origin)
      year_start <- as.Date(sprintf("%d-01-01", year))
      year_end   <- as.Date(sprintf("%d-12-31", year))
      start_idx  <- as.integer(year_start - origin) + 1
      end_idx    <- as.integer(year_end - origin) + 1
      
      nc.put.var.subset.by.axes(nc, variable, hybrid, 
                                list(T = start_idx:end_idx))
    }
  }
}

# run_hybrid_weather <- function(weather_path) {
#   nc <- nc_open(weather_path, write = TRUE)
#   origin      <- as.Date(ncatt_get(nc, "time")$time_origin)
#   max_offset  <- max(nc$dim$time$vals)
#   recent_date <- origin + max_offset
#   oldest_year <- as.integer(format(origin, "%Y"))
#   recent_year <- as.integer(format(recent_date, "%Y"))
#   years       <- seq(oldest_year, recent_year - 1)
#   start_fill_date <- recent_date + 1
#   base_year <- recent_year
#   end_fill_date <- as.Date(sprintf("%d-12-31", base_year))
#   
#   if (start_fill_date <= end_fill_date) 
#     make_hybrid_weather(nc, years, start_fill_date, end_fill_date)
#   
#   nc_close(nc)
# }

replace_base_doy_from_fill <- function(base_path, fill_path, doy_range, fill_year) {
  base_nc <- nc_open(base_path, write = TRUE)
  fill_nc <- nc_open(fill_path)
  on.exit({ nc_close(base_nc); nc_close(fill_nc) })
  
  origin_base <- as.Date(ncatt_get(base_nc, "time")$time_origin)
  base_dates <- as.Date(base_nc$dim$time$vals, origin = origin_base)
  base_years <- sort(unique(as.integer(format(base_dates, "%Y"))))
  
  start_doy <- doy_range[1]
  end_doy <- doy_range[2]
  
  for (v in names(base_nc$var)) {
    fill_var <- extract_variable_from_year(fill_nc, v, fill_year)
    
    for (yr in base_years[!base_years %in% fill_year]) {
      base_var <- extract_variable_from_year(base_nc, v, yr)
      hybrid <- combine_year_var(base_var, fill_var, start_doy, end_doy)
      y0 <- as.Date(sprintf("%d-01-01", yr))
      y1 <- y0 + dim(hybrid)[[3]] - 1
      start_idx <- as.integer(y0 - origin_base) + 1
      end_idx <- as.integer(y1 - origin_base) + 1
      
      nc.put.var.subset.by.axes(base_nc, v, hybrid,list(T = start_idx:end_idx))
    }
  }
  
  invisible(TRUE)
}

run_hybrid_weather <- function(weather_path, fill_year = lubridate::year(Sys.Date())) {
  nc0 <- nc_open(weather_path)
  origin <- as.Date(ncatt_get(nc0, "time")$time_origin)
  max_offset <- max(nc0$dim$time$vals)
  nc_close(nc0)
  recent_date <- origin + max_offset
  oldest_year <- as.integer(format(origin, "%Y"))
  recent_year <- as.integer(format(recent_date, "%Y"))
  # fill_years <- seq(oldest_year, recent_year - 1)
  start_fill_date <- as.Date(sprintf("%d-01-01", recent_year))
  end_fill_date <- recent_date
  # start_fill_date <- recent_date + 1
  # end_fill_date <- as.Date(sprintf("%d-12-31", recent_year))
  doy_range <- c(
    as.integer(format(start_fill_date, "%j")),
    as.integer(format(end_fill_date, "%j"))
  )

  replace_base_doy_from_fill(weather_path, weather_path, doy_range, fill_year)
}

run_hybrid_with_forecast_weather <- function(weather_path, forecast_path) {
  fc_nc <- nc_open(forecast_path)
  on.exit(nc_close(fc_nc))
  
  fc_origin <- as.Date(ncatt_get(fc_nc, "time")$time_origin)
  fc_times  <- fc_nc$dim$time$vals
  fc_dates  <- fc_origin + fc_times
  
  var_names <- names(fc_nc$var)
  
  complete_times <- sapply(seq_along(fc_times), function(t_idx) {
    all(sapply(var_names, function(vn) {
      vals <- ncvar_get(fc_nc, vn,
                        start = c(1, 1, t_idx),
                        count = c(-1, -1, 1))
      !any(is.na(vals))
    }))
  })
  
  if (!any(complete_times)) {
    stop("No forecast time step has all cells non-NA.")
  }
  
  start_fc <- fc_dates[min(which(complete_times))]
  end_fc   <- fc_dates[max(which(complete_times))]
  
  fill_year <- unique(as.integer(format(fc_dates, "%Y")))
  if (length(fill_year) != 1) {
    stop("Forecast spans multiple years: ", paste(fill_year, collapse = ", "))
  }
  doy_range <- as.integer(c(format(start_fc, "%j"), format(end_fc, "%j")))
  
  replace_base_doy_from_fill(weather_path, forecast_path, doy_range, fill_year)
}


# run_hybrid_weather('/Users/arthur/Documents/SoybeanFarmingSystems/era5-weather-import/merged_ERA5_2000-01-01–2025-04-17 copy.nc')

# run_hybrid_with_forecast_weather(weather_path = '/Users/arthur/Documents/SoybeanFarmingSystems/era5-weather-import/merged_ERA5_2000-01-01–2025-04-17 copy.nc',
#                                  forecast_path = '/Users/arthur/Documents/SoybeanFarmingSystems/era5-weather-import/aifs_20250621_03h_daily copy.nc')
