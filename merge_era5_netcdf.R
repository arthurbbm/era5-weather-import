library(ncdf4)

# Set directories
in_dir <- "./data"
out_dir <- dirname(in_dir)

run_merge_era5_netcdf <- function(in_dir, out_dir) {
  # Find NetCDF files
  nc_files <- list.files(in_dir, pattern = "\\.nc$", full.names = TRUE)
  stopifnot(length(nc_files) > 1)
  
  # Read dimensions from the first file
  nc0 <- nc_open(nc_files[1])
  lon_dim <- ncdim_def("longitude", nc0$dim$longitude$units, nc0$dim$longitude$vals)
  lat_dim <- ncdim_def("latitude", nc0$dim$latitude$units, nc0$dim$latitude$vals)
  time_units <- nc0$dim$valid_time$units
  nc_close(nc0)
  
  # Mapping of Copernicus variable names
  copernicus_dictionary <- list(
    maximum_2m_temperature = "tmax",
    minimum_2m_temperature = "tmin",
    sum_total_precipitation = "rain",
    sum_surface_solar_radiation_downwards = "srad",
    mean_2m_dewpoint_temperature = "dewp",
    mean_2m_temperature = "t2m"
  )
  
  # Extract metadata from each file
  file_info <- lapply(nc_files, function(path) {
    nc <- nc_open(path)
    on.exit(nc_close(nc), add = TRUE)
    
    orig_var <- names(nc$var)[1]
    missval <- nc$var[[orig_var]]$missval
    orig_units <- nc$var[[orig_var]]$units        # <<< grab original units
    
    base_name <- sub("\\..*$", "", basename(path))
    parts <- strsplit(base_name, "_")[[1]]
    key <- paste(parts[-c(1,2)], collapse = "_")
    short <- copernicus_dictionary[[key]]
    longname <- tools::toTitleCase(gsub("_", " ", key))
    
    t_vec <- ncvar_get(nc, "valid_time")
    start <- as.Date(strsplit(ncatt_get(nc, "valid_time")$units, " ")[[1]][3])
    
    list(
      path = path,
      short = short,
      orig_var = orig_var,
      longname = longname,
      units = orig_units,
      times = t_vec,
      start = start,
      missval = missval
    )
  })
  
  conv_map <- list(
    srad = list(fun = function(x) x / 1e6, units = "MJ m-2"),
    rain = list(fun = function(x) x * 1000, units = "mm"),
    tmax = list(fun = function(x) x - 273.15, units = "degree_Celsius"),
    tmin = list(fun = function(x) x - 273.15, units = "degree_Celsius"),
    t2m = list(fun = function(x) x - 273.15, units = "degree_Celsius"),
    dewp = list(fun = function(x) x - 273.15, units = "degree_Celsius")
  )
  
  # Order files chronologically
  file_info <- file_info[order(sapply(file_info, `[[`, "start"))]
  
  # Build global time dimension
  starts <- as.Date(sapply(file_info, `[[`, "start"))
  all_dates <- as.Date(unlist(mapply(function(info, st) seq.Date(from = st,
                                                                  by = "days",
                                                                  length.out = length(info$times)),
                                     file_info, starts, SIMPLIFY = FALSE)))
  origin <- min(starts)
  time_vals <- as.integer(all_dates - origin)
  unique_times <- unique(time_vals)
  
  time_dim <- ncdim_def("time", time_units, unique_times)
  
  # Define variables
  vars <- unique(sapply(file_info, `[[`, "short"))
  var_defs <- lapply(vars, function(name) {
    info <- file_info[[ match(name, sapply(file_info, `[[`, "short")) ]]
    # pick converted unit if in map, otherwise fall back to original
    unitstr <- if (name %in% names(conv_map)) conv_map[[name]]$units
    else info$units
    
    ncvar_def(
      name,
      units = unitstr,
      dim = list(lon_dim, lat_dim, time_dim),
      missval = info$missval,
      longname = info$longname
    )
  })
  
  start_date <- format(min(all_dates), "%Y-%m-%d")
  end_date <- format(max(all_dates), "%Y-%m-%d")
  output_file <- file.path(out_dir, sprintf("merged_ERA5_%s–%s.nc", start_date, end_date))
  out_nc <- nc_create(output_file, var_defs)
  ncatt_put(out_nc, "time", "time_origin", as.character(origin))
  
  # Write data for each variable
  for (name in vars) {
    infos <- Filter(function(x) x$short == name, file_info)
    tcur <- 1
    
    for (info in infos) {
      nc <- nc_open(info$path)
      on.exit(nc_close(nc), add = TRUE)
      
      dat <- ncvar_get(nc, info$orig_var)
      nt <- dim(dat)[3]
      
      # if there's a conversion for this short name, do it:
      if (name %in% names(conv_map)) {
        dat <- conv_map[[name]]$fun(dat)
      }
      
      ncvar_put(
        out_nc, name,
        vals = dat,
        start = c(1, 1, tcur),
        count = c(-1, -1, nt)
      )
      
      tcur <- tcur + nt
      nc_close(nc)
    }
  }
  
  title_text <- sprintf("ERA5 variables %s–%s in one file", start_date, end_date)
  
  ncatt_put(out_nc, 0, "title", title_text)
  nc_close(out_nc)
  message(basename(output_file), " created successfully")
  
  return(output_file)
}


# run_merge_era5_netcdf(in_dir, out_dir)
