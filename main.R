library(ecmwfr)
library(future)
library(future.apply)
library(assertthat)
library(stringr)

source("merge_era5_netcdf.R")


month_request_generator <- function(year, month, variable, stat = "mean", 
                                    area = c(40.613687151061505, -95.76804054400543, 
                                             35.99538225441254, -89.09913230512305)) {
  assertthat::assert_that(length(variable) == 1)
  assertthat::assert_that(length(stat) == 1)
  month <- stringr::str_pad(month, width = 2, pad = "0")
  list(
    product_type = "reanalysis",
    dataset_short_name = "derived-era5-single-levels-daily-statistics",
    target = paste0("ERA5_", year, "_", stat, "_", variable, ".nc"),
    data_format = "netcdf",
    variable = variable,
    year = year,
    month = month,
    day = stringr::str_pad(1:31, width = 2, pad = "0"),
    daily_statistic = sapply(stat, function(x) paste0("daily_", x), USE.NAMES = FALSE),
    time_zone = "utc-06:00",
    frequency = "1_hourly",
    area = area
  )
}


args = commandArgs(trailingOnly=T)

if (length(args) == 0) stop("No arguments provided. Please provide the json file with the start month, 
                            end month, start year, end year, and bounding box.")
file.exists(args[1]) || stop("The provided file does not exist. Please provide a valid json file.")
json_file <- args[[1]]
json_data <- jsonlite::read_json(json_file, simplifyDataFrame = TRUE)

uid <- json_data$uid
key <- json_data$key
start_month <- json_data$start_month
end_month <- json_data$end_month
start_year <- json_data$start_year
end_year <- json_data$end_year
bounding_box <- json_data$bounding_box

month <- start_month:end_month
years <- start_year:end_year
bbox <- c(bounding_box$north, bounding_box$west,
          bounding_box$south, bounding_box$east)

wf_set_key(key = key, user = uid)

# Define the download path variable
download_path <- "./data"

# Create a list of requests for each variable for each year
requests <- list()
for (yr in years) {
  requests[[paste0("mean_dewpoint_temperature_", yr)]] <- month_request_generator(
    month = month, 
    year = yr, 
    variable = "2m_dewpoint_temperature", 
    stat = "mean",
    area = bbox
  )
  requests[[paste0("mean_temperature_", yr)]] <- month_request_generator(
    month = month, 
    year = yr, 
    variable = "2m_temperature", 
    stat = "mean",
    area = bbox
  )
  requests[[paste0("minimum_temperature_", yr)]] <- month_request_generator(
    month = month, 
    year = yr, 
    variable = "2m_temperature", 
    stat = "minimum",
    area = bbox
  )
  requests[[paste0("maximum_temperature_", yr)]] <- month_request_generator(
    month = month, 
    year = yr, 
    variable = "2m_temperature", 
    stat = "maximum",
    area = bbox
  )
  requests[[paste0("sum_precipitation_", yr)]] <- month_request_generator(
    month = month, 
    year = yr, 
    variable = "total_precipitation", 
    stat = "sum",
    area = bbox
  )
  requests[[paste0("sum_radiation_", yr)]] <- month_request_generator(
    month = month,
    year = yr,
    variable = "surface_solar_radiation_downwards",
    stat = "sum",
    area = bbox
  )
}


requests <- Filter(function(req) {
  target_file <- file.path(download_path, req$target)
  
  if (file.exists(target_file) &&
      req$year != as.integer(format(Sys.Date(), "%Y"))) {
    message("Skipping request for: ", req$target, " (not current year and already downloaded)")
    return(FALSE)  }
  TRUE
}, requests)

Set up futures to run downloads asynchronously.
plan(multisession, workers = future::availableCores() - 1)

# Helper function that submits a single request and polls until the file is available.
download_async <- function(req, uid, download_path) {
  message("Submitting request for: ", req$target)
  req_info <- ecmwfr::wf_request(request = req, user = uid, transfer = FALSE, path = download_path)
  downloaded <- NULL
  # Poll every 5 seconds until the file is available
  while (is.null(downloaded)) {
    Sys.sleep(5)
    downloaded <- tryCatch({
      ecmwfr::wf_transfer(url = req_info$get_url(), user = uid, path = download_path, filename = req$target)
      ecmwfr::wf_delete(url = req_info$get_url(), user = uid)
    }, error = function(e) {
      message("Data for ", req$target, " not yet available. Waiting...")
      NULL
    })
  }
  # Return both the target and file name
  list(target = req$target, file = downloaded)
}

# Launch the downloads asynchronously for each request.
results <- future_lapply(requests, download_async, uid = uid, download_path = download_path)

# Map the downloaded file paths to their target names.
downloaded_files <- setNames(lapply(results, `[[`, "file"), sapply(results, `[[`, "target"))

print(downloaded_files)

# Merge the downloaded NetCDF files into a single NetCDF file
weather_path <- run_merge_era5_netcdf(in_dir = download_path, out_dir = dirname(download_path))

if (json_data$make_hybrid_weather) {
  source("hybrid_weather.R")
  run_hybrid_weather(weather_path)
}
