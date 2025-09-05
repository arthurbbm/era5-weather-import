#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 0) REQUIREMENTS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
library(parallel)
library(gdalUtilities)
library(jsonlite)
library(ncdf4)
library(dplyr)
library(lubridate)


build_forecast_schedule <- function(start_dt,
                                    end_dt   = Sys.time(),
                                    base_url = "https://storage.googleapis.com/ecmwf-open-data",
                                    tz       = "UTC",
                                    temp_dir = "temp") {
  tz            <- if (tz %in% c("CST","CDT")) "America/Chicago" else tz
  step_h        <- hours(3)
  start_dt <- as.POSIXct(start_dt, tz = tz)
  end_dt   <- as.POSIXct(end_dt, tz = tz)
  if (end_dt < start_dt) stop("end_dt must be ≥ start_dt")
  now_run_time <- floor_date(now(tz), "12 hours") - hours(12) # Add delay for current run to assure "oper" is available
  last_run_date <- format(now_run_time, "%Y%m%d")
  last_run_hour <- format(now_run_time, "%H")
  
  tibble(valid_time = seq(start_dt + step_h, end_dt, paste(hour(step_h), "hour"))) %>%
    mutate(
      valid_start     = valid_time - step_h,
      normal_run_time = floor_date(valid_start, "6 hours"),
      normal_run_hour = format(normal_run_time, "%H"),
      normal_stream   = if_else(normal_run_hour %in% c("06","18"), "scda", "oper"),
      is_future       = valid_start > now_run_time,
      run_date        = if_else(is_future, last_run_date, format(normal_run_time, "%Y%m%d")),
      run_hour        = if_else(is_future, last_run_hour, normal_run_hour),
      run_time        = if_else(is_future, now_run_time, normal_run_time),
      stream          = if_else(is_future, "oper", normal_stream),
      lead_h          = as.integer(difftime(valid_time, run_time, units = "hours")),
      url             = sprintf("%s/%s/%sz/ifs/0p25/%s/%s%s0000-%dh-%s-fc.grib2",
                                base_url, run_date, run_hour, stream,
                                run_date, run_hour, lead_h, stream),
      grib_path       = sprintf("./temp/aifs_%s_%sh_%03dh.grib2",
                                run_date, run_hour, lead_h)
    ) %>%
    filter(lead_h <= 144 | lead_h %% 6 == 0) %>%
    select(valid_time, run_date, run_hour, stream, lead_h, url, grib_path)
  
}

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 1) COMPUTE TODAY’S “LAGGED” UTC RUN (one cycle behind)
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
TEMP_DIR = "temp"
FORECAST_DIR = "forecast"
forecast_days_advance <- days(14)

run_hours <- c(0, 6, 12, 18)
utc_now   <- as.POSIXct(Sys.time(), tz="UTC")
today_utc <- as.Date(utc_now)
h_now     <- as.integer(format(utc_now, "%H"))

eligible_idx <- which(run_hours <= h_now)
idx           <- if (length(eligible_idx)==0) 1 else max(eligible_idx)
lag_idx       <- idx - 1
if (lag_idx < 1) {
  lag_idx  <- length(run_hours)
  today_utc <- today_utc - 1
}

run_hour <- sprintf("%02d", run_hours[lag_idx])
run_date <- format(today_utc, "%Y%m%d")
message("Using AIFS run: ", run_date, " ", run_hour, " UTC")

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 2) BUILD URLS, FILENAMES & ADVANCE HOURS VECTOR
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
start_dt <- format(as_date(today_utc) - days(6), "%Y-%m-%d %H:%M:%S")
end_dt <- format(as_date(today_utc) + forecast_days_advance, "%Y-%m-%d %H:%M:%S")

df <- build_forecast_schedule(
  start_dt = start_dt,
  end_dt = end_dt,
  tz = "America/Chicago",
  temp_dir = TEMP_DIR
)

urls <- df$url
grb_paths <- df$grib_path
nc4_paths <- file.path(
  FORECAST_DIR,
  sprintf("aifs_%s_%sh_%03dh.nc",
          df$run_date, df$run_hour, df$lead_h)
)

dir.create(TEMP_DIR)
dir.create(FORECAST_DIR)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 3) (optional) DOWNLOAD ALL GRIB2 FILES IN PARALLEL
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
options(timeout = 1200)
setTimeLimit(elapsed = 1200, transient = TRUE)
mclapply(seq_along(urls), function(i) {
  download.file(urls[i], destfile = grb_paths[i], mode = "wb")
}, mc.cores = detectCores())

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 4) CONVERT EACH TO NETCDF4, IN THE DESIRED ORDER
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
band_map <- list(
  srad = list(params = c(4, 7),
              pos    = c(1, 2)),
  rain = list(params = c(1, 193),
              pos    = c(1, 2)),
  t2m  = list(params = c(0, 0, 2, 3),
              pos    = c(1, 2, 12, 24)),
  dewp = list(params = c(0, 6, 2),
              pos    = c(1, 2, 12))
)
var_order <- names(band_map)   # c("srad","rain","t2m","dewp")

mclapply(seq_along(grb_paths), function(i) {
  # 4.1) read GDAL/GRIB metadata
  json_txt <- system2("gdalinfo",
                      args   = c("-json", grb_paths[i]),
                      stdout = TRUE)
  info     <- fromJSON(paste(json_txt, collapse="\n"))
  md       <- info$bands$metadata[[1]]
  
  # 4.2) pull out the raw PDS strings
  raw_pds <- if (is.data.frame(md)) {
    md$GRIB_PDS_TEMPLATE_ASSEMBLED_VALUES
  } else {
    vapply(md, `[[`, character(1), "GRIB_PDS_TEMPLATE_ASSEMBLED_VALUES")
  }
  
  # 4.3) split & convert into integer vectors
  pds_list <- lapply(strsplit(raw_pds, " "), as.integer)
  
  # 4.4) for each variable, find the band whose PDS vector
  #      matches at the requested positions
  b_idx <- vapply(var_order, function(var) {
    spec <- band_map[[var]]
    idx <- which(vapply(pds_list, function(p) {
      length(p) >= max(spec$pos) &&
        all(p[spec$pos] == spec$params)
    }, logical(1)))
    if (length(idx) != 1) {
      stop("Couldn’t find unique band for ", var)
    }
    idx
  }, integer(1))
  
  # 4.5) translate in that exact sequence
  gdal_translate(
    src_dataset = grb_paths[i],
    dst_dataset = nc4_paths[i],
    of          = "netCDF",
    b           = b_idx,
    co          = c("FORMAT=NETCDF4", "COMPRESS=DEFLATE")
  )
}, mc.cores = detectCores())


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 5) MERGE + RENAME + (only srad) CONVERT
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# conv_map: conversions for non-srad variables
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
conv_map <- list(
  rain = list(fun   = function(x) x * 1e3,          units = "mm"),
  t2m  = list(fun   = NULL,                        units = "degree_Celsius"),
  dewp = list(fun   = NULL,                        units = "degree_Celsius")
)

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# A) Gather filenames, parse run_dt and lead_h
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
fns        <- basename(nc4_paths)
run_dt_vec <- as.POSIXct(
  sub("aifs_(\\d{8}_\\d{2})h_.*", "\\1", fns),
  format = "%Y%m%d_%H", tz = "UTC"
)
lead_h_vec <- as.integer(
  sub(".*_(\\d{3})h\\.nc$", "\\1", fns)
)
valid_time <- run_dt_vec + lead_h_vec * 3600

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# B) Determine overall start_time & adv_h, then sort everything
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
start_time <- min(valid_time)
adv_h      <- as.integer(difftime(valid_time, start_time, units="hours"))

ord        <- order(adv_h)
nc4_paths  <- nc4_paths[ord]
run_dt_vec <- run_dt_vec[ord]
lead_h_vec <- lead_h_vec[ord]
adv_h      <- adv_h[ord]

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# C) Prepare output filename and dimensions
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
run_date <- format(start_time, "%Y%m%d")
run_hour <- format(start_time, "%H")
outname  <- sprintf("aifs_%s_%sh_merged.nc", run_date, run_hour)
if (file.exists(outname)) file.remove(outname)

# inspect first file for lon/lat dims
ncin     <- nc_open(nc4_paths[1])
lon_vals <- ncin$dim$lon$vals; lon_units <- ncin$dim$lon$units
lat_vals <- ncin$dim$lat$vals; lat_units <- ncin$dim$lat$units
old_vars <- names(ncin$var)                 # e.g. c("crs","Band1",…)
nc_close(ncin)

time_units <- sprintf("hours since %s", format(start_time, "%Y-%m-%d %H:%M:%S"))
time_dim   <- ncdim_def("time",  units=time_units, vals=adv_h, unlim=TRUE)
lon_dim    <- ncdim_def("longitude", units=lon_units, vals=lon_vals)
lat_dim    <- ncdim_def("latitude",  units=lat_units, vals=lat_vals)

# define variables
crs_var <- ncvar_def("crs", units="", dim=list(), longname="crs", prec="integer")
new_vars <- c("crs","srad","rain","t2m","dewp")

# build var_defs
var_defs <- list(crs_var)
ncin2    <- nc_open(nc4_paths[1])
for (i in seq_along(old_vars)[-1]) {
  nv <- new_vars[i]
  v0 <- ncin2$var[[ old_vars[i] ]]
  units <- if (nv == "srad") "MJ m-2" else conv_map[[nv]]$units
  var_defs <- c(var_defs, list(
    ncvar_def(
      name    = nv,
      units   = units,
      dim     = list(lon_dim, lat_dim, time_dim),
      missval = v0$missval
    )
  ))
}
nc_close(ncin2)

# create merged file
ncw <- nc_create(outname, vars=var_defs, force_v4=TRUE)
ncvar_put(ncw, "time", adv_h)
ncatt_put(ncw, "time", "calendar",    "gregorian")
ncatt_put(ncw, "time", "axis",        "T")
ncatt_put(ncw, "time", "time_origin", format(start_time, "%Y-%m-%d %H:%M:%S"))

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# D) Loop through each file, compute 3h or 6h increments for rain & srad, write slices
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# track previous raw cumulative for each run_dt-hour & var
prev_raw <- list()

for (ti in seq_along(nc4_paths)) {
  fname <- basename(nc4_paths[ti])
  step_h <- as.integer(sub(".*_(\\d{3})h\\.nc$", "\\1", fname))
  accum_h <- if (step_h <= 144) 3 else 6
  
  # label to separate runs: "YYYYMMDD_HH"
  run_label <- format(run_dt_vec[ti], "%Y%m%d_%H")
  
  ncin3 <- nc_open(nc4_paths[ti])
  for (i in seq_along(old_vars)[-1]) {
    nv        <- new_vars[i]
    raw_slice <- ncvar_get(ncin3, old_vars[i])
    
    if (nv %in% c("rain", "srad")) {
      key <- paste0(run_label, "_", nv)
      if (is.null(prev_raw[[key]])) {
        incr <- raw_slice
      } else {
        incr <- raw_slice - prev_raw[[key]]
      }
      prev_raw[[key]] <- raw_slice
      
      if (nv == "srad") {
        # convert W/m2 over 'accum_h' hours to MJ/m2
        slice <- incr * (accum_h * 3600) / 1e6
        slice <- incr / 1e6
      } else {
        slice <- incr * 1e3
      }
    } else if (!is.null(conv_map[[nv]]$fun)) {
      slice <- conv_map[[nv]]$fun(raw_slice)
    } else {
      slice <- raw_slice
    }
    
    ncvar_put(ncw, nv,
              vals  = slice,
              start = c(1, 1, ti),
              count = c(-1, -1, 1))
  }
  nc_close(ncin3)
}

nc_close(ncw)
message("Merged file written to: ", outname)


#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 6) AGGREGATE HOURLY → DAILY AND WRITE ALL FIELDS
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# A) OPEN MERGED HOURLY FILE
nc       <- nc_open(outname)
longitude      <- nc$dim$longitude$vals;    longitude_units <- nc$dim$longitude$units
latitude      <- nc$dim$latitude$vals;    latitude_units <- nc$dim$latitude$units
time     <- nc$dim$time$vals     # hours since run_origin
rain_h   <- ncvar_get(nc, "rain")
srad_h   <- ncvar_get(nc, "srad")
t2m_h    <- ncvar_get(nc, "t2m")
de_wp_h  <- ncvar_get(nc, "dewp")
nc_close(nc)

# B) BUILD DATE INDEX
run_origin <- as.POSIXct(paste0(run_date, run_hour), format="%Y%m%d%H", tz="UTC")
dates      <- as.Date(run_origin + time * 3600)
# full-year sequence
yrs        <- unique(year(dates))
start_date <- as.Date(paste0(min(yrs), "-01-01"))
end_date   <- as.Date(paste0(max(yrs), "-12-31"))
days       <- seq.Date(start_date, end_date, by = "day")
ndays      <- length(days)
nx <- length(longitude);  ny <- length(latitude)

# C) ALLOCATE DAILY ARRAYS
rain_d <- array(NA, c(nx, ny, ndays))
srad_d <- array(NA, c(nx, ny, ndays))
t2m_d  <- array(NA, c(nx, ny, ndays))
de_wp_d<- array(NA, c(nx, ny, ndays))
tmin_d <- array(NA, c(nx, ny, ndays))
tmax_d <- array(NA, c(nx, ny, ndays))

# D) FUNCTION TO COMPUTE DAILY STATS
compute_day <- function(d) {
  sel <- which(dates == days[d])
  if (length(sel) == 0) {
    return(list(
      rain = matrix(NA, nx, ny),
      srad = matrix(NA, nx, ny),
      t2m  = matrix(NA, nx, ny),
      dewp = matrix(NA, nx, ny),
      tmin = matrix(NA, nx, ny),
      tmax = matrix(NA, nx, ny)
    ))
  }
  rain  <- rowSums(rain_h[,,sel, drop=FALSE], dims=2, na.rm=TRUE)
  srad  <- rowSums(srad_h[,,sel, drop=FALSE], dims=2, na.rm=TRUE)
  t_slice   <- t2m_h[,,sel, drop=FALSE]
  dewp_slice<- de_wp_h[,,sel, drop=FALSE]
  t2m   <- apply(t_slice,   c(1,2), mean, na.rm=TRUE)
  dewp  <- apply(dewp_slice,c(1,2), mean, na.rm=TRUE)
  tmin  <- apply(t_slice,   c(1,2), min,  na.rm=TRUE)
  tmax  <- apply(t_slice,   c(1,2), max,  na.rm=TRUE)
  list(rain = rain, srad = srad, t2m = t2m,
       dewp = dewp, tmin = tmin, tmax = tmax)
}

# E) PARALLEL COMPUTATION
stats_list <- mclapply(seq_along(days), compute_day, mc.cores = detectCores())
for (d in seq_along(days)) {
  st <- stats_list[[d]]
  rain_d[,,d] <- st$rain
  srad_d[,,d] <- st$srad
  t2m_d[,,d]  <- st$t2m
  de_wp_d[,,d]<- st$dewp
  tmin_d[,,d] <- st$tmin
  tmax_d[,,d] <- st$tmax
}

# F) CROP TO BBOX
area    <- c(40.613687151061505, -95.76804054400543,
             35.99538225441254, -89.09913230512305)
latitude_max <- area[1]; longitude_min <- area[2]
latitude_min <- area[3]; longitude_max <- area[4]
longitude_idx <- which(longitude >= longitude_min & longitude <= longitude_max)
latitude_idx <- which(latitude >= latitude_min & latitude <= latitude_max)
longitude <- longitude[longitude_idx]
latitude <- sort(latitude[latitude_idx], decreasing = T)
rain_d  <- rain_d[longitude_idx, latitude_idx, , drop=FALSE]
srad_d  <- srad_d[longitude_idx, latitude_idx, , drop=FALSE]
t2m_d   <- t2m_d[ longitude_idx, latitude_idx, , drop=FALSE]
de_wp_d <- de_wp_d[longitude_idx, latitude_idx, , drop=FALSE]
tmin_d  <- tmin_d[longitude_idx, latitude_idx, , drop=FALSE]
tmax_d  <- tmax_d[longitude_idx, latitude_idx, , drop=FALSE]
nx <- length(longitude);  ny <- length(latitude)

# G) DEFINE DIMENSIONS AND VARIABLE METADATA
day_vals  <- as.numeric(days - days[1])
day_units <- sprintf("days since %s", days[1])
longitude_dim <- ncdim_def("longitude", units = longitude_units, vals = longitude)
latitude_dim <- ncdim_def("latitude", units = latitude_units, vals = latitude)
day_dim <- ncdim_def("time", units = day_units, vals = day_vals, unlim=TRUE)
vars_daily <- list(
  ncvar_def("rain", units = "mm", dim = list(longitude_dim, latitude_dim, day_dim), missval = NA, longname = "Daily total precipitation"),
  ncvar_def("srad", units = "MJ m-2", dim = list(longitude_dim, latitude_dim, day_dim), missval = NA, longname = "Daily total surface short-wave radiation"),
  ncvar_def("t2m",  units = "degree_Celsius", dim = list(longitude_dim, latitude_dim, day_dim), missval = NA, longname = "Daily mean 2-m temperature"),
  ncvar_def("dewp", units = "degree_Celsius", dim = list(longitude_dim, latitude_dim, day_dim), missval = NA, longname = "Daily mean 2-m dewpoint temperature"),
  ncvar_def("tmin", units = "degree_Celsius", dim = list(longitude_dim, latitude_dim, day_dim), missval = NA, longname = "Daily minimum 2-m temperature"),
  ncvar_def("tmax", units = "degree_Celsius", dim = list(longitude_dim, latitude_dim, day_dim), missval = NA, longname = "Daily maximum 2-m temperature")
)

# H) WRITE TO NETCDF
daily_file <- sprintf("aifs_%s_daily.nc", run_date)
ncw2 <- nc_create(daily_file, vars = vars_daily, force_v4 = TRUE)
ncvar_put(ncw2, "time", day_vals)
ncatt_put(ncw2, "time", "calendar",    "gregorian")
ncatt_put(ncw2, "time", "axis",        "T")
ncatt_put(ncw2, "time", "time_origin", as.character(days[1]))
ncvar_put(ncw2, "rain", rain_d)
ncvar_put(ncw2, "srad", srad_d)
ncvar_put(ncw2, "t2m",  t2m_d)
ncvar_put(ncw2, "dewp", de_wp_d)
ncvar_put(ncw2, "tmin", tmin_d)
ncvar_put(ncw2, "tmax", tmax_d)
nc_close(ncw2)
message("Daily‐aggregated, cropped file written to: ", daily_file)

# Delete Merged Weather
unlink(c(grb_paths, nc4_paths, outname))
