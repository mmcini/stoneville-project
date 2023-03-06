source("scripts/_functions.R")

days_difference <- function(dates) {
  days_diff <- c(0)
  for (i in seq(2, length(dates), 1)) {
    diff <- dates[i] - dates[i - 1]
    days_diff <- rbind(days_diff, diff)
  }
  return(days_diff)
}

#query_list <- s2_list(server = c("scihub", "gcloud"), tile = "15SXT", orbit = "126",
#                     time_interval = c("2018-07-02", "2019-07-01"))
save_path <- "D:/Arquivos/temp_rasters_stoneville/"
## safe_is_online("C:/Users/nello/DOCUME~1/SEN2R~1/lta_orders/lta_20230303_153543.json")
#s2_download(query_list, outdir = save_path)

#list <- safe_is_online("C:/Users/nello/DOCUME~1/SEN2R~1/lta_orders/lta_20230303_153543.json")

#query_list <- s2_list(server = c("scihub"), tile = "15SXT", orbit = "126",
#                      time_interval = c("2019-06-15"))
#save_path <- "D:/Arquivos/temp_rasters_stoneville/"
#s2_download(query_list, outdir = save_path)

## Loading cloud masks
file_names <- paste0(save_path, list.files(save_path, "MSK_CLOUDS_B00.gml$", recursive = T))
file_dates <- str_extract(file_names, "\\d{8}") %>%
              as.Date("%Y%m%d")

file_path_table <- tibble(dates = file_dates, path = file_names)

crs <- "+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs +type=crs"
area_of_interest <- st_read("GIS/area_outline/SWMRU_extent.shp")
area_of_interest <- st_transform(area_of_interest, crs)

list_of_masks <- list()
for (i in as.character(file_path_table$dates)) {
  paths <- file_path_table %>%
    filter(dates == i) %>%
    select(path)
  try(
    list_of_masks[[i]] <- st_intersection(st_read(paths[[1]]), area_of_interest),
    silent = F
  )
}

total_area <- st_area(area_of_interest)
dates_with_clouds <- names(list_of_masks)
cloud_masks_table <- tibble(dates = character(), area = numeric())
for (i in dates_with_clouds) {
  cloud_masks_table <- cloud_masks_table %>%
                       add_row(dates = i,
                               area = sum(as.numeric(st_area(list_of_masks[[i]]$extentOf))))
}
cloud_masks_table$dates <- as.Date(cloud_masks_table$dates, format= "%Y-%m-%d")
cloud_masks_table <- cloud_masks_table %>%
                     mutate(area_percent = (area / as.numeric(total_area)) * 100)

## Images with too many clouds
clouded_dates <- cloud_masks_table %>%
                 filter(area_percent > 0)

## Loading rasters
file_names <- paste0(save_path, list.files(save_path, "B\\d\\d(_10m)?\\.jp2$", recursive = T))
file_bands <- str_extract(file_names, "B\\d{2}")
file_dates <- str_extract(file_names, "\\d{8}") %>%
              as.Date("%Y%m%d")

file_path_table <- tibble(dates = file_dates, bands = file_bands, path = file_names)

date_range <- c("2018-01-01", "2020-07-01") ## period of raster images to select

used_files <- file_path_table %>%
              filter(between(dates, as.Date(date_range[1]), as.Date(date_range[2]))) %>%
              filter(!dates %in% clouded_dates$dates, bands == "B08") %>% ## filtering clouds
              arrange(dates)

used_files <- used_files %>%
              filter(dates != as.Date("2018-06-05"))

used_files <- used_files %>%
              mutate(days_diff = days_difference(used_files$dates))

list_of_bricks <- list()
for (i in unique(as.character(used_files$dates))) {
  paths <- used_files %>%
           filter(dates == i, bands == "B08") %>%
           select(path)
  list_of_bricks[[i]] <- mask(stack(paths[[1]]), area_of_interest)
  crs(list_of_bricks[[i]]) <- crs
  len <- length(names(list_of_bricks[[i]])) ## some bands have more than 1 layer, idk why!!
  if (len > 1) {band_names <- rep("B8", len)} else {band_names <- "B8"} ## this fixes the multilayer problem
  names(list_of_bricks[[i]]) <- band_names
}

stats_table <- tibble(dates = character(), means = numeric(), sds = numeric(), ranges = numeric())
for (i in names(list_of_bricks)) {
  raster <- list_of_bricks[[i]]
  mean <- global(as(raster, "SpatRaster"), "mean", na.rm = T)
  sd <- global(as(raster, "SpatRaster"), "sd", na.rm = T)
  range <- global(as(raster, "SpatRaster"), "range", na.rm = T)
  stats_table <- stats_table %>%
                 add_row(dates = i, means = mean$mean, sds = sd$sd, ranges = range$range)
}

stats_table$dates <- as.Date(stats_table$dates)

stats_table <- read_csv("tables/nir_test.csv")
stats_table <- stats_table %>%
               filter(means != 0)

ggplot(stats_table, aes(x = dates, y = means)) +
       geom_point() + geom_line() + geom_smooth(se = F) +
       geom_label(aes(label = format.Date(dates, format = "%d"))) +
       scale_x_date(date_breaks = "1 month",
                    date_labels = "%b/%y")


## Temporary manual selection of data
bare_soil_dates <- stats_table

dates_range1 <- c("2018-03-22", "2018-05-16")

dates_range2 <- c("2018-11-17", "2019-01-31")

dates_range3 <- c("2019-03-01", "2019-05-15")

used_bare_soil_files1 <- used_files %>%
                        filter(dates %in% bare_soil_dates$dates) %>%
                        filter(between(dates, as.Date(dates_range1[1]), as.Date(dates_range1[2])))
used_bare_soil_files2 <- used_files %>%
                        filter(dates %in% bare_soil_dates$dates) %>%
                        filter(between(dates, as.Date(dates_range2[1]), as.Date(dates_range2[2])))
used_bare_soil_files3 <- used_files %>%
                        filter(dates %in% bare_soil_dates$dates) %>%
                        filter(between(dates, as.Date(dates_range3[1]), as.Date(dates_range3[2])))
used_bare_soil_files <- rbind(used_bare_soil_files1, used_bare_soil_files2, used_bare_soil_files3)

list_of_bricks <- list()
for (i in unique(as.character(used_bare_soil_files$dates))) {
  paths <- used_files %>%
           filter(dates == i, bands == "B08") %>%
           select(path)
  raster <- crop(stack(paths[[1]]), area_of_interest)
  list_of_bricks[[i]] <- mask(raster, area_of_interest) ##CROP BEFORE DOING THIS!!
  crs(list_of_bricks[[i]]) <- crs
  len <- length(names(list_of_bricks[[i]])) ## some bands have more than 1 layer, idk why!!
  if (len > 1) {band_names <- rep("B8", len)} else {band_names <- "B8"} ## this fixes the multilayer problem
  names(list_of_bricks[[i]]) <- band_names
}

bare_soil_rasters <- stack(list_of_bricks)

bare_soil_rasters_mean <- app(as(bare_soil_rasters, "SpatRaster"), fun = "mean", rm.ma = T)

bare_soil_rasters_mean <- rast("GIS/band_mean/bare_soils_B8_mean.tif")

plot(bare_soil_rasters_mean)