source("scripts/_functions.R")

#query_list <- s2_list(server = c("scihub", "gcloud"), tile = "15SXT", orbit = "126",
#                     time_interval = c("2018-07-02", "2019-07-01"))
save_path <- "GIS/sentinel2_data/"
## safe_is_online("C:/Users/nello/DOCUME~1/SEN2R~1/lta_orders/lta_20230303_153543.json")
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
file_names <- paste0(save_path, list.files(save_path, "B\\d\\d(_[1,2]0m)?\\.jp2$", recursive = T))
file_bands <- str_extract(file_names, "B\\d{2}")
file_dates <- str_extract(file_names, "\\d{8}") %>%
              as.Date("%Y%m%d")

file_path_table <- tibble(dates = file_dates, bands = file_bands, path = file_names)

date_range <- c("2018-01-01", "2020-07-01") ## period of raster images to select

used_files <- file_path_table %>%
              filter(between(dates, as.Date(date_range[1]), as.Date(date_range[2]))) %>%
              filter(!dates %in% clouded_dates$dates,
                     bands == "B08" |
                     bands == "B11" |
                     bands == "B12") %>% ## filtering clouds
              arrange(dates)

used_files <- used_files %>%
              distinct(dates, bands, .keep_all = T) %>% 
              filter(dates != as.Date("2018-06-05"))

reference_raster <- rast(used_files$path[1]) %>%
                    crop(area_of_interest) %>%
                    mask(area_of_interest)

rast(list_of_stacks[[1]], nlyrs = 3)

band_names <- unique(used_files$bands)
list_of_stacks <- list()
for (i in unique(as.character(used_files$dates))) {
  paths <- used_files %>%
           filter(dates == i) %>%
           select(path)
  raster <- resample_bands(paths$path, reference_raster, area_of_interest) 
  list_of_stacks[[i]] <- mask(raster, area_of_interest)
  crs(list_of_stacks[[i]]) <- crs
  names(list_of_stacks[[i]]) <- band_names
}

combined_rasters <- rast(list_of_stacks)
indices <- rep(c(1, 2, 3), nlyr(combined_rasters) / 3) ## indices of groups of layers

stats_table <- list()
bands <- unique(used_files$bands)
bands <- rep(bands, nlyr(combined_rasters) / length(bands))
stats_table <- tibble(dates = used_files$dates, band = bands,
                      means = global(combined_rasters, "mean", na.rm = T)$mean,
                      sds = global(combined_rasters, "sd", na.rm = T)$sd) %>%
               filter(means > 0)
stats_table$dates <- as.Date(stats_table$dates)

ggplot(stats_table, aes(x = dates, y = means)) +
       geom_point() + geom_line() + geom_smooth(se = F) +
       geom_label(aes(label = format.Date(dates, format = "%d"))) +
       facet_wrap(.~band) +
       scale_x_date(date_breaks = "1 month",
                    date_labels = "%b/%y")

## Manual selection of data
bare_soil_dates <- stats_table$dates

dates_range1 <- c("2018-03-22", "2018-05-16")

dates_range2 <- c("2018-11-17", "2019-01-31")

dates_range3 <- c("2019-03-01", "2019-05-15")

used_bare_soil_files1 <- used_files %>%
                         filter(dates %in% bare_soil_dates) %>%
                         filter(between(dates, as.Date(dates_range1[1]), as.Date(dates_range1[2])))
used_bare_soil_files2 <- used_files %>%
                         filter(dates %in% bare_soil_dates) %>%
                         filter(between(dates, as.Date(dates_range2[1]), as.Date(dates_range2[2])))
used_bare_soil_files3 <- used_files %>%
                         filter(dates %in% bare_soil_dates) %>%
                         filter(between(dates, as.Date(dates_range3[1]), as.Date(dates_range3[2])))
used_bare_soil_files <- rbind(used_bare_soil_files1, used_bare_soil_files2, used_bare_soil_files3)

band_names <- unique(used_bare_soil_files$bands)
list_of_stacks <- list()
for (i in unique(as.character(used_bare_soil_files$dates))) {
  paths <- used_files %>%
           filter(dates == i) %>%
           select(path)
  raster <- resample_bands(paths$path, reference_raster, area_of_interest) 
  list_of_stacks[[i]] <- mask(raster, area_of_interest)
  crs(list_of_stacks[[i]]) <- crs
  names(list_of_stacks[[i]]) <- band_names
}

combined_rasters <- rast(list_of_stacks)
indices <- rep(c(1, 2, 3), nlyr(combined_rasters) / 3) ## indices of groups of layers

band_means <- tapp(combined_rasters, indices, "mean")
names(band_means) <- c("B8", "B11", "B12")
band_vars <- tapp(combined_rasters, indices, var)
names(band_vars) <- c("B8", "B11", "B12")

write_paths_means <- paste0("GIS/bare_soil_temporal_bands/bare_soils_",
                           unique(used_bare_soil_files$bands),
                           "_means.tif")
write_paths_vars <- paste0("GIS/bare_soil_temporal_bands/bare_soils_",
                          unique(used_bare_soil_files$bands),
                          "_vars.tif")

writeRaster(band_means, write_paths_means)
writeRaster(band_vars, write_paths_vars)
