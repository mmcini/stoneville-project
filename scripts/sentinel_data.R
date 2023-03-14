source("scripts/_functions.R")

#query_list <- s2_list(server = c("gcloud"), tile = "15SXT", orbit = "126", level = "L2A",
#                      time_interval = c("2019-07-01", "2020-07-10"))
save_path <- "GIS/sentinel2_data"
#s2_download(query_list, outdir = save_path, abort = F)

## Cloud masks
## Some files have cloud percentage rasters, some have gml, some dont?
## Comparing dates and searching which scenes have raster or shape cloud mask files
all_file_names <- list.files(save_path, pattern = "SAFE", full.names = T)
all_file_dates <- str_extract(all_file_names, "\\d{8}") %>%
                  as.Date("%Y%m%d")
all_file_path_table <- tibble(dates = all_file_dates, path = all_file_names) %>%
                       distinct(dates, .keep_all = T)

gml_file_names <- list.files(save_path, "MSK_.*_B00.gml$", full.names = T, recursive = T)
gml_file_dates <- str_extract(gml_file_names, "\\d{8}") %>%
                  as.Date("%Y%m%d")
gml_file_path_table <- tibble(dates = gml_file_dates, path = gml_file_names) %>%
                       distinct(dates, .keep_all = T)

jp2_file_names <- list.files(save_path, "MSK_CLDPRB_20m", full.names = T, recursive = T)
jp2_file_dates <- str_extract(jp2_file_names, "\\d{8}") %>%
                  as.Date("%Y%m%d")
jp2_file_path_table <- tibble(dates = jp2_file_dates, path = jp2_file_names) %>%
                       distinct(dates, .keep_all = T)

## Finding out which files have glms and which have rasters
joined_file_path_table <- full_join(all_file_path_table, gml_file_path_table, by = c("dates" = "dates")) %>%
                          full_join(jp2_file_path_table, by = c("dates" = "dates"))
names(joined_file_path_table) <- c("dates", "all_files", "glm", "raster")

## Selecting files separately
raster_files <- joined_file_path_table %>%
                filter(is.na(glm)) %>%
                select(dates, raster)
glm_files <- joined_file_path_table %>%
             select(dates, glm) %>%
             drop_na()

## loading area of interest
crs <- "+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs +type=crs"
area_of_interest <- st_read("GIS/area_outline/SWMRU_extent.shp")
area_of_interest <- st_transform(area_of_interest, crs)

## creating glm cloud masks
list_of_masks_from_glm <- list()
for (i in as.character(glm_files$dates)) {
  paths <- glm_files %>%
    filter(dates == i) %>%
    select(glm)
  try(
  list_of_masks_from_glm[[i]] <- st_intersection(st_read(paths[[1]]), area_of_interest, quite = T), 
  silent = T
  )
}

## creating raster cloud masks
list_of_masks_from_raster <- list()
for (i in as.character(raster_files$dates)) {
  cloud_probability <- 80
  paths <- raster_files %>%
    filter(dates == i) %>%
    select(raster)
  ## tries to do it
  attempt <- tryCatch(
    raster <- rast(paths$raster) %>%
            crop(area_of_interest) %>%
            mask(area_of_interest),
    error = function (e) e
    )
  ## if class is not "error", continue
  if (!inherits(attempt, "error")) {
    raster[raster < cloud_probability] <- NA
    raster[raster >= cloud_probability] <- 1
    shape <- as.polygons(raster)
    list_of_masks_from_raster[[i]] <- st_as_sf(shape)
  }
}

## Merging masks and calculating cloud cover percentage
list_of_masks <- c(list_of_masks_from_glm, list_of_masks_from_raster)
total_area <- st_area(area_of_interest)
dates_with_clouds <- names(list_of_masks)
cloud_masks_table <- tibble(dates = character(), area = numeric())
for (i in dates_with_clouds) {
  ## if this doesnt work (does not come from glms)
  attempt <- tryCatch(
    cloud_masks_table <- cloud_masks_table %>%
                         add_row(dates = i,
                                 area = sum(as.numeric(st_area(list_of_masks[[i]]$extentOf)))),
    error = function (e) e
  )
  ## do this (comes from rasters)
  if (inherits(attempt, "error")) {
    cloud_masks_table <- cloud_masks_table %>%
                         add_row(dates = i,
                                 area = sum(as.numeric(st_area(list_of_masks[[i]]))))
  }
}
cloud_masks_table$dates <- as.Date(cloud_masks_table$dates, format= "%Y-%m-%d")
cloud_masks_table <- cloud_masks_table %>%
                     mutate(area_percent = (area / as.numeric(total_area)) * 100)

## Images with too many clouds
clouded_dates <- cloud_masks_table %>%
                 filter(area_percent > 0)

## Loading rasters
file_names <- list.files(save_path, "B\\d\\d(_[1,2]0m)?\\.jp2$", full.names = T, recursive = T)
file_bands <- str_extract(file_names, "B\\d{2}")
file_dates <- str_extract(file_names, "\\d{8}") %>%
              as.Date("%Y%m%d")

file_path_table <- tibble(dates = file_dates, bands = file_bands, path = file_names)

used_files <- file_path_table %>%
              filter(!dates %in% clouded_dates$dates,
                     bands == "B02" |
                     bands == "B03" |
                     bands == "B04" |
                     bands == "B08" |
                     bands == "B11" |
                     bands == "B12") %>% ## filtering clouds
              arrange(dates)

used_files <- used_files %>% # selecting only one of each band (some have 10 m, 20 m, etc.)
              distinct(dates, bands, .keep_all = T)

reference_raster <- rast(used_files$path[1]) %>%
                    crop(area_of_interest) %>%
                    mask(area_of_interest)

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

number_of_bands <- length(unique(used_files$bands))
combined_rasters <- rast(list_of_stacks)
indices <- rep(1:number_of_bands, nlyr(combined_rasters) / 3) ## indices of groups of layers

stats_table <- list()
bands <- unique(used_files$bands)
bands <- rep(bands, nlyr(combined_rasters) / length(bands))
stats_table <- tibble(dates = used_files$dates, band = bands,
                      means = global(combined_rasters, "mean", na.rm = T)$mean,
                      sds = global(combined_rasters, "sd", na.rm = T)$sd) %>%
               filter(means > 0)
stats_table$dates <- as.Date(stats_table$dates)

## Removing outliers (possible mask problems)
outliers <- boxplot(stats_table$means)$out
stats_table <- stats_table %>%
               filter(!means %in% outliers)
write_csv(stats_table, "tables/sentine2_bands_stats/bands_stats.csv")

plot_data <- stats_table %>%
             filter(band == "B08")
ggplot(plot_data, aes(x = dates, y = means)) +
       geom_point() + geom_line() + geom_smooth(se = F) +
       geom_label(aes(label = format.Date(dates, format = "%d"))) +
       facet_wrap(.~band) +
       scale_x_date(date_breaks = "1 year",
                    date_labels = "%y")

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

## Correcting L1C images ###########################################################################
l1c_files <-str_extract(used_bare_soil_files$path, pattern = "GIS.*L1C.*SAFE") %>%
            unique()
sen2cor(l1c_prodlist = l1c_files, outdir = "GIS/sentinel2_data/processed_data")

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
indices <- names(combined_rasters) %>% ## indices of groups of layers
           str_extract(pattern = "\\d$") %>%
           as.numeric()

band_means <- tapp(combined_rasters, indices, "mean")
names(band_means) <- band_names
band_vars <- tapp(combined_rasters, indices, var)
names(band_vars) <- band_names

write_paths_means <- paste0("GIS/bare_soil_temporal_bands/bare_soils_",
                           unique(used_bare_soil_files$bands),
                           "_means.tif")
write_paths_vars <- paste0("GIS/bare_soil_temporal_bands/bare_soils_",
                          unique(used_bare_soil_files$bands),
                          "_vars.tif")

writeRaster(band_means, write_paths_means, overwrite = T)
writeRaster(band_vars, write_paths_vars, overwrite = T)
