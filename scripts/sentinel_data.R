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
glm_files <- joined_file_path_table %>%
             filter(is.na(raster)) %>%
             select(dates, glm)
raster_files <- joined_file_path_table %>%
                select(dates, raster) %>%
                drop_na()
nrow(glm_files) + nrow(raster_files)

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
B08_seasonal_var <- ggplot(plot_data, aes(x = dates, y = means)) +
                    geom_point() + geom_line() +
                    geom_smooth(method = "lm", formula = y~poly(x, 15), se = F) +
                    geom_label(aes(label = format.Date(dates, format = "%d/%m"))) +
                    facet_wrap(.~band) +
                    scale_x_date(date_breaks = "1 year",
                                 date_labels = "%y")

## Comparing to RGB images to identify which ones are bare
rgb_file_names <- list.files(save_path, pattern = "TCI_10m.jp2", full.names = T, recursive = T)
rgb_file_dates <- str_extract(rgb_file_names, "\\d{8}") %>%
                  as.Date("%Y%m%d")
rgb_file_path_table <- tibble(dates = rgb_file_dates, path = rgb_file_names) %>%
                       distinct(dates, .keep_all = T)

list_of_rgb <- list()
for (i in unique(as.character(rgb_file_path_table$dates))) {
  paths <- rgb_file_path_table %>%
           filter(dates == i) %>%
           select(path)
  raster <- crop(rast(paths$path), area_of_interest) 
  list_of_rgb[[i]] <- mask(raster, area_of_interest)
  crs(list_of_rgb[[i]]) <- crs
  names(list_of_rgb[[i]]) <- c("Red", "Green", "Blue")
}

plot_dates <- c("2018-12-22", "2019-07-25", "2020-01-21", "2020-07-14", "2020-12-26",
                "2021-07-24", "2021-11-26", "2022-07-19", "2023-01-05")
rgb_plot_list <- list()
count <- 1
for (i in plot_dates) {
  plot <- ggplot() + ylab("") + xlab("") +
          ggtitle(paste("Date:", i)) +
          geom_spatial_rgb(data = list_of_rgb[[i]] %>% as.data.frame(xy = T),
                           aes(x = x, y = y, r = Red, g = Green, b = Blue)) +
          coord_fixed(ratio = 1) +
          theme(panel.grid = element_blank(),
               axis.text.y = element_text(angle = 90,hjust = 0.5))
  rgb_plot_list[[count]] <- plot
  count <- count + 1
}
B08_seasonal_var / wrap_plots(rgb_plot_list, ncol = 9) + plot_layout(nrow = 2)

## Testing the quality of selected bare soil images
bare_soil_dates <- stats_table %>%
                   filter(means < 3000 & means > 0) %>%
                   select(dates) %>%
                   unique()

## Randomly iterating and selecting bad data
manual_bad_raster_selection <- c("2019-09-08", "2022-07-04", "2022-06-14", "2020-07-40", "2019-06-10",
                                 "2020-01-16", "2021-07-19", "2020-05-05", "2021-04-30", "2018-12-27",
                                 "2020-02-10", "2022-06-24", "2020-07-14", "2021-12-21", "2022-04-20",
                                 "2019-06-05", "2022-04-20", "2021-07-04", "2020-08-18", "2019-09-18",
                                 "2022-02-24", "2022-08-08", "2022-06-19", "2022-06-29", "2021-07-29",
                                 "2022-08-13", "2021-02-19", "2019-07-20", "2020-07-04", "2020-09-07",
                                 "2021-03-01", "2021-09-17", "2019-06-20", "2022-07-19", "2019-07-05",
                                 "2021-07-14", "2019-06-20", "2022-07-19", "2021-05-30", "2019-07-05",
                                 "2021-07-14", "2022-07-09", "2022-06-04", "2019-07-25", "2022-09-12",
                                 "2020-06-14", "2021-07-24", "2021-11-01", "2022-09-27", "2021-10-22",
                                 "2021-11-26", "2022-01-10", "2020-03-11", "2021-02-04", "2019-09-13",
                                 "2019-01-31", "2022-09-22", "2020-03-01", "2019-05-11", "2022-10-17",
                                 "2019-08-29", "2021-06-14", "2022-10-07", "2021-12-01", "2022-10-02",
                                 "2021-05-15")

rgb_file_path_table <- rgb_file_path_table %>% filter(dates %in% bare_soil_dates$dates) %>%
                       filter(!dates %in% as.Date(manual_bad_raster_selection))

list_of_rgb <- list()
for (i in unique(as.character(rgb_file_path_table$dates))) {
  paths <- rgb_file_path_table %>%
           filter(dates == i) %>%
           select(path)
  raster <- crop(rast(paths$path), area_of_interest) 
  list_of_rgb[[i]] <- mask(raster, area_of_interest)
  crs(list_of_rgb[[i]]) <- crs
  names(list_of_rgb[[i]]) <- c("Red", "Green", "Blue")
}

plot_dates <- names(list_of_rgb)
rgb_plot_list <- list()
for (i in plot_dates) {
  plot <- ggplot() + ylab("") + xlab("") +
          ggtitle(paste("Date:", i)) +
          geom_spatial_rgb(data = list_of_rgb[[i]] %>% as.data.frame(xy = T),
                           aes(x = x, y = y, r = Red, g = Green, b = Blue)) +
          coord_fixed(ratio = 1) +
          theme(panel.grid = element_blank(),
               axis.text.y = element_text(angle = 90,hjust = 0.5))
  rgb_plot_list[[i]] <- plot
}
random_dates <- sample(plot_dates, 27, replace = F)
wrap_plots(rgb_plot_list[random_dates], ncol = 9)

definitive_bare_soil_table <- used_files %>%
                              filter(as.character(dates) %in% plot_dates)
definitive_bare_soil_dates <- definitive_bare_soil_table$dates %>% as.character() %>% unique()
list_of_bare_soils <- list_of_stacks[definitive_bare_soil_dates]

number_of_bands <- length(unique(definitive_bare_soil_table$bands))
bare_soil_combined_rasters <- rast(list_of_bare_soils)

stats_table <- list()
bands <- unique(definitive_bare_soil_table$bands)
bands <- rep(bands, nlyr(combined_rasters) / length(bands))
stats_table <- tibble(dates = definitive_bare_soil_table$dates, band = bands,
                      means = global(combined_rasters, "mean", na.rm = T)$mean,
                      sds = global(combined_rasters, "sd", na.rm = T)$sd) %>%
               filter(means > 0)
stats_table$dates <- as.Date(stats_table$dates)

## Removing outliers (possible mask problems)
outliers <- boxplot(stats_table$means)$out
stats_table <- stats_table %>%
               filter(!means %in% outliers)
write_csv(stats_table, "tables/sentine2_bands_stats/bare_soil_bands_stats.csv")

plot_data <- stats_table %>%
             filter(band == "B08")
B08_bare_soil_seasonal_var <- ggplot(plot_data, aes(x = dates, y = means)) +
                              geom_point() + geom_line() +
                              geom_smooth(method = "lm", formula = y~poly(x, 15), se = F) +
                              geom_label(aes(label = format.Date(dates, format = "%d/%m"))) +
                              facet_wrap(.~band) +
                              scale_x_date(date_breaks = "1 year",
                                           date_labels = "%y")

## Calculating mean and var rasters

list_of_stat_plots <- list()
for (i in band_names) {
  ## All files (except clouds)
  indices <- names(combined_rasters) %>% ## indices of groups of layers
             str_extract(pattern = "\\d$") %>%
             as.numeric()
  
  band_means <- tapp(combined_rasters, indices, "mean")
  names(band_means) <- band_names
  band_vars <- tapp(combined_rasters, indices, "var")
  names(band_vars) <- band_names
  
  all_data_mean_plot <- ggplot() + ylab("") + xlab("") +
                        ggtitle(paste("All data - n =", used_files$dates %>% unique() %>% length())) +
                        geom_raster(data = band_means %>% as.data.frame(xy = T),
                                         aes_string(x = "x", y = "y", fill = i)) +
                        scale_fill_distiller(type = "div", palette = "Spectral", na.value = NA) +
                        coord_fixed(ratio = 1) +
                        theme(panel.grid = element_blank(),
                             axis.text.y = element_text(angle = 90,hjust = 0.5))
  
  all_data_dens <- ggplot(band_means %>% as.data.frame()) +
                   ylab("") + xlab("") + ggtitle("All data") +
                   geom_density(aes_string(x = i))
  
  ## Bare soils
  indices <- names(bare_soils_combined_rasters) %>% ## indices of groups of layers
             str_extract(pattern = "\\d$") %>%
             as.numeric()
  
  bare_soil_band_means <- tapp(bare_soil_combined_rasters, indices, "mean")
  names(bare_soil_band_means) <- band_names
  bare_soil_band_vars <- tapp(bare_soil_combined_rasters, indices, "var")
  names(bare_soil_band_vars) <- band_names
  
  bare_soil_mean_plot <- ggplot() + ylab("") + xlab("") +
                         ggtitle(paste("Bare soil data - n =", definitive_bare_soil_dates %>% length())) +
                         geom_raster(data = bare_soil_band_means %>% as.data.frame(xy = T),
                                          aes_string(x = "x", y = "y", fill = i)) +
                         scale_fill_distiller(type = "div", palette = "Spectral", na.value = NA) +
                         coord_fixed(ratio = 1) +
                         theme(panel.grid = element_blank(),
                              axis.text.y = element_text(angle = 90,hjust = 0.5))
  
  bare_soil_dens <- ggplot(bare_soil_band_means %>% as.data.frame()) +
                    ylab("") + xlab("") + ggtitle("Bare soil data") +
                    geom_density(aes_string(x = i))
  
  list_of_stat_plots[[i]] <- (all_data_mean_plot + bare_soil_mean_plot) / ( all_data_dens + bare_soil_dens)
}

list_of_stat_plots[[1]]
list_of_stat_plots[[2]]
list_of_stat_plots[[3]]
list_of_stat_plots[[4]]
list_of_stat_plots[[5]]
list_of_stat_plots[[6]]

write_paths_means <- paste0("GIS/bare_soil_temporal_bands/bare_soils_",
                           unique(used_bare_soil_files$bands),
                           "_means.tif")
write_paths_vars <- paste0("GIS/bare_soil_temporal_bands/bare_soils_",
                          unique(used_bare_soil_files$bands),
                          "_vars.tif")

writeRaster(band_means, write_paths_means, overwrite = T)
writeRaster(band_vars, write_paths_vars, overwrite = T)
