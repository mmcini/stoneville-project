source("scripts/_functions.R")

# Loading the data #################################################################################
raw_soil_data <- read_sf("GIS/area_outline/swmru_datasites_point_proj_utm15.shp")
raw_soil_data <- raw_soil_data %>%
                 select(grep("PPM_.", colnames(raw_soil_data)), "OM_PERCENT", "CEC", "WATER_PH")
crs <- CRS("+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs +type=crs")
raw_soil_data <- st_transform(raw_soil_data, crs)
area_of_interest <- st_read("GIS/area_outline/SWMRU_extent.shp")
area_of_interest <- st_transform(area_of_interest, crs)

files <- list.files("GIS/sentinel2_data", pattern = ".*B\\d{2}.*tiff$")
dates <- as.Date(str_extract(files, "\\d{4}-\\d{2}-\\d{2}"))
bands <- str_extract(files, "B\\d{2}")
raster_files <- tibble(files = paste0("GIS/sentinel2_data/", files),
                       dates = dates,
                       bands = bands)

## Getting raster paths
selected_dates <- read_csv("tables/nir_test.csv")
file_names <- paste0(save_path, list.files(save_path, "B\\d\\d(_[1,2]\\dm)?\\.jp2$", recursive = T))
file_bands <- str_extract(file_names, "B\\d{2}")
file_dates <- str_extract(file_names, "\\d{8}") %>%
              as.Date("%Y%m%d")

file_path_table <- tibble(dates = file_dates, bands = file_bands, path = file_names)

dem <- raster("GIS/DEM/DEM_1m_small_clip_USGS_3DEP_stoneville.tif")
dem <- projectRaster(dem, crs = crs)

## Comparing single rasters and the bare soil mean #################################################
## Selecting rasters (based on sentinel_data.R)
selected_files <- file_path_table %>%
                  filter(file_path_table$dates %in% selected_dates$dates) %>% 
                  filter(dates == as.Date("2019-01-26") | dates == as.Date("2018-12-17"),
                         bands == "B08") # selecting two bare soil images for comparison

single_bare_soil_raster_2018 <- raster(selected_files$path[1]) %>%
                                crop(area_of_interest) %>%
                                mask(area_of_interest)
single_bare_soil_raster_2019 <- raster(selected_files$path[2]) %>%
                                crop(area_of_interest) %>%
                                mask(area_of_interest)
mean_bare_soil_raster <- raster("GIS/band_mean/bare_soils_B8_mean.tif")

single_data_2018 <- single_bare_soil_raster_2018 %>%
                    as.data.frame(xy = T)
names(single_data_2018) <- c("x", "y", "single_B8")
single_data_2019 <- single_bare_soil_raster_2019 %>%
                    as.data.frame(xy = T)
names(single_data_2019) <- c("x", "y", "single_B8")
mean_data <- mean_bare_soil_raster %>%
             as.data.frame(xy = T)
names(mean_data) <- c("x", "y", "Mean 2018-2019")

single_mean_comparison_data <- cbind(mean_data,
                                     "2018-12-17" = single_data_2018$single_B8,
                                     "2019-02-26" = single_data_2019$single_B8) %>%
                               pivot_longer("Mean 2018-2019":"2019-02-26", names_to = "variables",
                                            values_to = "Reflectance")

bare_soil_raster <- ggplot() + ylab("") + xlab("") +
                    geom_raster(data = single_mean_comparison_data, aes(x = x, y = y, fill = Reflectance)) +
                    facet_wrap(. ~ variables) +
                    scale_fill_distiller(type = "div", palette = "Spectral", na.value = NA) +
                    theme(panel.grid = element_blank())

bare_soil_dens <- ggplot(single_mean_comparison_data) +
                  ylab("") + xlab("") +
                  geom_density(aes(x = Reflectance)) +
                  facet_wrap(. ~ variables)

bare_soil_raster + bare_soil_dens + plot_layout(nrow = 2)

## Loading rasters for correlation analysis ########################################################
selected_files <- file_path_table %>%
                  filter(file_path_table$dates %in% selected_dates$dates) %>% 
                  filter(dates == as.Date("2018-12-17")) %>% # selecting one bare soil image
                  slice(c(1:4, 8:12)) %>%
                  arrange(bands)

## Band 8 raster to use as reference to resample others
reference_raster <- raster(selected_files$path[7])

unique_dates <- unique(selected_files$dates)
list_of_bricks <- list()
for(i in as.character(unique_dates)) {
  paths <- selected_files %>%
           filter(dates == i) %>%
           select(path)
  rasters <- resample_bands(paths$path, reference_raster, area_of_interest)
  list_of_bricks[[i]] <- mask(rasters, area_of_interest)
  crs(list_of_bricks[[i]]) <- crs
}

names(list_of_bricks$`2018-12-17`) <- selected_files$bands

export <- rast(list_of_bricks$`2018-12-17`)
writeRaster(export, "GIS/sentinel2_data.tif")

writeRaster

# Extracting data from rasters #####################################################################
bands_2017 <- list_of_bricks[["2017-05-01"]] %>%
              terra::mask(area_of_interest)
ndmi <- ndmi_sentinel_2(bands_2017)
ndwi <- ndwi_sentinel_2(bands_2017)
resamp_dem <- terra::resample(dem, bands_2017)

cropped_bands_2017 <- stack(bands_2017, resamp_dem, ndmi, ndwi)
names(cropped_bands_2017) <- c(names(bands_2017), "elevation", "ndmi", "ndwi")

raster_soil_data <- cbind(raw_soil_data, raster::extract(cropped_bands_2017, raw_soil_data))

# Descriptive stats ################################################################################
pivoted_data <- raw_soil_data %>%
                select(grep("PPM_.", colnames(raw_soil_data)), "OM_PERCENT", "CEC", "WATER_PH") %>%
                pivot_longer(cols = "PPM_P":"WATER_PH",
                             names_to = "variables",
                             values_to = "values")

ggplot(pivoted_data, aes(y = values)) +
  geom_boxplot() +
  facet_wrap(variables ~ ., scales = "free")

ggplot(pivoted_data, aes(x = values)) +
  geom_density() +
  facet_wrap(variables ~ ., scales = "free")

## Correlations
cor <- raster_soil_data %>%
        as_tibble() %>% 
        select(-geometry) %>%
        drop_na() %>%
        cor(method = "spearman")
cor_pvalue <- raster_soil_data %>%
              as_tibble() %>%
              select(-geometry) %>%
              drop_na() %>%
              cor.mtest(method = "spearman", conf.level = 0.95) # calculates p-value

par(family = "Times New Roman")
corrplot::corrplot(cor, method = "color", order = "hclust",
                   p.mat = cor_pvalue$p, addrect = 2,
                   tl.col = "black", pch.cex = 1, sig.level = c(0.001, 0.01, 0.05),
                   insig = "label_sig")

target_vars <- c("OM_PERCENT", "CEC", "PPM_P", "PPM_K", "PPM_MG",
                 "PPM_CA", "PPM_S", "PPM_ZN", "WATER_PH")
cor_plot_data <- cor %>%
  as_tibble() %>%
  add_column(variables = rownames(cor)) %>%
  select("variables", target_vars) %>%
  filter(!variables %in% target_vars) %>%
  mutate(variables = factor(variables, levels = variables, ordered = T))

pvalue_plot_data <- cor_pvalue$p %>%
                    as_tibble() %>%
                    mutate(across(everything(), assign_significance)) %>%
                    add_column(variables = rownames(cor)) %>%
                    select("variables", target_vars) %>%
                    filter(!variables %in% target_vars) %>%
                    mutate(variables = factor(variables, levels = variables, ordered = T))

list_of_cor_plots <- list()
for (i in names(cor_plot_data[-1])) {
  plot <- ggplot(cor_plot_data, aes_string(x = "variables", y = i)) +
    geom_histogram(stat = "identity") + geom_label(aes(label = pvalue_plot_data[[i]])) +
    geom_hline(yintercept = -0.5) +
    scale_y_continuous(limits = c(-1, 0.1), breaks = seq(-1, 1, 0.25))
  list_of_cor_plots[[i]] <- ggplot_gtable(ggplot_build(plot)) # to solve lazy evaluation  
}
ggarrange(plotlist = list_of_cor_plots)

list_of_cor_plots <- list()
for (i in names(cor_plot_data[-1])) {
  data <- raster_soil_data %>% as_tibble() %>% select(-geometry)
  test <- cbind(data[,!names(data) %in% names(cor_plot_data[-1])], OM = data[[i]])
  prcc <- epi.prcc(test) %>%
    mutate(var = factor(var, levels = var, ordered = T))
  plot <- ggplot(prcc, aes_string(x = "var", y = "est")) + ylab(i) +
    geom_histogram(stat = "identity") + geom_label(aes(label = assign_significance(prcc$p.value))) +
    scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.5, 0.5, 0.25))
  list_of_cor_plots[[i]] <- ggplot_gtable(ggplot_build(plot)) # to solve lazy evaluation  
}
ggarrange(plotlist = list_of_cor_plots)

# Modeling #########################################################################################
target_vars <- c("OM_PERCENT", "CEC", "PPM_K")
explanatory_vars <- c("elevation", "ndmi", "ndwi", "B3", "B8", "B11")
data <- raster_soil_data %>%
        as_tibble()

## Cross validation models and scores
cv_models_list <- build_cv_models(data, target_vars, explanatory_vars, seed = 200)
model_scores_cv_figures(cv_models_list)
model_scores_tables(cv_models_list, "rf_scores")

## Models using all data available
models_list <- build_models(data, target_vars, explanatory_vars, seed = 200)

# Predicting rasters ###############################################################################
df_cropped_bands_2017 <- cropped_bands_2017 %>%
                         as.data.frame(xy = T) %>%
                         drop_na()

preproc <- preProcess(df_cropped_bands_2017[-c(1, 2)], method = c("center", "scale"))
preprocessed_cropped_bands_2017 <- predict(preproc, df_cropped_bands_2017) %>%
                                   rasterFromXYZ(crs = crs)

build_rasters(preprocessed_cropped_bands_2017, models_list, "rf_model_outputs")

# grid sampling tests ##############################################################################
data <- raster_soil_data %>%
        as("Spatial")
crs(data) <- crs

sample_progression <- c(seq(20, 100, 10), seq(150, 500, 50), seq(600, 1000, 100))
regular_models <- build_sampling_models(data, target_vars, explanatory_vars, area_of_interest,
                                        "regular", preprocessed_cropped_bands_2017,
                                        n_samples = sample_progression, crs = crs)

## Checking results
regular_results <- read_csv("tables/regular_grid_results/model_scores_table.csv")

list_of_plots <- list()
for (i in unique(regular_results$variable)) {
  regular_plot_data <- regular_results %>%
                       filter(variable == i)
  
  plot <- ggplot(regular_plot_data, aes(x = n_samples, y = R2_valid)) +
    geom_line() + ggtitle(i) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  list_of_plots[[i]] <- plot
}
regular_r2_plots <- ggarrange(plotlist = list_of_plots, ncol = 1)

## Comparing with ground truth
comparison_plots <- raster_comparisons("regular")

vars <- str_extract(names(comparison_plots), "(?<=samples_).*") %>%
        unique()
for (i in vars) {
  path <- paste0("figures/regular_comparison/", i, "_comparisons.png")
  wrap_plots(plotlist = comparison_plots[grepl(i, names(comparison_plots))], ncol = 1)
  ggsave(path, device = "png", bg = "white", units = "mm",
         dpi = 300, width = 300, height = 100 * length(sample_progression),
         limitsize = F)
}

## Comparison summary
regular_summary <- read_csv("tables/regular_grid_results/comparison_summary.csv")

list_of_plots <- list()
for (i in unique(regular_summary$var)) {
  regular_plot_data <-  regular_summary %>%
    filter(var == i)
  
  plot <- ggplot(regular_plot_data, aes(x = n, y = ov)) +
    geom_line() + geom_smooth(se = F) + geom_hline(yintercept = 0.95) + ggtitle(i) +
    scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.05))
  list_of_plots[[i]] <- plot
}
regular_plot_comparison <- ggarrange(plotlist = list_of_plots, ncol = 1)

# random sampling tests ############################################################################
data <- raster_soil_data %>%
        as("Spatial")
crs(data) <- crs

sample_progression <- c(seq(20, 100, 10), seq(150, 500, 50), seq(600, 1000, 100))
random_models <- build_sampling_models(data, target_vars, explanatory_vars, area_of_interest,
                                       "random", preprocessed_cropped_bands_2017,
                                       n_samples = sample_progression, crs = crs)

## Checking results
random_results <- read_csv("tables/random_grid_results/model_scores_table.csv")

list_of_plots <- list()
for (i in unique(random_results$variable)) {
  random_plot_data <- random_results %>%
                       filter(variable == i)
  
  plot <- ggplot(random_plot_data, aes(x = n_samples, y = R2_valid)) +
    geom_line() + ggtitle(i) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  list_of_plots[[i]] <- plot
}
random_r2_plots <- ggarrange(plotlist = list_of_plots, ncol = 1)

## Comparing with ground truth
comparison_plots <- raster_comparisons("random")

vars <- str_extract(names(comparison_plots), "(?<=samples_).*") %>%
        unique()
for (i in vars) {
  path <- paste0("figures/random_comparison/", i, "_comparisons.png")
  wrap_plots(plotlist = comparison_plots[grepl(i, names(comparison_plots))], ncol = 1)
  ggsave(path, device = "png", bg = "white", units = "mm",
         dpi = 300, width = 300, height = 100 * length(sample_progression),
         limitsize = F)
}

## Comparison summary
random_summary <- read_csv("tables/random_grid_results/comparison_summary.csv")

list_of_plots <- list()
for (i in unique(random_summary$var)) {
  random_plot_data <-  random_summary %>%
                       filter(var == i)
  
  plot <- ggplot(random_plot_data, aes(x = n, y = ov)) +
    geom_line() + geom_smooth(se = F) + geom_hline(yintercept = 0.95) + ggtitle(i) +
    scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.05))
  list_of_plots[[i]] <- plot
}
random_plot_comparison <- ggarrange(plotlist = list_of_plots, ncol = 1)

# cLHS sampling tests ##############################################################################
## Building the models
clhs_rasters <- stack(cropped_bands_2017$B7, cropped_bands_2017$B11, resamp_dem, ndmi) %>%
                mask(area_of_interest)

data <- raster_soil_data %>%
        as("Spatial")
crs(data) <- crs

sample_progression <- c(seq(20, 100, 10), seq(150, 500, 50), seq(600, 1000, 100))
clhs_models <- build_clhs_models(data, target_vars, explanatory_vars,
                                 clhs_rasters, preprocessed_cropped_bands_2017,
                                 n_samples = sample_progression, crs = crs)

## Checking results
clhs_results <- read_csv("tables/clhs_grid_results/model_scores_table.csv")

list_of_plots <- list()
for (i in unique(clhs_results$variable)) {
  clhs_plot_data <- clhs_results %>%
                    filter(variable == i)
  
  plot <- ggplot(clhs_plot_data, aes(x = n_samples, y = R2_valid)) +
    geom_line() + ggtitle(i) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  list_of_plots[[i]] <- plot
}
clhs_r2_plots <- ggarrange(plotlist = list_of_plots, ncol = 1)

## Comparing with ground truth
comparison_plots <- raster_comparisons("clhs")

vars <- str_extract(names(comparison_plots), "(?<=samples_).*") %>%
        unique()
for (i in vars) {
  path <- paste0("figures/clhs_comparison/", i, "_comparisons.png")
  wrap_plots(plotlist = comparison_plots[grepl(i, names(comparison_plots))], ncol = 1)
  ggsave(path, device = "png", bg = "white", units = "mm",
         dpi = 300, width = 300, height = 100 * length(sample_progression),
         limitsize = F)
}

## Comparison summary
clhs_summary <- read_csv("tables/clhs_grid_results/comparison_summary.csv")

list_of_plots <- list()
for (i in unique(clhs_summary$var)) {
  clhs_plot_data <-  clhs_summary %>%
    filter(var == i)
  
  plot <- ggplot(clhs_plot_data, aes(x = n, y = ov)) +
    geom_line() + geom_smooth(se = F) + geom_hline(yintercept = 0.95) + ggtitle(i) +
    scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.05))
  list_of_plots[[i]] <- plot
}
clhs_plot_comparison <- ggarrange(plotlist = list_of_plots, ncol = 1)

# Comparing all sampling methods ###################################################################
p1 <- regular_plot_comparison + plot_annotation(title = "Regular")
p2 <- random_plot_comparison + plot_annotation(title = "Random")
p3 <- clhs_plot_comparison + plot_annotation(title = "cLHS")
wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3)

p1 <- regular_r2_plots + plot_annotation(title = "Regular")
p2 <- random_r2_plots + plot_annotation(title = "Random")
p3 <- clhs_r2_plots + plot_annotation(title = "cLHS")
wrap_elements(p1) + wrap_elements(p2) + wrap_elements(p3)