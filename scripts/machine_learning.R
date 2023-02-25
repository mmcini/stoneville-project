source("scripts/_functions.R")

# Loading the data #################################################################################
raw_soil_data <- read_sf("GIS/shapes/swmru_datasites_point_proj_utm15.shp")
raw_soil_data <- raw_soil_data %>%
                 select(grep("PPM_.", colnames(raw_soil_data)), "OM_PERCENT", "CEC", "WATER_PH")
crs <- CRS("+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs +type=crs")
raw_soil_data <- st_transform(raw_soil_data, crs)

files <- list.files("GIS/rasters/sentinel2", pattern = ".*B\\d{2}.*tiff$")
dates <- as.Date(str_extract(files, "\\d{4}-\\d{2}-\\d{2}"))
bands <- str_extract(files, "B\\d{2}")
raster_files <- tibble(files = paste0("GIS/rasters/sentinel2/", files),
                       dates = dates,
                       bands = bands)

## Assigns bricks to a list based on dates
unique_dates <- unique(raster_files$dates)
list_of_bricks <- list()
for(i in as.character(unique_dates)) {
  paths <- raster_files %>%
           filter(dates == i) %>%
           select(files)
  list_of_bricks[[i]] <- stack(paths[[1]])
  crs(list_of_bricks[[i]]) <- crs
  names(list_of_bricks[[i]]) <- paste0("B", 1:nlayers(list_of_bricks[[i]]))
}

dem <- raster("GIS/rasters/DEM/DEM_1m_small_clip_USGS_3DEP_stoneville.tif")
dem <- projectRaster(dem, crs = crs)

# Extracting data from rasters #####################################################################
area_of_interest <- st_read("GIS/shapes/SWMRU_extent.shp")
area_of_interest <- st_transform(area_of_interest, crs)

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

for (i in names(cor_plot_data[-1])) {
  print(paste(pvalue_plot_data[[i]]))
}

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
target_vars <- c("OM_PERCENT")
explanatory_vars <- c("elevation", "ndmi", "B3", "B8", "B11")
data <- raster_soil_data %>%
        as_tibble()

## Cross validation models and scores
cv_models_list <- build_cv_models(data, target_vars, explanatory_vars, seed = 200)
model_scores_cv_figures(cv_models_list)
model_scores_tables(cv_models_list, "rf_scores", return = T)

## Models using all data available
models_list <- build_models(data, target_vars, explanatory_vars, seed = 200)

# Predicting rasters ###############################################################################
df_cropped_bands_2017 <- cropped_bands_2017 %>%
  as.data.frame(xy = T) %>%
  drop_na()

preproc <- preProcess(df_cropped_bands_2017[-c(1, 2)], method = c("center", "scale"))
preprocessed_cropped_bands_2017 <- predict(preproc, df_cropped_bands_2017) %>%
  rasterFromXYZ(crs = crs)

build_rasters(preprocessed_cropped_bands_2017, models_list, "model_outputs")

# grid sampling tests ##############################################################################
data <- raster_soil_data %>%
        as("Spatial")
crs(data) <- crs

sample_progression <- c(50)
regular_models <- build_sampling_models(data, target_vars, explanatory_vars, area_of_interest,
                                        "regular", preprocessed_cropped_bands_2017,
                                        n_samples = sample_progression, crs = crs)

# random sampling tests ############################################################################
data <- raster_soil_data %>%
        as("Spatial")
crs(data) <- crs

sample_progression <- c(50)
regular_models <- build_sampling_models(data, target_vars, explanatory_vars, area_of_interest,
                                        "random", preprocessed_cropped_bands_2017,
                                        n_samples = sample_progression, crs = crs)

# cLHS sampling tests ##############################################################################
## Building the models
clhs_rasters <- stack(cropped_bands_2017$B3, cropped_bands_2017$B6,
                      cropped_bands_2017$B7, cropped_bands_2017$B8,
                      resamp_dem) %>%
                mask(area_of_interest)

plot(clhs_rasters)

data <- raster_soil_data %>%
        as("Spatial")
crs(data) <- crs

sample_progression <- c(seq(100, 2100, 400), 2145)
clhs_models <- build_clhs_models(data, target_vars, explanatory_vars,
                                 clhs_rasters, preprocessed_cropped_bands_2017,
                                 n_samples = sample_progression, crs = crs)

## Checking results
clhs_results <- read_csv("tables/clhs_tests/model_scores_table.csv")

list_of_plots <- list()
for (i in unique(clhs_results$variable)) {
  clhs_plot_data <- clhs_results %>%
                    filter(variable == i)
  
  plot <- ggplot(clhs_plot_data, aes(x = n_samples, y = R2_valid)) +
    geom_line() + ggtitle(i)
  list_of_plots[[i]] <- plot
}
ggarrange(plotlist = list_of_plots)
