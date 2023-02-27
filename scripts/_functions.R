library(conflicted)
library(ggspatial)
library(tidyverse)
library(extrafont)
library(patchwork)
library(yardstick)
library(corrplot)
library(ggrepel)
library(ggpubr)
library(raster)
library(terra)
library(caret)
library(rgdal)
library(gstat)
library(tune)
library(epiR)
library(clhs)
library(sf)
library(sp)

conflicts_prefer(dplyr::select(), dplyr::filter())

theme_set(theme_bw() + theme(text = element_text(family = "Times New Roman")))

## Functions #######################################################################################
assign_significance <- function(x) {
  case_when(
    x > 0.05 ~ "ns",
    x < 0.05 & x > 0.01 ~ "*",
    x < 0.01 & x > 0.001 ~ "**",
    x < 0.001 ~ "***"
  )
}

ndmi_sentinel_2 <- function(raster, band_8 = "B8", band_11 = "B11") {
  return(calc(raster, function (x) {(x[band_8] - x[band_11]) / (x[band_8] + x[band_11])}))
}

ndwi_sentinel_2 <- function(raster, band_8 = "B8", band_12 = "B12") {
  return(calc(raster, function (x) {(x[band_8] - x[band_12]) / (x[band_8] + x[band_12])}))
}

annotate_valid_scores <- function(data, r2, rmse, y_coord) {
  ## data must have pred and obs columns for relative positioning
  annotate("text", label = c(r2, rmse), x = min(data$obs),
           y = c(y_coord, y_coord * 0.90),
           vjust = 0, hjust = 0, family = "Times New Roman")
}

annotate_coord <- function(y_limit, x_limit, coord_scale) {
  # gets the the highest value from x and y ranges to use as annotates y coord
  y_coord <- min(y_limit) + (max(y_limit) - min(y_limit)) * coord_scale
  x_coord <- min(x_limit) + (max(x_limit) - min(x_limit)) * coord_scale
  coords <- c(y_coord, x_coord)
  return(max(coords))
}

mean_from_folds <- function(data, type = "RMSE", folds_column = "Resample") {
  folds <- unique(data[[folds_column]])
  values <- c()
  count <- 1
  if (type == "RMSE") {
    for (fold in folds) {
      set <- data %>%
        filter(Resample == fold)
      values[count] <- caret::RMSE(pred = set$pred, obs = set$obs)
      count <- count + 1
    }
  } else if (type == "R2") {
    for (fold in folds) {
      set <- data %>%
        filter(Resample == fold)
      values[count] <- caret::R2(pred = set$pred, obs = set$obs)
      count <- count + 1
    }
  } else { return(print("Type must be RMSE or R2")) }
  return(mean(values))
}

validation_plot <- function(cv_set, valid_set, variable = "", dataset = "",
                            model = "", group_by = NULL, coord_scale = 0.9) {
  
  prediction_plot_layout <- list(geom_abline(slope = 1), coord_obs_pred(), theme_bw(),
                                 theme(text = element_text(family = "Times New Roman"),
                                       legend.title = element_blank()))
  
  if (!is.null(dev.list())) {dev.off()} # refreshes device if not null
  data <- list(cv_set, valid_set)
  has_group <- !is.null(group_by)
  plots <- list()
  count <- 1
  for (set in data) {
    if (count == 1) { # calculates the mean of each fold for cv dataset (count = 1)
      rmse <- mean_from_folds(set, type = "RMSE")
      r2 <- mean_from_folds(set, type = "R2")
    } else {
      rmse <- RMSE(pred = set$pred, obs = set$obs)
      r2 <- caret::R2(pred = set$pred, obs = set$obs)
    }
    rmse_text <- paste("RMSE: ", round(rmse, 2))
    r2_text <- paste("R2: ", round(r2, 2))
    if (has_group) { # adds colors to points if grouped
      set$country <- droplevels(set$country)
      colors <- scale_color_manual(drop = T, limits = levels(set$country), values = country_colors)
      colored_points <- list(geom_point(aes_string(color = group_by)), colors)
    }
    plots[[count]] <- ggplot(set, aes(x = obs, y = pred)) +
      {if (has_group) {colored_points}
        else {geom_point(color = "gray")}} +
      xlab(paste("Observed", variable)) +
      ylab(paste("Predicted", variable)) +
      geom_smooth(aes(color = null), method = "lm",se = F,
                  linetype = "dashed", col = "black") +
      prediction_plot_layout
    coord <- annotate_coord(layer_scales(plots[[1]])$y$get_limits(),
                            layer_scales(plots[[1]])$x$get_limits(),
                            coord_scale)
    plots[[count]] <- plots[[count]] +
      annotate_valid_scores(set, r2_text, rmse_text, coord)
    count <- count + 1
  }
  plots[[1]] <- plots[[1]] + ggtitle(paste(dataset, "-", model))
  plots[[2]] <- plots[[2]] + ggtitle(paste(dataset, "-", model))
  return(plots)
}

importance_plot_layout <- list(geom_bar(stat = "identity", width = 0.2),coord_flip(),
                               xlab("Variables"), ylab("Importance (%)"),
                               theme_bw(), theme(text = element_text("Times New Roman")))

build_cv_models <- function(data, target_vars, explan_vars, seed = 200) {
  models_list <- list()
  data_list <- list()
  bar <- txtProgressBar(max = length(target_vars), style = 3)
  
  for (i in seq_len(length(target_vars))){
    explan_data <- data %>%
      select(all_of(explan_vars))
    target_data <- data %>%
      select(target_vars[i])
    
    ## Preprocesing data
    preproc <- preProcess(explan_data, method = c("center", "scale"))
    preprocessed_vars <- predict(preproc, explan_data)
    
    modeling_data <- preprocessed_vars %>%
      add_column("{target_vars[i]}" := target_data[[target_vars[i]]]) %>%
      add_column(geometry = data$geometry)
    
    ## Data partition
    set.seed(seed)
    partition_index <- createDataPartition(modeling_data[[target_vars[i]]], p = 0.8, list = F)
    train_data <- modeling_data %>%
      slice(partition_index)
    valid_data <- modeling_data %>%
      slice(-partition_index)
    data_list[["train"]] <- train_data; data_list[["valid"]] <- valid_data
    
    ## Training
    control <- trainControl(method = "cv", number = 10, savePredictions = T, search = "random")
    param_grid <- expand.grid(mtry = seq(2, 20, 2))
    variables <- paste(explan_vars, collapse = " + ")
    formula <- as.formula(paste(target_vars[i], "~", variables))
    set.seed(seed)
    rf_model <- train(formula, data = train_data, method = "rf",
                      trControl = control, tuneGrid = param_grid)
    
    models_list[[target_vars[i]]] <- list(model = rf_model, data = data_list)
    
    setTxtProgressBar(bar, i)
  }
  close(bar)
  return(models_list)
}

build_models <- function(data, target_vars, explan_vars, seed = 200) {
  models_list <- list()
  data_list <- list()
  bar <- txtProgressBar(max = length(target_vars), style = 3)
  
  for (i in seq_len(length(target_vars))){
    explan_data <- data %>%
      select(all_of(explan_vars))
    target_data <- data %>%
      select(target_vars[i])
    
    ## Preprocesing data
    preproc <- preProcess(explan_data, method = c("center", "scale"))
    preprocessed_vars <- predict(preproc, explan_data)
    
    modeling_data <- preprocessed_vars %>%
      add_column("{target_vars[i]}" := target_data[[target_vars[i]]]) %>%
      add_column(geometry = data$geometry)
    
    ## Training
    control <- trainControl(method = "cv", number = 10, savePredictions = T, search = "random")
    param_grid <- expand.grid(mtry = seq(2, 20, 2))
    variables <- paste(explan_vars, collapse = " + ")
    formula <- as.formula(paste(target_vars[i], "~", variables))
    set.seed(seed)
    rf_model <- train(formula, data = modeling_data, method = "rf",
                      trControl = control, tuneGrid = param_grid)
    models_list[[target_vars[i]]] <- rf_model
    
    setTxtProgressBar(bar, i)
  }
  close(bar)
  return(models_list)
}

model_scores_cv_figures <- function(models) {
  var_names <- names(models)
  for (i in seq_len(length(models))) {
    
    ## Extracting model info
    model <- models[[i]]$model
    train_data <- models[[i]]$data$train
    valid_data <- models[[i]]$data$valid
    var <- var_names[i]
    
    ## Scores
    mtry_param <- model$bestTune$mtry
    model_cv <- model$pred %>%
      filter(mtry == mtry_param) %>% # best model
      arrange(rowIndex)
    model_valid <- predict(model, newdata = valid_data) %>%
      as_tibble() %>%
      add_column(obs = valid_data[[var]]) %>%
      rename(pred = value)
    
    path <- paste0("figures/rf_outputs/", var, "_scores.png")
    model_plots <- validation_plot(model_cv, model_valid, variable = var,
                                   dataset = "Sentinel-2", model = "RF")
    ggarrange(plotlist = model_plots, ncol = 2, common.legend = T, legend = "bottom")
    ggsave(path, device = "png", bg = "white", units = "mm", dpi = 300, width = 200, height = 150)
    
    ## Residuals
    model_resid <- model_cv %>%
      mutate(resid = pred - obs)
    
    path <- paste0("figures/rf_outputs/", var, "_residuals.png")
    ggplot(model_resid, aes(y = resid, x = obs)) +
      ylab("Residuals") + xlab(paste("Observed", var)) +
      geom_point(color = "gray") +
      geom_abline(slope = 0)
    ggsave(path, device = "png", bg = "white", units = "mm", dpi = 300, width = 200, height = 150)
    
    ## Importance
    model_importance <- varImp(model)$importance %>%
      rownames_to_column("variables") %>%
      arrange(Overall) %>%
      mutate(variables = factor(variables,levels = variables))
    
    path <- paste0("figures/rf_outputs/", var, "_importance.png")
    ggplot(model_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
      ggtitle(paste(var, "RF Variable Importance"))
    ggsave(path, device = "png", bg = "white", units = "mm", dpi = 300, width = 200, height = 150)
  }
}

model_scores_tables <- function(models, folder, return = F) {
  var_names <- names(models)
  model_scores <- tibble(variable = character(), n_samples = numeric(),
                         n_cv = numeric(), n_valid = numeric(),
                         rpiq_cv = numeric(), rpiq_valid = numeric(),
                         rpd_cv = numeric(), rpd_valid = numeric(),
                         RMSE_cv = numeric(), RMSE_valid = numeric(),
                         R2_cv = numeric(), R2_valid = numeric())
  
  for (i in seq_len(length(models))) {
    ## Extracting model info
    model <- models[[i]]$model
    train_data <- models[[i]]$data$train
    valid_data <- models[[i]]$data$valid
    var <- var_names[i]
    
    ## Scores
    mtry_param <- model$bestTune$mtry
    model_cv <- model$pred %>%
      filter(mtry == mtry_param) %>% # best model
      arrange(rowIndex)
    model_valid <- predict(model, newdata = valid_data) %>%
      as_tibble() %>%
      add_column(obs = valid_data[[var]]) %>%
      rename(pred = value)
    
    model_scores <- model_scores %>%
      add_row(variable = var, n_samples = nrow(model_cv) + nrow(model_valid),
              n_cv = nrow(model_cv), n_valid = nrow(model_valid),
              rpiq_cv = rpiq_vec(model_cv$obs, model_cv$pred),
              rpiq_valid = rpiq_vec(model_valid$obs, model_valid$pred),
              rpd_cv = rpd_vec(model_cv$obs, model_cv$pred),
              rpd_valid = rpd_vec(model_valid$obs, model_valid$pred),
              RMSE_cv = mean_from_folds(model_cv, type = "RMSE"),
              R2_cv = mean_from_folds(model_cv, type = "R2"),
              RMSE_valid = caret::RMSE(model_valid$pred, model_valid$obs),
              R2_valid = caret::R2(model_valid$pred, model_valid$obs))
  }
  
  ## controls if data should be written or return as an object
  if (!return) {
    path <- paste0("tables/", folder, "/", "model_scores_table.csv")
    write_csv(model_scores, path)
  } else {
    return(model_scores)
  }
}

build_rasters <- function(rasters, models, folder, file_suffix = "_rf_map.tif") {
  var_names <- names(models)
  for (i in seq_len(length(models))) {
    var <- var_names[i]
    path <- paste0("GIS/", folder, "/", var, file_suffix)
    model <- models[[i]]
    
    map <- raster::predict(rasters, model)
    writeRaster(map, path, overwrite= T)
  }
}

build_sampling_rasters <- function(rasters, models, folder, file_suffix = "_rf_map.tif") {
  var_names <- names(models)
  for (i in seq_len(length(models))) {
    var <- var_names[i]
    path <- paste0("GIS/", folder, "/", var, file_suffix)
    model <- models[[i]]$model
    
    map <- raster::predict(rasters, model)
    writeRaster(map, path, overwrite= T)
  }
}

closest_points <- function(known_points, sampled_points) {
  ## returns the known points closest to those sampled by clhs
  dist <- spDists(sampled_points, known_points) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(min = which.min(c_across(everything())))
  return(known_points[dist$min, ])
}

build_clhs_models <- function(data, target_vars, explan_vars,
                              clhs_input, predict_rasters, n_samples = c(), crs, seed = 200) {
  set.seed(seed)
  seeds <- runif(length(n_samples))
  list_of_models <- list()
  model_scores <- tibble(variable = character(), n_samples = numeric(),
                       n_cv = numeric(), n_valid = numeric(),
                       rpiq_cv = numeric(), rpiq_valid = numeric(),
                       rpd_cv = numeric(), rpd_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric())
  
  for (i in seq_len(length(n_samples))) {
    iteration <- paste0("n_samples_", n_samples[i])
    set.seed(seeds[i])
    clhs_samples <- clhs(clhs_input, size = n_samples[i], simple = F)$sampled_data
    crs(clhs_samples) <- crs
    points <- closest_points(data, clhs_samples)
    points_table <- points %>%
                    as_tibble()
    
    print(paste("Iteration", i, "of", length(n_samples)))
    cv_models_list <- build_cv_models(points_table, target_vars, explan_vars, seed = 200)
    
    file_suffix <- paste0(iteration, ".tif")
    build_sampling_rasters(predict_rasters, cv_models_list, "clhs_raster_model_outputs", file_suffix)
    
    model_scores <- model_scores %>%
                    rbind(model_scores_tables(cv_models_list, "clhs_points", return = T))
    list_of_models[[iteration]] <- list(models_list = cv_models_list, points = points)
  }
  
  n_samples_names <- names(list_of_models)
  for (i in n_samples_names) {
    path <- paste0("GIS/clhs_points/", i, ".shp")
    points_sf <- st_as_sf(list_of_models[[i]]$points)
    st_write(points_sf, path, delete_dsn = T)
    
  }
  path <- paste0("tables/clhs_tests/model_scores_table.csv")
  write_csv(model_scores, path)
  return(list_of_models)
}

build_sampling_models <- function(data, target_vars, explan_vars, area, type, 
                                  predict_rasters, n_samples = c(), crs, seed = 200) {
  set.seed(seed)
  seeds <- runif(length(n_samples))
  list_of_models <- list()
  model_scores <- tibble(variable = character(), n_samples = numeric(),
                       n_cv = numeric(), n_valid = numeric(),
                       rpiq_cv = numeric(), rpiq_valid = numeric(),
                       rpd_cv = numeric(), rpd_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric())
  
  for (i in seq_len(length(n_samples))) {
    iteration <- paste0("n_samples_", n_samples[i])
    set.seed(seeds[i])
    samples <- area %>% 
            as("Spatial") %>%
            spsample(n = n_samples[i], type = type)
    points <- closest_points(data, samples)
    points_table <- points %>%
                    as_tibble()
    
    print(paste("Iteration", i, "of", length(n_samples)))
    cv_models_list <- build_cv_models(points_table, target_vars, explan_vars, seed = 200)
    
    file_suffix <- paste0("_", iteration, ".tif")
    build_sampling_rasters(predict_rasters, cv_models_list, paste0(type, "_raster_model_outputs"), file_suffix)
    
    model_scores <- model_scores %>%
                    rbind(model_scores_tables(cv_models_list, paste0(type, "_points"), return = T))
    list_of_models[[iteration]] <- list(models_list = cv_models_list, points = points)
  }
  
  n_samples_names <- names(list_of_models)
  for (i in n_samples_names) {
    path <- paste0("GIS/", type, "_points/", i, ".shp")
    points_sf <- st_as_sf(list_of_models[[i]]$points)
    st_write(points_sf, path, delete_dsn = T)
    
  }
  path <- paste0("tables/", type, "_grid_results/model_scores_table.csv")
  write_csv(model_scores, path)
  return(list_of_models)
}

raster_obs_pred <- function(obs, pred) {
  ## Convert rasters into a table with obs and pred columns
  obs <- obs %>%
         as.data.frame(xy = T)
  names(obs) <- c("x", "y", "obs")
  pred <- pred %>%
          as.data.frame()
  names(pred) <- "pred"
  return(cbind(obs, pred))
}

raster_comparisons <- function(type) {
  files <- list.files("GIS/regular_raster_model_outputs/")
  sample_numbers <- str_extract(files, "\\d+") %>%
                    as.numeric() %>%
                    unique()
  vars <- str_extract(files, ".+?(?=_n)") %>%
          unique()
  
  list_of_plots <- list()
  for (i in seq_len(length(sample_numbers))) {
    for (j in seq_len(length(vars))) {
      path_obs <- paste0("GIS/rf_model_outputs/", vars[j], "_rf_map.tif")
      path_pred <- paste0("GIS/", type, "_raster_model_outputs/",
                          vars[j], "_n_samples_", sample_numbers[i], ".tif")
      obs <- rast(path_obs)
      pred <- rast(path_pred)
      obs_pred_data <- raster_obs_pred(obs, pred)
      summary_stats <-tibble(R2 = R2(obs_pred_data$pred, obs_pred_data$obs),
                             RMSE = RMSE(obs_pred_data$pred, obs_pred_data$obs),
                             obs_iqr = IQR(obs_pred_data$obs),
                             pred_iqr = IQR(obs_pred_data$pred),
                             obs_var = var(obs_pred_data$obs),
                             pred_var = var(obs_pred_data$pred))
      annotation <- c(
        paste("R2:", round(summary_stats$R2, 2)),
        paste("RMSE:", round(summary_stats$RMSE, 2)),
        paste("Ref. IQR:", round(summary_stats$obs_iqr, 2)),
        paste("Pred. IQR:", round(summary_stats$pred_iqr, 2)),
        paste("Ref. variance:", round(summary_stats$obs_var, 2)),
        paste("Pred. variance:", round(summary_stats$pred_var, 2))
      )
      pivoted_raster_data <- raster_obs_pred(obs, pred) %>%
                             pivot_longer(c("obs", "pred"), names_to = "variables", values_to = "values")
      rasters_plot <- ggplot() +
        geom_sf(data = area_of_interest) + xlab("") + ylab("") +
        geom_raster(data = pivoted_raster_data, aes(x = x, y = y, fill = values)) +
        facet_wrap(. ~ variables,
                   labeller = labeller(variables = c(obs = "Reference - 2145 samples",
                                                     pred = paste("Predicted -", sample_numbers[i], "samples")))) +
        scale_fill_distiller(type = "div", palette = "Spectral") +
        guides(fill = guide_colorbar(title = "")) +
        annotation_scale(data = data.frame(variables = c("pred"))) +
        annotation_north_arrow(data = data.frame(variables = c("pred")), ## 'data =' plots in only one facet
                               height = unit(1, "cm"), width = unit(1, "cm"),
                               pad_x = unit(0.1, "cm"), pad_y = unit(0.5, "cm"),
                               which_north = T,
                               style = north_arrow_fancy_orienteering()) +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
      
      density_plot <- ggplot(pivoted_raster_data) +
        ylab("") + xlab("") + ggtitle(paste(vars[j], "-", sample_numbers[i], "samples")) +
        geom_density(aes(x = values, fill = variables), alpha = 0.5) +
        scale_fill_discrete(labels = c("Reference", "Predicted")) +
        theme(legend.title = element_blank())
      
      density_values <- pivoted_raster_data$values[which(pivoted_raster_data$variables == "obs")] %>%
                        density()
      max_y <- max(density_values$y) * 1.2 ## getting max and min y values for annotation
      min_y <- min(density_values$y)
      range <- (max_y - min_y) * 0.8
      density_plot <- density_plot +
        annotate("text", x = max(pivoted_raster_data$values) * 0.8,
                 y = seq(max_y, max_y - range / 2, -range / 10),
                 label = annotation, family= "Times New Roman")
      
      if (i == length(sample_numbers)) {
        rasters_plot <- rasters_plot + theme(legend.position = "bottom")
        density_plot <- density_plot + theme(legend.position = "bottom")
      } else {
        rasters_plot <- rasters_plot + theme(legend.position = "none")
        density_plot <- density_plot + theme(legend.position = "none")
      }
      
      patch <- density_plot + rasters_plot + plot_layout()
      item_name <- paste0(sample_numbers[i], "_samples_", vars[j])
      list_of_plots[[item_name]] <- patch
    }
  }
  return(list_of_plots)
}

## Performs Z-transform back and forth
z_backtrans_info <- function(x) {
  return(list(mean = mean(x), sd = sd(x)))
}

z_trans <- function(x) {
  return((x - mean(x)) / sd(x))
}

z_backtrans <- function(x, backtransform_info) {
  return((x * backtransform_info[["sd"]] + backtransform_info[["mean"]]))
}