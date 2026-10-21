# =============================================================================
# Script: 04_national_ml_attribution_grouped_cv.R
# Project: Climate–Power Outage Relationships Across the United States
# Author: Danlu Cai
#
# Description:
# Evaluate XGBoost and Random Forest associations with county-level power
# outage metrics at daily, monthly and annual resolutions (2015–2022).
# Compare unweighted degree days, weighted degree days and their combination
# using 5-fold county-grouped cross-validation. Export model performance,
# permutation/Gain importance, annual diagnostics and figures.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(xgboost)
  library(ranger)
  library(patchwork)
  library(tibble)
})

# ---- USER CONFIGURATION AND PATH CHECKS --------------------------------------
INPUT_DIR <- Sys.getenv(
  "MERGED_DATA_INPUT_DIR",
  unset = "data/processed/merged"
)

OUTPUT_DIR <- Sys.getenv(
  "ML_ATTRIBUTION_OUTPUT_DIR",
  unset = "results/national_ml_attribution"
)

CSV_DIR  <- file.path(OUTPUT_DIR, "csv_results")
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
dir.create(CSV_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(INPUT_DIR)) {
  stop("Input directory not found: ", INPUT_DIR)
}

SEED <- 42L
N_FOLDS <- 5L
N_TREES_RF <- 200L
N_ROUNDS_XGB <- 100L

CATEGORY_COLORS <- c(
  "Temperature & HVAC Load" = "#ff7f00", 
  "Storm"                   = "#e31a1c",
  "Temporal"                = "#31a354", 
  "Exposure & Geography"   = "#636363",
  "Other"                   = "#999999"
)

configs <- list(
  Daily = list(
    file = "Digicare_Merged_Daily_2015_2022_CONUS_spatial.csv",
    metrics = c("daily_coh", "total_outage_hours"),
    unweighted = c("HDD", "CDD"), 
    weighted = c("weighted_CDD", "weighted_HDD"),
    base = c("month", "heat_ext", "cold_ext", "T2_C", "event_count",
             "max_magnitude", "event_types", "state")
  ),
  Monthly = list(
    file = "Digicare_Merged_Monthly_2015_2022_CONUS_spatial.csv",
    metrics = c("monthly_coh", "monthly_outage_hours"),
    unweighted = c("monthly_CDD", "monthly_HDD"),
    weighted = c("weighted_CDD", "weighted_HDD"),
    base = c("month", "heat_days", "cold_days", "mean_T2C",
             "monthly_event_count", "monthly_max_mag", "storm_days", "state")
  ),
  Yearly = list(
    file = "Digicare_Merged_Yearly_2015_2022_CONUS_spatial.csv",
    metrics = c("annual_coh", "annual_outage_hours"),
    unweighted = c("annual_CDD", "annual_HDD"),
    weighted = c("weighted_CDD", "weighted_HDD"),
    base = c("avg_T2C", "annual_event_count", "annual_max_mag",
             "annual_heat_days", "annual_cold_days", "annual_storm_days", "state")
  )
)

category_for <- function(features) {
  dplyr::case_when(
    grepl("CDD|HDD|T2|heat|cold", features, ignore.case = TRUE) ~ "Temperature & HVAC Load",
    grepl("storm|event|magnitude", features, ignore.case = TRUE) ~ "Storm",
    grepl("^month$|^year$", features, ignore.case = TRUE) ~ "Temporal",
    grepl("state|population|density", features, ignore.case = TRUE) ~ "Exposure & Geography",
    TRUE ~ "Other"
  )
}

prepare_model_data <- function(df, metric, predictors) {
  county_col <- intersect(c("fips", "fips_code", "county", "GEOID"), names(df))[1]
  if (is.na(county_col)) stop("No county identifier found")
  if (!metric %in% names(df)) stop("Missing outcome: ", metric)
  
  predictors <- intersect(predictors, names(df))
  predictors <- setdiff(predictors, c(county_col, metric, "county_area", "County_Area", "area", "Area"))
  if (!length(predictors)) return(NULL)
  
  x <- df[, predictors, drop = FALSE]
  
  categorical <- vapply(x, function(v) is.character(v) || is.factor(v), logical(1))
  categorical[names(categorical) %in% c("state", "State", "event_types")] <- TRUE
  
  for (nm in names(x)) {
    if (categorical[[nm]]) {
      x[[nm]] <- factor(as.character(x[[nm]]))
    } else {
      x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
    }
  }
  
  target_raw <- suppressWarnings(as.numeric(df[[metric]]))
  groups <- as.character(df[[county_col]])
  valid <- !is.na(groups) & nzchar(groups) & is.finite(target_raw) & target_raw >= 0 &
    stats::complete.cases(x)
  
  if (sum(valid) < 50L) return(NULL)
  
  x <- x[valid, , drop = FALSE]
  
  keep_cols <- vapply(x, function(v) length(unique(v)) > 1L, logical(1))
  x <- x[, keep_cols, drop = FALSE]
  if (!ncol(x)) return(NULL)
  
  predictors <- names(x)
  target <- log10(target_raw[valid] + 1)
  groups <- groups[valid]
  
  if (length(unique(groups)) < N_FOLDS || stats::sd(target) == 0) return(NULL)
  
  terms_formula <- stats::reformulate(sprintf("`%s`", gsub("`", "", names(x), fixed = TRUE)))
  mm <- stats::model.matrix(terms_formula, data = x)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (!ncol(mm)) return(NULL)
  
  assign_idx <- attr(stats::model.matrix(terms_formula, data = x), "assign")
  term_labels <- attr(stats::terms(terms_formula), "term.labels")
  original <- gsub("^`|`$", "", term_labels)
  feature_source <- original[assign_idx[colnames(stats::model.matrix(terms_formula, data = x)) != "(Intercept)"]]
  names(feature_source) <- colnames(mm)
  
  list(x = mm, y = target, county = groups, source = feature_source,
       original_predictors = predictors, n = length(target))
}

run_cv <- function(prepared, model_type, seed = SEED) {
  set.seed(seed)
  counties <- sample(unique(prepared$county))
  county_fold <- setNames(rep(seq_len(N_FOLDS), length.out = length(counties)), counties)
  folds <- unname(county_fold[prepared$county])
  
  x <- prepared$x
  y <- prepared$y
  feature_names <- colnames(x)
  fold_stats <- vector("list", N_FOLDS)
  importance <- matrix(0, nrow = N_FOLDS, ncol = length(feature_names),
                       dimnames = list(NULL, feature_names))
  
  for (k in seq_len(N_FOLDS)) {
    tr <- which(folds != k)
    te <- which(folds == k)
    if (!length(tr) || !length(te)) stop("Empty cross-validation fold")
    
    if (model_type == "XGBoost") {
      model <- xgboost::xgboost(
        data = x[tr, , drop = FALSE], label = y[tr],
        nrounds = N_ROUNDS_XGB, objective = "reg:squarederror",
        eta = 0.1, max_depth = 6, nthread = 1, verbosity = 0
      )
      pred <- as.numeric(stats::predict(model, x[te, , drop = FALSE]))
      imp <- xgboost::xgb.importance(model = model)
      if (nrow(imp)) {
        matched <- intersect(imp$Feature, feature_names)
        importance[k, matched] <- imp$Gain[match(matched, imp$Feature)]
      }
    } else {
      # Use matrix interface (x, y) to bypass ranger's formula parser
      model <- ranger::ranger(
        x = x[tr, , drop = FALSE],
        y = y[tr],
        num.trees = N_TREES_RF,
        importance = "permutation",
        seed = seed + k
      )
      pred <- as.numeric(predict(model, data = x[te, , drop = FALSE])$predictions)
      imp <- model$variable.importance
      importance[k, intersect(names(imp), feature_names)] <- imp[intersect(names(imp), feature_names)]
    }
    
    tss <- sum((y[te] - mean(y[te]))^2)
    r2 <- if (tss > 0) 1 - sum((y[te] - pred)^2) / tss else NA_real_
    fold_stats[[k]] <- tibble::tibble(
      Fold = k, 
      Test_Counties = length(unique(prepared$county[te])),
      Test_Rows = length(te), 
      R2 = r2,
      MAE = mean(abs(y[te] - pred))
    )
  }
  
  fold_stats <- dplyr::bind_rows(fold_stats)
  
  source <- prepared$source[colnames(importance)]
  raw_imp <- vapply(unique(source), function(s) rowSums(importance[, source == s, drop = FALSE]),
                    numeric(N_FOLDS))
  if (is.null(dim(raw_imp))) raw_imp <- matrix(raw_imp, ncol = 1,
                                               dimnames = list(NULL, unique(source)))
  
  signed_mean <- colMeans(raw_imp)
  display_imp <- pmax(raw_imp, 0)
  denom <- rowSums(display_imp)
  share <- sweep(display_imp, 1, ifelse(denom > 0, denom, 1), "/") * 100
  
  imp_df <- tibble::tibble(
    Feature = colnames(raw_imp),
    Importance = colMeans(share),
    SD = apply(share, 2, stats::sd),
    Raw_Mean = signed_mean,
    Category = category_for(colnames(raw_imp))
  ) %>% dplyr::arrange(dplyr::desc(Importance))
  
  list(performance = fold_stats, importance = imp_df,
       r2 = mean(fold_stats$R2, na.rm = TRUE),
       mae = mean(fold_stats$MAE, na.rm = TRUE))
}

get_year <- function(df) {
  if ("year" %in% names(df)) return(suppressWarnings(as.integer(df$year)))
  if ("Year" %in% names(df)) return(suppressWarnings(as.integer(df$Year)))
  if ("date" %in% names(df)) return(as.integer(format(as.Date(df$date), "%Y")))
  if ("month_stamp" %in% names(df)) return(as.integer(format(as.Date(df$month_stamp), "%Y")))
  stop("No year identifier found")
}

plot_importance <- function(imp, subtitle) {
  top <- head(imp, 8)
  ggplot2::ggplot(top, ggplot2::aes(x = stats::reorder(Feature, Importance),
                                    y = Importance, fill = Category)) +
    ggplot2::geom_col(alpha = 0.75) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, Importance - SD),
                                        ymax = Importance + SD), width = 0.25) +
    ggplot2::coord_flip() + 
    ggplot2::scale_fill_manual(values = CATEGORY_COLORS) +
    ggplot2::labs(x = NULL, y = "Mean importance share (%)", subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_yearly <- function(yearly, category) {
  if (!nrow(yearly)) return(NULL)
  
  a <- ggplot2::ggplot(yearly, ggplot2::aes(year, R2)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::labs(x = NULL, y = "County-grouped CV R²", title = "Year-specific performance") +
    ggplot2::theme_minimal(base_size = 12)
  
  if (!nrow(category)) return(a)
  
  b <- ggplot2::ggplot(category, ggplot2::aes(x = factor(year), y = Importance, fill = Category)) +
    ggplot2::geom_col() + 
    ggplot2::scale_fill_manual(values = CATEGORY_COLORS) +
    ggplot2::labs(x = "Year", y = "Importance share (%)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  
  a / b
}

# ---- MAIN EXECUTION LOOP -----------------------------------------------------
master <- list()

for (scale_name in names(configs)) {
  cfg <- configs[[scale_name]]
  path <- file.path(INPUT_DIR, cfg$file)
  if (!file.exists(path)) {
    warning("Skipping missing dataset: ", path)
    next
  }
  
  message("Processing ", scale_name, ": ", path)
  df <- readr::read_csv(path, show_col_types = FALSE, guess_max = 10000)
  df$Analysis_Year <- get_year(df)
  
  density <- intersect(c("pop_density", "Pop_Density", "density"), names(df))[1]
  population <- intersect(c("population", "Pop", "pop"), names(df))[1]
  exposure <- if (!is.na(density)) density else if (!is.na(population)) population else character()
  
  runs <- list(
    Run1_Unweighted = unique(c(cfg$base, cfg$unweighted, exposure)),
    Run2_Weighted   = unique(c(cfg$base, cfg$weighted)),
    Run3_Combined   = unique(c(cfg$base, cfg$unweighted, cfg$weighted, exposure))
  )
  
  for (metric in cfg$metrics) {
    if (!metric %in% names(df)) {
      warning("Missing outcome: ", metric)
      next
    }
    
    for (run_name in names(runs)) {
      predictors <- intersect(runs[[run_name]], names(df))
      if (!length(predictors)) next
      
      full <- prepare_model_data(df, metric, predictors)
      if (is.null(full)) {
        warning("Insufficient usable observations: ", scale_name, " ", metric, " ", run_name)
        next
      }
      
      for (model_type in c("XGBoost", "RF")) {
        prefix <- paste(scale_name, metric, model_type, run_name, sep = "_")
        message("  ", prefix)
        
        overall <- run_cv(full, model_type)
        readr::write_csv(overall$importance, file.path(CSV_DIR, paste0("FeatureImportance_", prefix, ".csv")))
        readr::write_csv(overall$performance, file.path(CSV_DIR, paste0("CV_Folds_", prefix, ".csv")))
        
        yearly <- list()
        cats <- list()
        
        for (yr in sort(unique(stats::na.omit(df$Analysis_Year)))) {
          subset <- df[!is.na(df$Analysis_Year) & df$Analysis_Year == yr, , drop = FALSE]
          prepared <- prepare_model_data(subset, metric, predictors)
          if (is.null(prepared)) next
          
          result <- run_cv(prepared, model_type, seed = SEED)
          yearly[[as.character(yr)]] <- tibble::tibble(
            year = yr, R2 = result$r2, MAE = result$mae,
            Valid_Rows = prepared$n, Counties = length(unique(prepared$county))
          )
          cats[[as.character(yr)]] <- result$importance %>%
            dplyr::group_by(Category) %>%
            dplyr::summarise(Importance = sum(Importance), .groups = "drop") %>%
            dplyr::mutate(year = yr)
        }
        
        yearly_df <- dplyr::bind_rows(yearly)
        category_df <- dplyr::bind_rows(cats)
        
        readr::write_csv(yearly_df, file.path(CSV_DIR, paste0("YearlyTrends_", prefix, ".csv")))
        readr::write_csv(category_df, file.path(CSV_DIR, paste0("YearlyCategoryImportance_", prefix, ".csv")))
        
        top_plot <- plot_importance(
          overall$importance,
          sprintf("County-grouped CV R² = %.3f; MAE = %.3f", overall$r2, overall$mae)
        )
        bottom_plot <- plot_yearly(yearly_df, category_df)
        figure <- if (is.null(bottom_plot)) top_plot else top_plot / bottom_plot
        
        ggplot2::ggsave(
          file.path(PLOT_DIR, paste0("Attribution_", prefix, ".png")),
          plot = figure, width = 11, height = if (is.null(bottom_plot)) 6 else 14,
          dpi = 300, bg = "white"
        )
        
        master[[prefix]] <- tibble::tibble(
          Scale = scale_name, Metric = metric,
          Model = model_type, Run = run_name,
          Num_Features = length(predictors),
          Valid_Rows = full$n,
          Overall_R2 = overall$r2,
          Overall_MAE = overall$mae
        )
      }
    }
  }
}

if (length(master)) {
  readr::write_csv(dplyr::bind_rows(master), file.path(CSV_DIR, "Master_ML_Performance_Summary.csv"))
  message("Analysis complete. Results saved to: ", OUTPUT_DIR)
} else {
  warning("No model results generated. Check input filenames, columns and missing values.")
}
