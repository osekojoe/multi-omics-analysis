# ==============================================================================
# Script: fabia_full_script_with_safe_heatmap.R
# Description: Simulate data, fit FABIA, and evaluate score & loading recovery.
#              Includes robust heatmap helper to avoid 'breaks not unique' errors.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup & Helper Functions
# ------------------------------------------------------------------------------

# Null coalescing operator helper
`%||%` <- function(x, y) if (!is.null(x)) x else y

# Package management
required_pkgs <- c("fabia", "pheatmap", "ggplot2", "dplyr", "tidyr", "clue", "SUMO")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))

if (length(missing_pkgs)) {
  message("Please install missing packages: ", paste(missing_pkgs, collapse = ", "))
}

# Load libraries
suppressPackageStartupMessages({
  library(fabia)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(clue)
  library(SUMO)
})

# ------------------------------------------------------------------------------
# 2. Data Simulation
# ------------------------------------------------------------------------------

# Simulation 1
sim1 <- simulateMultiOmics(
  vector_features = c(4000, 2500),
  n_samples = 100,
  n_factors = 3,
  snr = 2.0,
  signal.samples = c(3, 1),
  signal.features = list(c(4.5, 0.5), c(4.5, 0.5)),
  factor_structure = "mixed",
  num.factor = "multiple",
  seed = 123
)

# Simulation 2
sim2 <- simulateMultiOmics(
  vector_features = c(4000, 3500),
  n_samples = 100,
  n_factors = 3,
  snr = 1.0,
  signal.samples = c(3, 1),
  signal.features = list(c(2.5, 0.5), c(3, 2.5)),
  factor_structure = "mixed",
  num.factor = "multiple",
  seed = 123
)

# Combine omics from different simulations into one sim_object
sim_object <- sim1
sim_object$omics$omic2 <- sim2$omics$omic2
sim_object[["list_betas"]][[2]] <- sim2[["list_betas"]][[2]]

# Prepare final data list
simX_list <- sim_object$omics
names(simX_list) <- paste0("omic", seq_along(simX_list))

# ------------------------------------------------------------------------------
# 3. True Data Visualization & Helpers
# ------------------------------------------------------------------------------

# Extract true scores
true_scores <- do.call(cbind, lapply(sim_object$list_alphas, as.numeric))
colnames(true_scores) <- names(sim_object$list_alphas)
rownames(true_scores) <- rownames(sim_object$omics[[1]]) %||% seq_len(nrow(true_scores))

# Visualize True Scores (Optional)
df_scores_true <- as.data.frame(true_scores) %>%
  mutate(Sample = rownames(true_scores)) %>%
  pivot_longer(cols = -Sample, names_to = "Factor", values_to = "Score")

p_true <- ggplot(df_scores_true, aes(x = Sample, y = Score, color = Factor, group = Factor)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~Factor, scales = "free_y") +
  labs(title = "True Factor Scores (Simulated Data)", x = "Sample", y = "Score") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

print(p_true)

# Helper: Extract true loadings per block
get_true_loadings_per_block <- function(sim_object) {
  true_loadings_list <- list()
  for (i in seq_along(sim_object$list_betas)) {
    block_betas <- sim_object$list_betas[[i]]
    block_mats <- lapply(block_betas, function(b) matrix(as.numeric(b), ncol = 1))
    block_combined <- do.call(cbind, block_mats)
    
    # Generate systematic rownames
    rownames(block_combined) <- paste0(
      names(sim_object$omics)[i] %||% paste0("omic", i),
      "_feature_",
      seq_len(nrow(block_combined))
    )
    true_loadings_list[[i]] <- block_combined
  }
  names(true_loadings_list) <- names(sim_object$list_betas) %||% 
    paste0("omic", seq_along(sim_object$list_betas))
  
  return(true_loadings_list)
}

true_loadings <- get_true_loadings_per_block(sim_object)

# ------------------------------------------------------------------------------
# 4. Preprocessing Helpers
# ------------------------------------------------------------------------------

preprocess_block <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 0 || nrow(X) == 0) return(X)
  
  # Remove constant columns
  col_var <- apply(X, 2, var, na.rm = TRUE)
  if (any(col_var == 0, na.rm = TRUE)) X <- X[, col_var != 0, drop = FALSE]
  
  # Remove constant rows
  row_var <- apply(X, 1, var, na.rm = TRUE)
  if (any(row_var == 0, na.rm = TRUE)) X <- X[row_var != 0, , drop = FALSE]
  
  # Scale (center + scale) and impute Inf/NA with 0
  if (nrow(X) > 0 && ncol(X) > 0) {
    X <- tryCatch(scale(X), error = function(e) X)
    X[!is.finite(X)] <- 0
  }
  X
}

ensure_block_colnames <- function(Xlist, block_names = NULL) {
  block_names <- block_names %||% names(Xlist) %||% paste0("omic", seq_along(Xlist))
  for (i in seq_along(Xlist)) {
    if (is.null(colnames(Xlist[[i]]))) {
      colnames(Xlist[[i]]) <- paste0(block_names[i], "_", seq_len(ncol(Xlist[[i]])))
    } else {
      nn <- colnames(Xlist[[i]])
      nn <- make.unique(nn)
      colnames(Xlist[[i]]) <- paste0(block_names[i], "_", nn)
    }
  }
  names(Xlist) <- block_names
  Xlist
}

# ------------------------------------------------------------------------------
# 5. Robust Heatmap Helper
# ------------------------------------------------------------------------------

#' Safe Pheatmap
#' Wraps pheatmap to handle numerical precision issues causing 'breaks are not unique'
safe_pheatmap <- function(mat,
                          display_numbers = TRUE,
                          number_format = function(x) sprintf("%.2f", x),
                          fontsize_number = 10,
                          number_color = "black",
                          main = NULL,
                          ncolors = 50,
                          color_palette = NULL,
                          na.value = 0,
                          ...) {
  
  # Ensure input is matrix
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = 1, ncol = 1)
  mat <- as.matrix(mat)
  
  # Replace NA
  if (all(is.na(mat))) mat[] <- na.value
  mat[is.na(mat)] <- na.value
  
  # Prepare display numbers
  disp_nums <- NULL
  if (isTRUE(display_numbers)) {
    disp_nums <- matrix(number_format(mat), nrow = nrow(mat), ncol = ncol(mat))
    rownames(disp_nums) <- rownames(mat)
    colnames(disp_nums) <- colnames(mat)
  }
  
  # Setup Palette
  if (is.null(color_palette)) {
    color_palette <- colorRampPalette(c("white", "steelblue"))(ncolors)
  } else {
    if (is.function(color_palette)) color_palette <- color_palette(ncolors)
    if (length(color_palette) < ncolors) {
      color_palette <- grDevices::colorRampPalette(color_palette)(ncolors)
    }
  }
  
  # Numeric range with epsilon padding to ensure strict breaks
  rng_min <- min(mat, na.rm = TRUE)
  rng_max <- max(mat, na.rm = TRUE)
  
  if (isTRUE(all.equal(rng_min, rng_max))) {
    eps <- max(1e-8, abs(rng_min) * 1e-8, 1e-8)
    rng_min2 <- rng_min - eps
    rng_max2 <- rng_max + eps
  } else {
    pad <- max(1e-12, (rng_max - rng_min) * 1e-8)
    rng_min2 <- rng_min - pad
    rng_max2 <- rng_max + pad
  }
  
  # Create strictly increasing breaks
  breaks <- seq(from = rng_min2, to = rng_max2, length.out = (ncolors + 1))
  
  # Render
  pheatmap::pheatmap(
    mat,
    color = color_palette,
    breaks = breaks,
    main = main,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = disp_nums,
    number_color = number_color,
    fontsize_number = fontsize_number,
    ...
  )
}

# ------------------------------------------------------------------------------
# 6. Core: Fit FABIA
# ------------------------------------------------------------------------------

fit_fabia <- function(Xlist, p = 3, alpha = 0.05, seed = 123, 
                      center = 2, norm = 1, cyc = 1000, verbose = TRUE) {
  set.seed(seed)
  stopifnot(is.list(Xlist), length(Xlist) >= 1)
  
  # Preprocess
  Xlist_clean <- lapply(Xlist, preprocess_block)
  Xlist_clean <- ensure_block_colnames(Xlist_clean)
  
  # Concatenate features (FABIA expects Features x Samples)
  X_concat <- t(do.call(cbind, Xlist_clean))
  
  if (verbose) {
    message("FABIA input matrix: ", nrow(X_concat), " features x ", ncol(X_concat), " samples")
  }
  
  # Run FABIA
  fab_fit <- fabia::fabia(X_concat, p = p, alpha = alpha, cyc = cyc, center = center, norm = norm)
  
  # Extract results
  factor_scores <- t(fab_fit@Z)  # Samples x Factors
  loadings_all  <- fab_fit@L     # Features x Factors
  
  # Split loadings back into blocks
  block_sizes <- vapply(Xlist_clean, ncol, integer(1))
  cuts <- cumsum(block_sizes)
  
  loadings_per_block <- lapply(seq_along(block_sizes), function(i) {
    start <- if (i == 1) 1 else (cuts[i - 1] + 1)
    end   <- cuts[i]
    loadings_all[start:end, , drop = FALSE]
  })
  names(loadings_per_block) <- names(Xlist_clean)
  
  list(
    fabia_model = fab_fit,
    factor_scores = factor_scores,
    loadings_per_block = loadings_per_block,
    block_sizes = block_sizes
  )
}

# ------------------------------------------------------------------------------
# 7. Evaluation: Scores
# ------------------------------------------------------------------------------

evaluate_scores <- function(true_mat, est_mat, n_factors_eval = NULL, 
                            plots = TRUE, title_prefix = "") {
  
  true_mat <- as.matrix(true_mat)
  est_mat  <- as.matrix(est_mat)
  
  if (nrow(true_mat) != nrow(est_mat)) stop("Sample mismatch between true and estimated factors.")
  
  # Determine dimensions
  K_true <- ncol(true_mat)
  K_est <- ncol(est_mat)
  n_f <- n_factors_eval %||% min(K_true, K_est)
  
  trueK <- true_mat[, seq_len(n_f), drop = FALSE]
  estK  <- est_mat[, seq_len(min(K_est, n_f)), drop = FALSE]
  K <- ncol(trueK)
  
  # 1. Calculate correlations
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  if (is.null(dim(cor_mat))) cor_mat <- matrix(cor_mat, nrow = 1, ncol = 1)
  cor_mat[is.na(cor_mat)] <- 0
  
  # Name rows/cols
  rn <- colnames(trueK) %||% paste0("True_F", seq_len(ncol(trueK)))
  cn <- colnames(estK)  %||% paste0("Est_F", seq_len(ncol(estK)))
  rownames(cor_mat) <- rn; colnames(cor_mat) <- cn
  
  # 2. Match factors (Hungarian Algorithm / LSAP)
  abs_cor <- abs(cor_mat)
  cost <- 1 - abs_cor
  assignment <- clue::solve_LSAP(as.matrix(cost))
  cols <- as.integer(assignment)
  
  # 3. Align signs based on match
  matched_corrs <- numeric(K)
  matched_signs <- numeric(K)
  
  for (i in seq_len(K)) {
    j <- cols[i]
    cval <- cor_mat[i, j]
    if (is.na(cval)) cval <- 0
    
    matched_corrs[i] <- cval
    sgn <- sign(cval)
    if (sgn == 0) sgn <- 1
    matched_signs[i] <- sgn
  }
  
  est_matched_signed <- estK[, cols, drop = FALSE] %*% diag(matched_signs)
  colnames(est_matched_signed) <- paste0("Est_matched_", seq_len(K))
  
  # 4. Metrics
  mean_abs_cor <- mean(abs(matched_corrs), na.rm = TRUE)
  prop_ge_0.7 <- mean(abs(matched_corrs) >= 0.7, na.rm = TRUE)
  prop_ge_0.5 <- mean(abs(matched_corrs) >= 0.5, na.rm = TRUE)
  
  R2 <- numeric(K); MSE <- numeric(K)
  for (i in seq_len(K)) {
    fit <- lm(trueK[, i] ~ est_matched_signed[, i])
    ssr <- sum((predict(fit) - mean(trueK[, i]))^2)
    sst <- sum((trueK[, i] - mean(trueK[, i]))^2)
    R2[i]  <- if (sst == 0) NA else ssr / sst
    MSE[i] <- mean((trueK[, i] - est_matched_signed[, i])^2)
  }
  
  df_summary <- data.frame(
    true_factor = colnames(trueK) %||% paste0("True", seq_len(K)),
    matched_estimated_index = cols,
    pearson_corr = matched_corrs,
    abs_corr = abs(matched_corrs),
    sign = matched_signs,
    R2 = R2,
    MSE = MSE,
    stringsAsFactors = FALSE
  )
  
  # 5. Plotting
  if (plots) {
    tmp <- cor_mat
    if (length(unique(as.vector(tmp))) == 1) tmp <- tmp + 1e-8
    
    safe_pheatmap(tmp,
                  main = paste0(title_prefix, " Score Correlation (true x est)"),
                  display_numbers = TRUE,
                  number_color = "black",
                  fontsize_number = 10)
    
    df_plot <- data.frame(
      true = as.vector(trueK),
      est = as.vector(est_matched_signed),
      factor = factor(rep(seq_len(K), each = nrow(trueK)))
    )
    
    print(
      ggplot(df_plot, aes(x = true, y = est)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~factor, scales = "free") +
        labs(title = paste0(title_prefix, " True vs Estimated Scores"))
    )
  }
  
  list(
    correlation_matrix = cor_mat,
    assignment = cols,
    matched_signs = matched_signs,
    summary = df_summary,
    mean_abs_cor = mean_abs_cor,
    prop_abs_cor_ge_0.7 = prop_ge_0.7,
    prop_abs_cor_ge_0.5 = prop_ge_0.5,
    est_signed_matched = est_matched_signed,
    true_matched = trueK
  )
}

# ------------------------------------------------------------------------------
# 8. Evaluation: Loadings (Per Block)
# ------------------------------------------------------------------------------

evaluate_loadings_per_block <- function(true_loadings_list, est_loadings_list, 
                                        factor_assignment, matched_signs,
                                        plots = TRUE, title_prefix = "") {
  blocks <- names(true_loadings_list)
  res_block <- list()
  
  for (block in blocks) {
    true_mat <- as.matrix(true_loadings_list[[block]])
    est_mat  <- as.matrix(est_loadings_list[[block]])
    
    # Align by rownames if possible
    if (!is.null(rownames(true_mat)) && 
        !is.null(rownames(est_mat)) && 
        all(rownames(true_mat) %in% rownames(est_mat))) {
      est_mat <- est_mat[rownames(true_mat), , drop = FALSE]
    } else if (nrow(true_mat) == nrow(est_mat)) {
      warning("Row names differ but counts match; using index alignment for block: ", block)
    } else {
      stop("Feature mismatch in block ", block)
    }
    
    # Filter assignments to ensure validity
    factor_assignment_use <- factor_assignment
    matched_signs_use <- matched_signs
    max_factor <- ncol(est_mat)
    
    if (max(factor_assignment_use) > max_factor) {
      valid_idx <- which(factor_assignment_use <= max_factor)
      factor_assignment_use <- factor_assignment_use[valid_idx]
      matched_signs_use <- matched_signs_use[valid_idx]
    }
    
    # Match and Sign flip
    est_matched <- est_mat[, factor_assignment_use, drop = FALSE] %*% diag(matched_signs_use)
    colnames(est_matched) <- paste0("Est_matched_", seq_len(ncol(est_matched)))
    
    K_true_block <- ncol(true_mat)
    K_est_block  <- ncol(est_matched)
    K <- min(K_true_block, K_est_block)
    
    if (K == 0) {
      warning("No factors to compare in block ", block)
      res_block[[block]] <- list(summary = data.frame(), est_matched = est_matched)
      next
    }
    
    # Calculate Metrics
    corrs <- numeric(K); rmse <- numeric(K)
    for (k in seq_len(K)) {
      a <- true_mat[, k]
      bvec <- est_matched[, k]
      corrs[k] <- suppressWarnings(cor(a, bvec, use = "pairwise.complete.obs"))
      rmse[k] <- sqrt(mean((a - bvec)^2, na.rm = TRUE))
    }
    
    df <- data.frame(
      block = block,
      factor = paste0("Factor", seq_len(K)),
      pearson_corr = corrs,
      RMSE = rmse,
      stringsAsFactors = FALSE
    )
    
    if (plots) {
      # Correlation heatmap
      cor_block <- cor(true_mat[, seq_len(K), drop = FALSE], 
                       est_matched[, seq_len(K), drop = FALSE], 
                       use = "pairwise.complete.obs")
      
      if (is.null(dim(cor_block))) cor_block <- matrix(cor_block, nrow = 1, ncol = 1)
      cor_block[is.na(cor_block)] <- 0
      
      # Handle constant values for heatmap stability
      if (nrow(cor_block) > 1 && ncol(cor_block) > 1 && 
          length(unique(as.vector(cor_block))) == 1) {
        cor_block <- cor_block + 1e-8
      }
      
      rownames(cor_block) <- paste0("True_F", seq_len(nrow(cor_block)))
      colnames(cor_block) <- paste0("Est_F", seq_len(ncol(cor_block)))
      
      safe_pheatmap(cor_block,
                    main = paste0(title_prefix, " Loadings Correlation (", block, ")"),
                    display_numbers = TRUE,
                    number_color = "black",
                    fontsize_number = 8)
      
      # Scatterplot
      df_plot <- bind_rows(lapply(seq_len(K), function(k) {
        data.frame(
          feature = rownames(true_mat),
          true_loading = true_mat[, k],
          est_loading = est_matched[, k],
          factor = paste0("Factor", k),
          stringsAsFactors = FALSE
        )
      }))
      
      print(
        ggplot(df_plot, aes(x = true_loading, y = est_loading)) +
          geom_point(alpha = 0.4, size = 0.6) +
          geom_smooth(method = "lm", se = FALSE) +
          facet_wrap(~factor, scales = "free") +
          labs(title = paste0(title_prefix, " True vs Est Loadings (", block, ")"))
      )
    }
    
    res_block[[block]] <- list(summary = df, est_matched = est_matched)
  }
  
  res_block
}

# ------------------------------------------------------------------------------
# 9. Master Wrapper: Fit & Evaluate
# ------------------------------------------------------------------------------

fit_and_evaluate_fabia <- function(Xlist, sim_object = NULL, p = 3, alpha = 0.05, 
                                   seed = 123, plots = TRUE, n_factors_eval = NULL, 
                                   verbose = TRUE) {
  
  # 1. Fit Model
  fit_res <- fit_fabia(Xlist, p = p, alpha = alpha, seed = seed, verbose = verbose)
  
  eval_res <- NULL
  
  # 2. Evaluate against Truth (if provided)
  if (!is.null(sim_object) && !is.null(sim_object$list_alphas)) {
    mats <- sim_object$list_alphas
    true_mat <- do.call(cbind, lapply(mats, as.numeric))
    colnames(true_mat) <- names(mats)
    rn <- rownames(sim_object$omics[[1]]) %||% seq_len(nrow(true_mat))
    rownames(true_mat) <- rn
    
    # Evaluate Scores
    eval_res <- evaluate_scores(
      true_mat, 
      fit_res$factor_scores, 
      n_factors_eval = n_factors_eval, 
      plots = plots, 
      title_prefix = "FABIA"
    )
    
    # Evaluate Loadings
    if (!is.null(sim_object$list_betas)) {
      true_loadings_list <- get_true_loadings_per_block(sim_object)
      
      eval_res$loadings_per_block <- evaluate_loadings_per_block(
        true_loadings_list,
        fit_res$loadings_per_block,
        factor_assignment = eval_res$assignment,
        matched_signs = eval_res$matched_signs,
        plots = plots,
        title_prefix = "FABIA"
      )
    }
  }
  
  list(
    fit = fit_res,
    evaluation = eval_res
  )
}

# ------------------------------------------------------------------------------
# 10. Execution
# ------------------------------------------------------------------------------

res <- fit_and_evaluate_fabia(
  Xlist = simX_list, 
  sim_object = sim_object, 
  p = 3, 
  alpha = 0.05, 
  seed = 123, 
  plots = TRUE, 
  n_factors_eval = 3
)

# Standard FABIA plots
par(mfrow = c(2, 2))
plot(res[["fit"]][["factor_scores"]], main = "FABIA scores omic1", cex = 0.3)
plot(res[["fit"]][["loadings_per_block"]][["omic1"]], main = "FABIA loadings omic1", cex = 0.3)
plot(res[["fit"]][["loadings_per_block"]][["omic2"]], main = "FABIA loadings omic2", cex = 0.3)
par(mfrow = c(1, 1))

# Results Quick Access
# res$fit$factor_scores           # Estimated factor scores
# res$fit$loadings_per_block      # List of loadings per block
# res$evaluation$summary          # Summary metrics for matched factors