# ==============================================================================
# Description: Simulate data, fit FABIA, and evaluate score & loading recovery.
#
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
  snr = 0.05,
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
  snr = 0.05,
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


# --- PLOT: Raw Simulation Heatmap ---
plot_simData(sim_object = sim_object, type = "heatmap")


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
# 7. Evaluation Functions 
# ------------------------------------------------------------------------------


evaluate_scores <- function(true_mat, est_mat, n_factors_eval = NULL, plots = TRUE, title_prefix = "") {
  true_mat <- as.matrix(true_mat)
  est_mat  <- as.matrix(est_mat)
  
  K_true <- ncol(true_mat); K_est <- ncol(est_mat)
  n_f <- n_factors_eval %||% max(K_true, K_est) 
  
  # Ensure we look at the top factors requested
  trueK <- true_mat[, seq_len(min(K_true, n_f)), drop = FALSE]
  estK  <- est_mat[, seq_len(min(K_est, n_f)), drop = FALSE]
  
  # Basic Correlation
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  if (is.null(dim(cor_mat))) cor_mat <- matrix(cor_mat, nrow = 1, ncol = 1)
  cor_mat[is.na(cor_mat)] <- 0
  
  # Match via Hungarian Algorithm
  dim_max <- max(dim(cor_mat))
  cost_mat <- matrix(1, nrow = dim_max, ncol = dim_max) 
  nr <- nrow(cor_mat); nc <- ncol(cor_mat)
  cost_mat[1:nr, 1:nc] <- 1 - abs(cor_mat)
  
  assignment <- clue::solve_LSAP(cost_mat)
  cols_full <- as.integer(assignment)
  
  # Extract matches relevant to original matrix
  cols <- cols_full[1:nr] 
  cols[cols > nc] <- NA
  
  # Calculate matched signs/corrs
  matched_corrs <- numeric(nr)
  matched_signs <- numeric(nr)
  
  for (i in seq_len(nr)) {
    j <- cols[i]
    if (!is.na(j)) {
      cval <- cor_mat[i, j]
      matched_corrs[i] <- cval
      matched_signs[i] <- if (sign(cval) == 0) 1 else sign(cval)
    } else {
      matched_corrs[i] <- 0
      matched_signs[i] <- 1
    }
  }
  
  # Construct Signed Estimate Matrix
  est_matched_signed <- matrix(0, nrow=nrow(estK), ncol=nr)
  rownames(est_matched_signed) <- rownames(estK)
  colnames(est_matched_signed) <- paste0("Est_matched_", seq_len(nr))
  
  for(i in seq_len(nr)){
    j <- cols[i]
    if(!is.na(j)){
      est_matched_signed[, i] <- estK[, j] * matched_signs[i]
    }
  }
  
  # --- PLOTTING LOGIC (Restored) ---
  if (plots) {
    # 1. Heatmap of Score Correlations
    tmp <- cor_mat
    # Handle constant values for heatmap breaks
    if (length(unique(as.vector(tmp))) == 1) tmp <- tmp + 1e-8
    
    rownames(tmp) <- paste0("True_F", seq_len(nrow(tmp)))
    colnames(tmp) <- paste0("Est_F", seq_len(ncol(tmp)))
    
    safe_pheatmap(
      tmp,
      main = paste0(title_prefix, " Score Correlation (true x est)"),
      display_numbers = TRUE,
      number_color = "black",
      fontsize_number = 10,
      cellwidth = 30,  # Consistent size with loading heatmaps
      cellheight = 30,
      treeheight_row = 0,
      treeheight_col = 0
    )
    
    # 2. Scatter Plots (True vs Matched Estimated Scores)
    # Only plot the valid true factors (ignore phantom padding if it existed)
    K_plot <- ncol(trueK)
    df_plot <- data.frame(
      true = as.vector(trueK),
      est = as.vector(est_matched_signed[, 1:K_plot, drop=FALSE]),
      factor = factor(rep(paste0("Factor ", seq_len(K_plot)), each = nrow(trueK)))
    )
    
    print(
      ggplot(df_plot, aes(x = true, y = est)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~factor, scales = "free") +
        labs(title = paste0(title_prefix, " True vs Estimated Scores")) +
        theme_bw()
    )
  }
  
  # Return
  list(
    correlation_matrix = cor_mat,
    assignment = cols, 
    matched_signs = matched_signs,
    est_signed_matched = est_matched_signed,
    true_matched = trueK,
    unassigned_est_indices = setdiff(seq_len(ncol(estK)), cols[!is.na(cols)])
  )
}

evaluate_loadings_per_block <- function(true_loadings_list, est_loadings_list, 
                                        factor_assignment, matched_signs, unassigned_est_indices,
                                        plots = TRUE, title_prefix = "") {
  blocks <- names(true_loadings_list)
  res_block <- list()
  
  for (block in blocks) {
    true_mat <- as.matrix(true_loadings_list[[block]])
    est_mat  <- as.matrix(est_loadings_list[[block]])
    
    # --- 1. Alignment ---
    if (!is.null(rownames(true_mat)) && !is.null(rownames(est_mat)) && all(rownames(true_mat) %in% rownames(est_mat))) {
      est_mat <- est_mat[rownames(true_mat), , drop = FALSE]
    } 
    
    # --- 2. Determine Target Dimensions (Square p x p) ---
    p_model <- ncol(est_mat) 
    
    # --- 3. Construct "Est_Final" (Aligned + Unmatched) ---
    est_final <- matrix(0, nrow = nrow(est_mat), ncol = p_model)
    colnames(est_final) <- paste0("Est_F", seq_len(p_model))
    
    assigned_count <- length(factor_assignment)
    
    # A. Fill aligned columns
    for (i in seq_along(factor_assignment)) {
      est_idx <- factor_assignment[i]
      if (!is.na(est_idx) && est_idx <= p_model) {
        est_final[, i] <- est_mat[, est_idx] * matched_signs[i]
      }
    }
    
    # B. Fill remaining columns with Unmatched Estimated factors
    current_cols <- seq_len(p_model)
    filled_cols  <- seq_len(assigned_count)
    remaining_slots <- setdiff(current_cols, filled_cols)
    leftover_est_indices <- unassigned_est_indices
    
    if (length(remaining_slots) > 0 && length(leftover_est_indices) > 0) {
      cnt <- min(length(remaining_slots), length(leftover_est_indices))
      for (k in 1:cnt) {
        slot <- remaining_slots[k]
        est_idx <- leftover_est_indices[k]
        est_final[, slot] <- est_mat[, est_idx]
      }
    }
    
    # --- 4. Pad "True_Mat" with zeros (Phantoms) ---
    if (ncol(true_mat) < p_model) {
      missing_n <- p_model - ncol(true_mat)
      zeros <- matrix(0, nrow = nrow(true_mat), ncol = missing_n)
      colnames(zeros) <- paste0("Phantom_True_", seq_len(missing_n))
      true_mat_padded <- cbind(true_mat, zeros)
    } else {
      true_mat_padded <- true_mat
    }
    true_mat_padded <- true_mat_padded[, 1:p_model, drop=FALSE]
    
    # --- 5. Visualization ---
    if (plots) {
      # We suppress warnings here because 'true_mat_padded' contains columns of zeros.
      # cor() correctly warns that SD is zero for those columns. We expect this.
      cor_block <- suppressWarnings(cor(true_mat_padded, est_final, use = "pairwise.complete.obs"))
      
      # Replace NA (from zero-variance correlations) with 0
      cor_block[is.na(cor_block)] <- 0
      
      rownames(cor_block) <- paste0("True_F", seq_len(nrow(cor_block)))
      colnames(cor_block) <- paste0("Est_F", seq_len(ncol(cor_block)))
      
      # Force heatmap dimensions
      safe_pheatmap(
        cor_block,
        main = paste0(title_prefix, " Loadings Corr. (", block, ")"),
        display_numbers = TRUE,
        number_color = "black",
        fontsize_number = 9,
        cellwidth = 30,  
        cellheight = 30, 
        treeheight_row = 0, 
        treeheight_col = 0
      )
      
      # Scatter plots (only for valid factors)
      valid_k <- ncol(true_mat)
      df_plot <- bind_rows(lapply(seq_len(valid_k), function(k) {
        data.frame(
          feature = rownames(true_mat),
          true_loading = true_mat[, k],
          est_loading = est_final[, k],
          factor = paste0("Factor", k),
          stringsAsFactors = FALSE
        )
      }))
      
      print(
        ggplot(df_plot, aes(x = true_loading, y = est_loading)) +
          geom_point(alpha = 0.4, size = 0.6) +
          geom_smooth(method = "lm", se = FALSE) +
          facet_wrap(~factor, scales = "free") +
          labs(title = paste0(title_prefix, " True vs Est Loadings (", block, ")")) +
          theme_bw()
      )
    }
    
    res_block[[block]] <- list(summary = "See heatmap", est_final = est_final)
  }
  
  res_block
}

# ------------------------------------------------------------------------------
# 8. Main Wrapper: Fit & Evaluate
# ------------------------------------------------------------------------------

fit_and_evaluate_fabia <- function(Xlist, sim_object = NULL, p = 3, alpha = 0.05, seed = 123,
                                   plots = TRUE, n_factors_eval = NULL, verbose = TRUE) {
  
  # 1. Fit
  fit_res <- fit_fabia(Xlist, p = p, alpha = alpha, seed = seed, verbose = verbose)
  
  # 2. Evaluate (if sim_object provided)
  eval_res <- NULL
  if (!is.null(sim_object) && !is.null(sim_object$list_alphas)) {
    mats <- sim_object$list_alphas
    true_mat <- do.call(cbind, lapply(mats, as.numeric))
    colnames(true_mat) <- names(mats)
    rn <- rownames(sim_object$omics[[1]]) %||% seq_len(nrow(true_mat))
    rownames(true_mat) <- rn
    
    # Pass n_factors_eval (or default to p to ensure square heatmap)
    n_eval <- n_factors_eval %||% p
    
    # Score Evaluation
    eval_res <- evaluate_scores(
      true_mat,
      fit_res$factor_scores,
      n_factors_eval = n_eval,
      plots = plots,
      title_prefix = "FABIA"
    )
    
    # Loading Evaluation
    if (!is.null(sim_object$list_betas)) {
      true_loadings_list <- get_true_loadings_per_block(sim_object)
      
      eval_res$loadings_per_block <- evaluate_loadings_per_block(
        true_loadings_list,
        fit_res$loadings_per_block,
        factor_assignment = eval_res$assignment,
        matched_signs = eval_res$matched_signs,
        unassigned_est_indices = eval_res$unassigned_est_indices, # Argument ensures correct heatmap dims
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
# 9. Execution
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

