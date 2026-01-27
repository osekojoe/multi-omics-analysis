# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================

# Required packages 
required_pkgs <- c("SUMO", "FactoMineR", "fabia", "clue", "pheatmap", 
                   "ggplot2", "matrixStats", "reshape2", "e1071")
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install)) install.packages(to_install)

library(SUMO)
library(FactoMineR)
library(fabia)
library(clue)
library(pheatmap)
library(ggplot2)
library(matrixStats)
library(reshape2)
library(e1071) # Required for kurtosis calculation

# ==============================================================================
# 2. Data Simulation 
# ==============================================================================

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


# ==============================================================================
# 3. CORE ALGORITHM: Kurtosis-Weighted MFA-FABIA
# ==============================================================================

#' Run FABIA with feature weighting based on Kurtosis (Non-Normality)
#'
#' @param Xlist List of matrices (omics data).
#' @param gamma Tempering parameter for weights (0 = no weighting).
#' @param p FABIA parameter (number of biclusters).
#' @param alpha FABIA parameter (sparseness). Lower values = denser loadings.
#' @param block_equalize Boolean. Equalize block variance before MFA?
#' @param scale_vars Boolean. Scale variables?
#' @param seed Random seed.
#' @return List containing MFA results, weighted data, FABIA results, and weights.
mfa_weighted_fabia <- function(Xlist,
                               gamma = 0.5,
                               p = 3, 
                               alpha = 0.5, # Default lowered to allow denser loadings
                               block_equalize = TRUE,
                               scale_vars = TRUE,
                               seed = 1) {
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  stopifnot(is.list(Xlist), length(Xlist) >= 2)
  set.seed(seed)
  
  # --- Helper: clean NA/Inf per block ---
  clean_block <- function(M) {
    M <- as.matrix(M)
    all_bad <- apply(M, 2, function(z) all(!is.finite(z)))
    if (any(all_bad)) M <- M[, !all_bad, drop = FALSE]
    
    for (j in seq_len(ncol(M))) {
      z <- M[, j]
      bad <- !is.finite(z)
      if (any(bad)) {
        m <- mean(z[is.finite(z)], na.rm = TRUE)
        if (!is.finite(m)) m <- 0
        z[bad] <- m
        M[, j] <- z
      }
    }
    M
  }
  
  # --- 1) Preprocess ---
  Zlist <- lapply(Xlist, function(X) {
    X <- as.matrix(X)
    if (scale_vars) X <- scale(X) else X <- scale(X, scale = FALSE)
    X
  })
  
  if (block_equalize) {
    Zlist <- Map(function(Z, k) Z / sqrt(k), Zlist, lapply(Zlist, ncol))
  }
  
  Zlist <- lapply(Zlist, clean_block)
  
  names(Zlist) <- names(Zlist) %||% paste0("omic", seq_along(Zlist))
  for (g in seq_along(Zlist)) {
    if (is.null(colnames(Zlist[[g]]))) {
      colnames(Zlist[[g]]) <- paste0("V", seq_len(ncol(Zlist[[g]])))
    }
    colnames(Zlist[[g]]) <- paste0(names(Zlist)[g], "_", colnames(Zlist[[g]]))
  }
  
  Z   <- do.call(cbind, Zlist)              
  grp <- vapply(Zlist, ncol, integer(1))    
  
  # --- 2) MFA (Structure Preservation) ---
  # We still run MFA to preserve the structural information object
  # even though we switched weighting to Kurtosis.
  df_all <- as.data.frame(Z)
  mfa <- FactoMineR::MFA(
    df_all,
    group      = grp,
    type       = rep("s", length(grp)),
    ncp        = 5,
    name.group = names(Zlist),
    graph      = FALSE
  )
  
  # --- 3) Build feature weights (KURTOSIS FIX) ---
  # Replaces previous cos2 method.
  # Logic: High Kurtosis = Heavy Tails = Informative for FABIA (Sparse Signal)
  #        Low Kurtosis  = Gaussian    = Noise
  
  message("Calculating Kurtosis for feature weighting...")
  
  # Calculate Excess Kurtosis (Normal = 0)
  kurt_vals <- apply(Z, 2, e1071::kurtosis, type = 2, na.rm = TRUE)
  
  # FABIA loves positive excess kurtosis (heavy tails).
  # We floor negative kurtosis (platykurtic/flat) to 0 or a small epsilon
  # so we don't zero out weights completely but de-prioritize them significantly.
  w_raw <- pmax(0, kurt_vals) 
  
  # Normalize weights
  mean_w <- mean(w_raw, na.rm = TRUE)
  if (mean_w == 0) mean_w <- 1 # Safety check
  w <- w_raw / mean_w
  
  # Handle potential NAs
  w[!is.finite(w)] <- 0
  
  # Apply tempering
  w <- w^gamma
  
  # --- 4) Apply weights ---
  Z_w <- sweep(Z, 2, w, `*`)    
  
  cuts <- cumsum(grp)
  weighted_blocks <- lapply(seq_along(grp), function(i) {
    if (i == 1) Z_w[, 1:cuts[1], drop = FALSE]
    else Z_w[, (cuts[i - 1] + 1):cuts[i], drop = FALSE]
  })
  names(weighted_blocks) <- names(Zlist)
  
  # --- 5) FABIA on weighted data ---
  # Using standard fabia. 
  # Note: Lower alpha results in denser loadings (less sparsity).
  fab_fit <- fabia::fabia(as.matrix(t(Z_w)), p = p, alpha = alpha)
  
  factor_scores <- t(fab_fit@Z)    
  loadings_all  <- fab_fit@L       
  
  idx <- split(seq_len(nrow(loadings_all)), rep(seq_along(grp), times = grp))
  loadings_per_block <- lapply(idx, function(ix) loadings_all[ix, , drop = FALSE])
  names(loadings_per_block) <- names(Zlist)
  
  fab_res <- list(
    FABIA              = fab_fit,
    factor_scores      = factor_scores,
    loadings_per_block = loadings_per_block,
    group_sizes        = grp
  )
  
  list(
    mfa             = mfa,
    weighted        = Z_w,
    weighted_blocks = weighted_blocks,
    fabia           = fab_res,
    weights         = w,
    kurtosis_raw    = kurt_vals
  )
}

# ==============================================================================
# 4. EVALUATION HELPERS (UNCHANGED)
# ==============================================================================

get_true_factors_from_sim <- function(sim_object) {
  if (!is.null(sim_object$list_alphas) && is.list(sim_object$list_alphas)) {
    mats <- sim_object$list_alphas
    mat <- do.call(cbind, lapply(mats, as.numeric))
    colnames(mat) <- names(mats)
    rownames(mat) <- if (!is.null(rownames(sim_object$omics[[1]]))) rownames(sim_object$omics[[1]]) else seq_len(nrow(mat))
    return(mat)
  }
  stop("Could not find true factors in sim_object$list_alphas.")
}

evaluate_factor_recovery <- function(sim_object, res, n_factors = NULL, plot = TRUE, return_all = FALSE) {
  
  if (is.list(res) && !is.null(res$fabia) && !is.null(res$fabia$factor_scores)) {
    est <- as.matrix(res$fabia$factor_scores)
  } else if (!is.null(res$factor_scores)) {
    est <- as.matrix(res$factor_scores)
  } else {
    stop("Cannot find estimated factor scores.")
  }
  n <- nrow(est)
  true <- as.matrix(get_true_factors_from_sim(sim_object))
  
  K_est <- ncol(est); K_true <- ncol(true)
  if (is.null(n_factors)) n_factors <- min(K_est, K_true)
  
  estK <- est[, seq_len(min(ncol(est), n_factors)), drop = FALSE]
  trueK <- true[, seq_len(min(ncol(true), n_factors)), drop = FALSE]
  K <- ncol(trueK)
  
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  cor_mat_na <- cor_mat; cor_mat_na[is.na(cor_mat_na)] <- 0
  
  # Hungarian Assignment
  cost <- 1 - abs(cor_mat_na)
  if (any(cost < 0)) cost <- cost - min(cost)
  assignment <- clue::solve_LSAP(as.matrix(cost))
  assigned_cols <- as.integer(assignment)
  
  matched_corrs <- numeric(K); matched_signs <- numeric(K)
  for (i in seq_len(K)) {
    j <- assigned_cols[i]
    cval <- if (is.na(cor_mat[i, j])) 0 else cor_mat[i, j]
    matched_corrs[i] <- cval
    matched_signs[i] <- if (sign(cval) == 0) 1 else sign(cval)
  }
  
  est_matched_signed <- estK[, assigned_cols, drop = FALSE] %*% diag(matched_signs)
  colnames(est_matched_signed) <- paste0("Est_matched_", seq_len(K))
  
  # Metrics
  r2 <- numeric(K); mse <- numeric(K)
  for (i in seq_len(K)) {
    fit <- lm(trueK[, i] ~ est_matched_signed[, i])
    ssr <- sum((predict(fit) - mean(trueK[, i]))^2)
    sst <- sum((trueK[, i] - mean(trueK[, i]))^2)
    r2[i] <- if (sst == 0) NA else ssr / sst
    mse[i] <- mean((trueK[, i] - est_matched_signed[, i])^2)
  }
  
  if (plot) {
    pheatmap(cor_mat, main = "Corr. score (true x estimated factors)",
             cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
    
    df_list <- lapply(seq_len(K), function(i) {
      data.frame(sample = seq_len(n), true = trueK[, i], est = est_matched_signed[, i], factor = paste0("True_", i))
    })
    g <- ggplot(do.call(rbind, df_list), aes(x = true, y = est)) +
      geom_point(alpha = 0.6) + facet_wrap(~factor, scales = "free") +
      geom_smooth(method = "lm", se = FALSE) +
      labs(title = "True vs Estimated Factors (Matched & Signed)")
    print(g)
  }
  
  summary_tbl <- data.frame(
    true_factor = colnames(trueK), matched_idx = assigned_cols,
    corr = matched_corrs, R2 = r2, MSE = mse, stringsAsFactors = FALSE
  )
  
  res_list <- list(cor_mat = cor_mat, summary = summary_tbl, mean_abs_cor = mean(abs(matched_corrs)))
  if (return_all) {
    res_list$est_signed_matched <- est_matched_signed
    res_list$true_matched <- trueK
  }
  return(res_list)
}

evaluate_loading_recovery <- function(sim_object, fab_res, verbose = TRUE) {
  K <- length(sim_object$list_alphas)
  block_names <- names(sim_object$omics)
  block_sizes <- vapply(sim_object$omics, ncol, integer(1))
  p_total <- sum(block_sizes)
  
  true_all <- matrix(0, nrow = p_total, ncol = K)
  colnames(true_all) <- paste0("factor", seq_len(K))
  row_start <- 1
  for (b in seq_along(block_names)) {
    ncols <- block_sizes[b]
    row_end <- row_start + ncols - 1
    betas <- sim_object$list_betas[[b]]
    
    if (!is.null(betas) && is.list(betas)) {
      if (is.null(names(betas))) {
        for (j in seq_along(betas)) if (j <= K) true_all[row_start:row_end, j] <- as.numeric(betas[[j]])
      } else {
        for (nm in names(betas)) {
          m <- regmatches(nm, regexec("beta(\\d+)$", nm))
          if (length(m) && length(m[[1]]) >= 2) {
            target <- as.integer(m[[1]][2])
            if (target <= K) true_all[row_start:row_end, target] <- as.numeric(betas[[nm]])
          }
        }
      }
    } else if (!is.null(betas)) {
      true_all[row_start:row_end, 1] <- as.numeric(betas)
    }
    row_start <- row_end + 1
  }
  
  est_all <- as.matrix(fab_res$FABIA@L)
  K_est <- ncol(est_all); K_use <- min(K, K_est)
  
  cor_mat <- cor(true_all[, 1:K_use], est_all[, 1:K_use], use = "pairwise.complete.obs")
  abs_cor <- abs(cor_mat); abs_cor[!is.finite(abs_cor)] <- 0
  
  cost <- max(abs_cor) - abs_cor
  assignment <- clue::solve_LSAP(cost)
  assigned_pairs <- cbind(seq_len(nrow(cost)), as.integer(assignment))
  
  matched <- data.frame()
  est_all_flipped <- est_all
  for (i in seq_len(nrow(assigned_pairs))) {
    tr <- assigned_pairs[i, 1]; es <- assigned_pairs[i, 2]
    cc <- suppressWarnings(cor(true_all[, tr], est_all[, es]))
    if (!is.finite(cc)) cc <- 0
    if (cc < 0) { est_all_flipped[, es] <- -est_all_flipped[, es]; cc <- -cc }
    rmse <- sqrt(mean((true_all[, tr] - est_all_flipped[, es])^2))
    matched <- rbind(matched, data.frame(true_comp = tr, est_comp = es, cor = cc, rmse = rmse))
  }
  
  list(true_all = true_all, est_all_flipped = est_all_flipped, cor_mat = cor_mat, 
       matched = matched, block_sizes = block_sizes)
}


#' Per-Omic Loading Recovery Analysis
#' 
#' Breaks down global loading recovery into specific omics layers.
#' Handles zero-variance vectors (common in sparse estimation) robustly.
per_omic_loading_recovery <- function(eval_res, sim_object, plot_heatmaps = TRUE, plot_scatter = TRUE, scatter_n = 5) {
  
  # Extract global matrices from the previous evaluation result
  true_all <- eval_res$true_all
  est_all_flipped <- eval_res$est_all_flipped
  matched <- eval_res$matched
  block_sizes <- eval_res$block_sizes
  block_names <- names(block_sizes)
  
  # Define indices for each block
  cuts <- cumsum(block_sizes)
  idx_list <- lapply(seq_along(block_sizes), function(i) {
    if (i == 1) seq_len(cuts[1]) else (cuts[i - 1] + 1):cuts[i]
  })
  names(idx_list) <- block_names
  
  per_block_summary <- list()
  
  for (bname in block_names) {
    idx <- idx_list[[bname]]
    true_block <- true_all[idx, , drop = FALSE]
    est_block <- est_all_flipped[idx, , drop = FALSE]
    
    # --- ROBUST CORRELATION CALCULATION ---
    # 1. Identify columns with Zero SD (Constant vectors)
    sd_true <- apply(true_block, 2, sd, na.rm = TRUE)
    sd_est  <- apply(est_block, 2, sd, na.rm = TRUE)
    
    valid_true <- sd_true > 1e-12
    valid_est  <- sd_est > 1e-12
    
    # 2. Initialize Correlation Matrix with Zeros
    cormat_block <- matrix(0, nrow = ncol(true_block), ncol = ncol(est_block))
    rownames(cormat_block) <- colnames(true_block)
    colnames(cormat_block) <- colnames(est_block)
    
    # 3. Calculate correlation only for valid vectors
    if (any(valid_true) && any(valid_est)) {
      sub_cor <- stats::cor(
        true_block[, valid_true, drop = FALSE], 
        est_block[, valid_est, drop = FALSE], 
        use = "pairwise.complete.obs"
      )
      cormat_block[valid_true, valid_est] <- sub_cor
    }
    cormat_block[!is.finite(cormat_block)] <- 0
    
    # --- PLOTTING ---
    # 1. Heatmap
    if (plot_heatmaps) {
      pheatmap(cormat_block, 
               cluster_rows = FALSE, cluster_cols = FALSE, 
               display_numbers = TRUE,
               main = paste0("Correlation: ", bname))
    }
    
    # 2. Scatter Plots (Matched Factors Only)
    # Filter matches relevant to this block
    matched_in_block <- do.call(rbind, lapply(seq_len(nrow(matched)), function(i) {
      tcomp <- matched$true_comp[i]
      ecomp <- matched$est_comp[i]
      cor_val <- cormat_block[tcomp, ecomp]
      data.frame(true_comp = tcomp, est_comp = ecomp, cor_block = cor_val)
    }))
    
    if (plot_scatter && nrow(matched_in_block) > 0) {
      # Limit to top N scatter plots if requested
      to_plot <- matched_in_block
      if (!is.null(scatter_n)) to_plot <- to_plot[seq_len(min(nrow(to_plot), scatter_n)), ]
      
      for (ri in seq_len(nrow(to_plot))) {
        tc <- to_plot$true_comp[ri]; ec <- to_plot$est_comp[ri]
        
        df <- data.frame(true = true_block[, tc], est = est_block[, ec])
        
        p <- ggplot(df, aes(x = true, y = est)) +
          geom_point(alpha = 0.4) + 
          geom_smooth(method = "lm", se = FALSE, color = "blue") +
          labs(title = sprintf("%s: True Factor %d vs Est %d", bname, tc, ec),
               subtitle = sprintf("Correlation in this block: %.2f", to_plot$cor_block[ri]),
               x = "True Loading", y = "Estimated Loading") +
          theme_minimal()
        print(p)
      }
    }
    per_block_summary[[bname]] <- list(cormat = cormat_block, matched = matched_in_block)
  }
  return(per_block_summary)
}

# ==============================================================================
# 5. EXECUTION
# ==============================================================================

if (interactive()) {
  # --- A. Data Simulation ---
  # (Simulated data logic is preserved from Section 2)
  sim_object <- sim1
  sim_object$omics$omic2 <- sim2$omics$omic2
  sim_object[["list_betas"]][[2]] <- sim2[["list_betas"]][[2]]
  simX_list <- sim_object$omics
  
  # PLOT 1: Data Simulation Heatmap
  plot_simData(sim_object = sim_object, type = "heatmap")
  
  # --- B. Run Method (FIXED) ---
  # Parameters set for:
  # 1. Kurtosis Weighting (Internal)
  # 2. Dense Loadings (alpha = 0.01)
  # 3. p = 3 (number of factors in simulation)
  res_kurtosis <- mfa_weighted_fabia(
    simX_list, 
    gamma = 0.5,    # Standard tempering for kurtosis
    p = 3, 
    alpha = 0.01    # VERY LOW alpha to ensure non-zero (dense) loadings
  )
  
  # PLOT 2: Basic FABIA plots
  plot(res_kurtosis$fabia$FABIA, main = "FABIA Results (Kurtosis Weighted)")
  
  # PLOT 3: Factor Scores & Loadings
  factor_scores <- res_kurtosis$fabia$factor_scores
  factor_loadings <- res_kurtosis$fabia$loadings_per_block
  plot(factor_scores, main = "Factor Scores")
  
  # --- C. Factor Recovery Evaluation ---
  message("Evaluating Factors...")
  eval_factors <- evaluate_factor_recovery(sim_object, res_kurtosis, n_factors = 3, plot = TRUE, return_all = TRUE)
  print(eval_factors$summary)
  
  # --- D. Loading Recovery Evaluation (Global) ---
  message("Evaluating Loadings...")
  eval_loadings <- evaluate_loading_recovery(sim_object, res_kurtosis$fabia, verbose = TRUE)
  
  # Check density of results (Verification of non-zero request)
  sparsity_pct <- sum(abs(res_kurtosis$fabia$FABIA@L) > 1e-4) / prod(dim(res_kurtosis$fabia$FABIA@L)) * 100
  message(sprintf("Loading Density (%% non-zero): %.2f%%", sparsity_pct))
  
  # PLOT 4: Global Loading Correlation Heatmap
  pheatmap(eval_loadings$cor_mat, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           display_numbers = TRUE, 
           main = "Global Loading Correlation")
  
  # --- E. Per-Omic Loading Recovery (NEW) ---
  message("Plotting Per-Omic Correlations...")
  
  # PLOT 5 & 6: Specific Heatmaps and Scatters for Omic1 and Omic2
  per_block_res <- per_omic_loading_recovery(
    eval_res = eval_loadings, 
    sim_object = sim_object,
    plot_heatmaps = TRUE, 
    plot_scatter = TRUE, 
    scatter_n = 3 # Plot top 3 matched factors per block
  )
  
  # Check results
  print(per_block_res$omic1$matched)
  print(per_block_res$omic2$matched)
}


## --- Feature Selection Overlap (Jaccard Index)
evaluate_feature_jaccard <- function(sim_object, res, top_n = 500) {
  # 1. Get True High-Weight Features
  # Extract true betas from simulation
  eval_load <- evaluate_loading_recovery(sim_object, res$fabia, verbose=FALSE)
  true_L <- eval_load$true_all
  est_L  <- eval_load$est_all_flipped
  
  # Use the matched pairs from previous evaluation
  matches <- eval_load$matched
  
  jaccards <- numeric(nrow(matches))
  
  for(i in 1:nrow(matches)) {
    t_idx <- matches$true_comp[i]
    e_idx <- matches$est_comp[i]
    
    # Get indices of top N features by absolute weight
    top_true <- order(abs(true_L[, t_idx]), decreasing = TRUE)[1:top_n]
    top_est  <- order(abs(est_L[, e_idx]), decreasing = TRUE)[1:top_n]
    
    # Calculate Intersection over Union
    intersection <- length(intersect(top_true, top_est))
    union_set    <- length(union(top_true, top_est))
    
    jaccards[i] <- intersection / union_set
  }
  
  names(jaccards) <- paste0("Factor_", matches$true_comp)
  return(jaccards)
}

# Usage:
j_scores <- evaluate_feature_jaccard(sim_object, res_kurtosis, top_n = 300)
barplot(j_scores, main="Feature Selection Overlap (Jaccard)", ylim=c(0,1))


## Classification Accuracy (AUC-ROC)

library(pROC)

evaluate_sample_auc <- function(sim_object, res, factor_idx_matched) {
  # 1. Get True Binary Labels (Active vs Inactive samples)
  # In your simulation, 'list_alphas' contains the true factor activity
  true_factors <- get_true_factors_from_sim(sim_object)
  
  # 2. Get Estimated Scores
  est_scores <- res$fabia$factor_scores
  
  aucs <- numeric()
  
  for(i in seq_along(factor_idx_matched)) {
    # Match true factor i with estimated factor j (from your Hungarian result)
    est_col <- factor_idx_matched[i]
    
    # Create binary label: Signal is present if true value is non-zero (or above threshold)
    # Note: simulateMultiOmics usually creates continuous factors, 
    # but we can threshold them to see "Active" samples if the design is sparse.
    # Alternatively, use the raw true values for correlation, but here we want classification.
    true_binary <- abs(true_factors[, i]) > 1e-3 
    
    if(length(unique(true_binary)) > 1) {
      # Calculate ROC
      roc_obj <- pROC::roc(true_binary, est_scores[, est_col], quiet=TRUE)
      aucs[i] <- roc_obj$auc
    } else {
      aucs[i] <- NA # Skip if factor is active in ALL or NO samples
    }
  }
  return(aucs)
}

# Usage:
matched_indices <- eval_factors$summary$matched_idx 
auc_scores <- evaluate_sample_auc(sim_object, res_kurtosis, matched_indices)
print(auc_scores)


# ---- Reconstruction Error (RMSE)
evaluate_reconstruction <- function(Xlist, res) {
  # 1. Combine Original Data (Centered/Scaled)
  # We must replicate the preprocessing used inside the function
  Zlist <- lapply(Xlist, function(x) scale(as.matrix(x)))
  Zlist <- Map(function(Z, k) Z / sqrt(k), Zlist, lapply(Zlist, ncol))
  X_orig <- do.call(cbind, Zlist)
  
  # 2. Reconstruct Data from Model
  # X_hat = Scores * Loadings^T
  Scores <- as.matrix(res$fabia$factor_scores)
  Loadings <- as.matrix(res$fabia$FABIA@L)
  
  X_hat <- Scores %*% t(Loadings)
  
  # 3. Calculate Residuals
  Residuals <- X_orig - X_hat
  RMSE <- sqrt(mean(Residuals^2))
  
  # Optional: Variance Explained
  var_total <- sum(X_orig^2)
  var_resid <- sum(Residuals^2)
  R2_global <- 1 - (var_resid / var_total)
  
  return(list(RMSE = RMSE, R2_global = R2_global))
}

# Usage:
recon_stats <- evaluate_reconstruction(simX_list, res_kurtosis)
cat(sprintf("Global Variance Explained: %.2f%%\n", recon_stats$R2_global * 100))


# ==============================================================================
# 8. PER-FACTOR SNR TREND ANALYSIS
# ==============================================================================

# 1. Define range: 0.01 to 2.0 with steps of 0.05
snr_values <- seq(from = 0.01, to = 2.0, by = 0.4)

# Initialize storage
factor_results <- data.frame(
  SNR = numeric(),
  Omic = character(),
  Factor = character(),
  Best_Correlation = numeric(),
  stringsAsFactors = FALSE
)

message(sprintf("Analyzing Factor-Specific Recovery across %d SNR steps...", length(snr_values)))
pb <- txtProgressBar(min = 0, max = length(snr_values), style = 3)

for (i in seq_along(snr_values)) {
  s <- snr_values[i]
  setTxtProgressBar(pb, i)
  
  # --- A. Re-simulate Data ---
  # We use the same structure as before
  sim1_s <- simulateMultiOmics(
    vector_features = c(4000, 2500), n_samples = 100, n_factors = 3,
    snr = s, 
    signal.samples = c(3, 1), signal.features = list(c(4.5, 0.5), c(4.5, 0.5)),
    factor_structure = "mixed", num.factor = "multiple", seed = 123
  )
  sim2_s <- simulateMultiOmics(
    vector_features = c(4000, 3500), n_samples = 100, n_factors = 3,
    snr = s,
    signal.samples = c(3, 1), signal.features = list(c(2.5, 0.5), c(3, 2.5)),
    factor_structure = "mixed", num.factor = "multiple", seed = 123
  )
  
  sim_obj_s <- sim1_s
  sim_obj_s$omics$omic2 <- sim2_s$omics$omic2
  sim_obj_s[["list_betas"]][[2]] <- sim2_s[["list_betas"]][[2]]
  simX_s <- sim_obj_s$omics
  names(simX_s) <- paste0("omic", seq_along(simX_s))
  
  # --- B. Run Sparse MFA-FABIA ---
  res_s <- suppressMessages(
    mfa_weighted_fabia(simX_s, gamma = 0.5, p = 3, alpha = 0.01, seed = 1)
  )
  
  # --- C. Get Correlation Matrices ---
  eval_load_s <- evaluate_loading_recovery(sim_obj_s, res_s$fabia, verbose = FALSE)
  per_omic_s <- per_omic_loading_recovery(
    eval_res = eval_load_s, sim_object = sim_obj_s, 
    plot_heatmaps = FALSE, plot_scatter = FALSE
  )
  
  # --- D. Extract Max Correlation PER FACTOR ---
  
  # Function to extract best match for each row (True Factor)
  extract_best_matches <- function(omic_res, omic_name, snr_val) {
    cormat <- omic_res$cormat
    # Iterate through True Factors (Rows of the correlation matrix)
    if (!is.null(cormat) && nrow(cormat) > 0) {
      for (f_idx in 1:nrow(cormat)) {
        # Find the single highest absolute correlation in this row
        best_val <- max(abs(cormat[f_idx, ]), na.rm = TRUE)
        
        # Append to global dataframe
        factor_results <<- rbind(factor_results, data.frame(
          SNR = snr_val,
          Omic = omic_name,
          Factor = rownames(cormat)[f_idx], # e.g., "factor1"
          Best_Correlation = best_val
        ))
      }
    }
  }
  
  # Run extraction for both omics
  extract_best_matches(per_omic_s$omic1, "Omic1", s)
  extract_best_matches(per_omic_s$omic2, "Omic2", s)
}
close(pb)

# ==============================================================================
# 9. PLOTTING
# ==============================================================================

# Clean factor names for plotting (optional)
factor_results$Factor <- factor(factor_results$Factor, levels = c("factor1", "factor2", "factor3"))

ggplot(factor_results, aes(x = SNR, y = Best_Correlation, color = Factor)) +
  # Facet by Omic so we can see differences between layers
  facet_wrap(~Omic) + 
  geom_smooth(method = "loess", se = FALSE, alpha = 0.2, size = 0.8) +
  #geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("factor1" = "#E41A1C", "factor2" = "#377EB8", "factor3" = "#4DAF4A")) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title = "Differential Factor Recovery Analysis",
    subtitle = "Tracking the single best bicluster match for each True Factor",
    x = "Signal-to-Noise Ratio (SNR)",
    y = "Max Correlation (Recovery Score)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )




##############################################################
# ==============================================================================
# 10. JACCARD INDEX TREND ANALYSIS 
# ==============================================================================

# 1. Define range and parameters
snr_values <- seq(from = 0.01, to = 2.0, by = 0.4)
top_n_features <- 400  # Number of top features to compare

# Initialize storage
jaccard_results <- data.frame(
  SNR = numeric(),
  Omic = character(),
  Factor = character(),
  Jaccard_Index = numeric(),
  stringsAsFactors = FALSE
)

message(sprintf("Analyzing Jaccard Overlap (Top %d) across %d SNR steps...", 
                top_n_features, length(snr_values)))
pb <- txtProgressBar(min = 0, max = length(snr_values), style = 3)

for (i in seq_along(snr_values)) {
  s <- snr_values[i]
  setTxtProgressBar(pb, i)
  
  # --- A. Re-simulate Data ---
  sim1_s <- simulateMultiOmics(
    vector_features = c(4000, 2500), n_samples = 100, n_factors = 3,
    snr = s, 
    signal.samples = c(3, 1), signal.features = list(c(4.5, 0.5), c(4.5, 0.5)),
    factor_structure = "mixed", num.factor = "multiple", seed = 123
  )
  sim2_s <- simulateMultiOmics(
    vector_features = c(4000, 3500), n_samples = 100, n_factors = 3,
    snr = s,
    signal.samples = c(3, 1), signal.features = list(c(2.5, 0.5), c(3, 2.5)),
    factor_structure = "mixed", num.factor = "multiple", seed = 123
  )
  
  sim_obj_s <- sim1_s
  sim_obj_s$omics$omic2 <- sim2_s$omics$omic2
  sim_obj_s[["list_betas"]][[2]] <- sim2_s[["list_betas"]][[2]]
  simX_s <- sim_obj_s$omics
  names(simX_s) <- paste0("omic", seq_along(simX_s))
  
  # --- B. Run Sparse MFA-FABIA ---
  res_s <- suppressMessages(
    mfa_weighted_fabia(simX_s, gamma = 0.5, p = 3, alpha = 0.01, seed = 1)
  )
  
  # --- C. ROBUST Helper to Get True Loadings per Block ---
  # This ensures that "beta2" from simulation always goes to Column 2 (Factor 2)
  # preventing color/label mismatch when Factor 1 is missing.
  get_true_block_robust <- function(sim_obj, block_idx, n_global_factors = 3) {
    betas_list <- sim_obj$list_betas[[block_idx]]
    
    # Determine number of features (rows) from the first available vector
    if (length(betas_list) > 0) {
      n_features <- length(betas_list[[1]])
    } else {
      return(NULL)
    }
    
    # Initialize matrix of ZEROS for ALL global factors (1 to 3)
    mat_out <- matrix(0, nrow = n_features, ncol = n_global_factors)
    colnames(mat_out) <- paste0("factor", 1:n_global_factors)
    
    # Populate columns based on names (e.g., "beta2" -> Column 2)
    if (is.list(betas_list)) {
      names_betas <- names(betas_list)
      for (nm in names_betas) {
        # Extract index from string "betaX"
        if (grepl("beta", nm)) {
          idx <- as.integer(sub("beta", "", nm))
        } else {
          idx <- as.integer(nm) # Fallback
        }
        
        # If valid index, fill the specific column
        if (!is.na(idx) && idx <= n_global_factors) {
          mat_out[, idx] <- as.numeric(betas_list[[nm]])
        }
      }
    } else if (is.matrix(betas_list)) {
      # If already a matrix, fill safely
      k <- min(ncol(betas_list), n_global_factors)
      mat_out[, 1:k] <- betas_list[, 1:k]
    }
    
    return(mat_out)
  }
  
  # --- D. Calculate Jaccard per Omic / per Factor ---
  
  process_jaccard <- function(omic_name, block_idx) {
    # 1. Get True Matrix (Robust)
    true_mat <- get_true_block_robust(sim_obj_s, block_idx, n_global_factors = 3)
    if (is.null(true_mat)) return(NULL)
    
    # Get Estimated Loadings
    if (!omic_name %in% names(res_s$fabia$loadings_per_block)) return(NULL)
    est_mat <- res_s$fabia$loadings_per_block[[omic_name]]
    
    # 2. Match and Calculate
    if (ncol(est_mat) > 0) {
      cormat <- cor(true_mat, est_mat, use="pairwise.complete.obs")
      cormat[is.na(cormat)] <- 0
      
      # Iterate through Global Factors (1, 2, 3)
      for (f in 1:ncol(true_mat)) {
        fname <- colnames(true_mat)[f]
        
        # SKIP empty true factors (if signal is zero, Jaccard is meaningless)
        # We record 0 to keep the plot lines consistent
        if (sd(true_mat[, f]) < 1e-6) {
          jaccard_results <<- rbind(jaccard_results, data.frame(
            SNR = s, Omic = omic_name, Factor = fname, Jaccard_Index = 0
          ))
          next
        }
        
        # Find best estimated match
        best_match_idx <- which.max(abs(cormat[f, ]))
        
        if (length(best_match_idx) > 0) {
          v_true <- true_mat[, f]
          v_est  <- est_mat[, best_match_idx]
          
          top_true <- order(abs(v_true), decreasing = TRUE)[1:top_n_features]
          top_est  <- order(abs(v_est), decreasing = TRUE)[1:top_n_features]
          
          intersect_len <- length(intersect(top_true, top_est))
          union_len     <- length(union(top_true, top_est))
          
          j_val <- if(union_len > 0) intersect_len / union_len else 0
          
          jaccard_results <<- rbind(jaccard_results, data.frame(
            SNR = s, Omic = omic_name, Factor = fname, Jaccard_Index = j_val
          ))
        }
      }
    }
  }
  
  process_jaccard("omic1", 1)
  process_jaccard("omic2", 2)
}
close(pb)

# ==============================================================================
# 11. PLOTTING JACCARD TREND
# ==============================================================================

# Ensure factors are ordered for the legend
jaccard_results$Factor <- factor(jaccard_results$Factor, levels = c("factor1", "factor2", "factor3"))

ggplot(jaccard_results, aes(x = SNR, y = Jaccard_Index, color = Factor)) +
  facet_wrap(~Omic) + 
  geom_smooth(method = "loess", se = FALSE, alpha = 0.2, size = 0.8) +
  # geom_line(alpha = 0.7, size = 1) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("factor1" = "#E41A1C", "factor2" = "#377EB8", "factor3" = "#4DAF4A")) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Feature Selection Overlap (Jaccard Index) vs. SNR",
    subtitle = sprintf("Comparison of Top %d loaded features (True vs Best Estimated Match)", top_n_features),
    x = "Signal-to-Noise Ratio (SNR)",
    y = "Jaccard Index (Intersection / Union)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )


