# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================

# Required packages 
required_pkgs <- c("SUMO", "FactoMineR", "fabia", "clue", "pheatmap", 
                   "ggplot2", "matrixStats", "reshape2", "e1071", "pROC")
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install)) install.packages(to_install)

library(SUMO)         # For simulateMultiOmics
library(FactoMineR)
library(fabia)
library(clue)
library(pheatmap)
library(ggplot2)
library(matrixStats)
library(reshape2)
library(e1071)        # Required for kurtosis calculation
library(pROC)

# ==============================================================================
# 2. DATA SIMULATION 
# ==============================================================================

message("--- Step 1: Simulating Data ---")

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
# plot_simData(sim_object = sim_object, type = "heatmap")


# ==============================================================================
# 3. CORE ALGORITHM: Kurtosis-Weighted MFA-FABIA (Refined)
# ==============================================================================



#' Run FABIA with feature weighting based on Kurtosis and MFA Normalization
#'
#' @param Xlist List of matrices (omics data).
#' @param gamma Tempering parameter for weights (0 = no weighting).
#' @param p FABIA parameter (number of biclusters).
#' @param alpha FABIA parameter (sparseness). Lower values = denser loadings.
#' @param block_equalize Boolean. Equalize block variance via Lambda1?
#' @param scale_vars Boolean. Scale variables?
#' @param seed Random seed.
#' @return List containing MFA results, weighted data, FABIA results, and weights.
mfa_weighted_fabia <- function(Xlist,
                               gamma = 0.5,
                               p = 3, 
                               alpha = 0.5,
                               block_equalize = TRUE,
                               scale_vars = TRUE,
                               seed = 1) {
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  stopifnot(is.list(Xlist), length(Xlist) >= 2)
  set.seed(seed)
  
  # --- Helper: Robust Cleaning ---
  # Remove zero-variance columns and handle Infs
  clean_and_prep <- function(M) {
    M <- as.matrix(M)
    vars <- matrixStats::colVars(M)
    # Remove constant columns (variance approx 0)
    M <- M[, vars > 1e-9, drop = FALSE]
    return(M)
  }
  
  # A. Robust Preprocessing
  Xlist_clean <- lapply(Xlist, clean_and_prep)
  
  # B. MFA Block Equalizing & Scaling
  # Logic matched to Script 1:
  # 1. Scale Variables
  # 2. Divide block by sqrt(First Eigenvalue)
  
  Zlist <- lapply(Xlist_clean, function(X) {
    # 1. Standard Scaling (Mean 0, Sd 1)
    if(scale_vars) {
      s <- scale(X)
      s[is.na(s)] <- 0
    } else {
      s <- scale(X, scale = FALSE)
    }
    
    # 2. MFA Normalization
    if(block_equalize) {
      # Calculate PCA for this block
      pca_block <- prcomp(s, center = FALSE, scale. = FALSE)
      # Get first eigenvalue (lambda 1)
      lambda1 <- pca_block$sdev[1]^2
      
      # Divide by sqrt(lambda1) to equalize inertia of the first dimension
      if(lambda1 > 1e-9) {
        s <- s / sqrt(lambda1)
      }
    }
    return(s)
  })
  
  # Rename features with prefix (Crucial for splitting later)
  names(Zlist) <- names(Zlist) %||% paste0("omic", seq_along(Zlist))
  for(n in names(Zlist)) colnames(Zlist[[n]]) <- paste0(n, "_", colnames(Zlist[[n]]))
  
  # Combine into Global Matrix Z
  Z <- do.call(cbind, Zlist)
  grp <- vapply(Zlist, ncol, integer(1)) 
  
  # --- C. Feature Weighting (Kurtosis) ---
  
  
  message("Calculating Kurtosis for feature weighting...")
  
  # Calculate Excess Kurtosis
  kurt_vals <- apply(Z, 2, e1071::kurtosis, type = 2, na.rm = TRUE)
  
  # Weighting Logic:
  # High kurtosis (heavy tails) -> High weight
  w_raw <- pmax(0, kurt_vals) 
  mean_w <- mean(w_raw, na.rm = TRUE)
  if (mean_w == 0) mean_w <- 1 
  w <- (w_raw / mean_w)^gamma
  w[!is.finite(w)] <- 0
  
  # Apply Weights
  Z_w <- sweep(Z, 2, w, `*`)
  
  # Create weighted blocks list (for return object)
  cuts <- cumsum(grp)
  weighted_blocks <- lapply(seq_along(grp), function(i) {
    if (i == 1) Z_w[, 1:cuts[1], drop = FALSE]
    else Z_w[, (cuts[i - 1] + 1):cuts[i], drop = FALSE]
  })
  names(weighted_blocks) <- names(Zlist)
  
  # --- D. Run FABIA ---
  
  
  # Transpose because FABIA expects samples in columns
  fab_fit <- fabia::fabia(as.matrix(t(Z_w)), p = p, alpha = alpha, random = 1)
  
  factor_scores <- t(fab_fit@Z)    
  loadings_all  <- fab_fit@L       
  
  # Split Loadings per Block based on prefix
  loadings_per_block <- list()
  for(n in names(Zlist)) {
    pattern <- paste0("^", n, "_")
    idx <- grep(pattern, rownames(loadings_all))
    if(length(idx) > 0) {
      loadings_per_block[[n]] <- loadings_all[idx, , drop=FALSE]
    }
  }
  
  fab_res <- list(
    FABIA              = fab_fit,
    factor_scores      = factor_scores,
    loadings_per_block = loadings_per_block,
    group_sizes        = grp
  )
  
  return(list(
    # mfa object is not strictly created here via FactoMineR anymore, 
    # but the manual normalization accomplishes the pre-processing.
    weighted        = Z_w,
    weighted_blocks = weighted_blocks,
    fabia           = fab_res,
    weights         = w,
    kurtosis_raw    = kurt_vals
  ))
}

# ==============================================================================
# 4. EVALUATION HELPERS 
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
  
  # Construct True Loading Matrix
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
  
  # Est loadings
  est_all <- as.matrix(fab_res$FABIA@L)
  K_est <- ncol(est_all); K_use <- min(K, K_est)
  
  # Align via Hungarian on GLOBAL loadings
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
    matched <- rbind(matched, data.frame(true_comp = tr, est_comp = es, cor = cc))
  }
  
  list(true_all = true_all, est_all_flipped = est_all_flipped, cor_mat = cor_mat, 
       matched = matched, block_sizes = block_sizes)
}


#' Per-Omic Loading Recovery Analysis
per_omic_loading_recovery <- function(eval_res, sim_object, plot_heatmaps = TRUE, plot_scatter = TRUE, scatter_n = 5) {
  
  true_all <- eval_res$true_all
  est_all_flipped <- eval_res$est_all_flipped
  matched <- eval_res$matched
  block_sizes <- eval_res$block_sizes
  block_names <- names(block_sizes)
  
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
    
    # Robust Correlation
    sd_true <- apply(true_block, 2, sd, na.rm = TRUE)
    sd_est  <- apply(est_block, 2, sd, na.rm = TRUE)
    
    valid_true <- sd_true > 1e-12
    valid_est  <- sd_est > 1e-12
    
    cormat_block <- matrix(0, nrow = ncol(true_block), ncol = ncol(est_block))
    rownames(cormat_block) <- colnames(true_block)
    colnames(cormat_block) <- colnames(est_block)
    
    if (any(valid_true) && any(valid_est)) {
      sub_cor <- stats::cor(
        true_block[, valid_true, drop = FALSE], 
        est_block[, valid_est, drop = FALSE], 
        use = "pairwise.complete.obs"
      )
      cormat_block[valid_true, valid_est] <- sub_cor
    }
    cormat_block[!is.finite(cormat_block)] <- 0
    
    if (plot_heatmaps) {
      pheatmap(cormat_block, 
               cluster_rows = FALSE, cluster_cols = FALSE, 
               display_numbers = TRUE,
               main = paste0("Correlation: ", bname))
    }
    
    matched_in_block <- do.call(rbind, lapply(seq_len(nrow(matched)), function(i) {
      tcomp <- matched$true_comp[i]
      ecomp <- matched$est_comp[i]
      cor_val <- cormat_block[tcomp, ecomp]
      data.frame(true_comp = tcomp, est_comp = ecomp, cor_block = cor_val)
    }))
    
    per_block_summary[[bname]] <- list(cormat = cormat_block, matched = matched_in_block)
  }
  return(per_block_summary)
}

# ==============================================================================
# 5. EXECUTION & EVALUATION
# ==============================================================================

if (interactive()) {
  
  message("--- Step 2: Running Method ---")
  
  # --- Run Method ---
  # Note: Alpha is set low (0.01) to encourage dense loadings in simulation
  res_kurtosis <- mfa_weighted_fabia(
    simX_list, 
    gamma = 0.5,    
    p = 3, 
    alpha = 0.01,
    block_equalize = TRUE # This now triggers the lambda1 normalization
  )
  
  # --- Factor Recovery ---
  message("--- Step 3: Evaluating Factors ---")
  eval_factors <- evaluate_factor_recovery(sim_object, res_kurtosis, n_factors = 3, plot = TRUE, return_all = TRUE)
  print(eval_factors$summary)
  
  # --- Loading Recovery (Global) ---
  message("--- Step 4: Evaluating Loadings ---")
  eval_loadings <- evaluate_loading_recovery(sim_object, res_kurtosis$fabia, verbose = TRUE)
  
  # --- Per-Omic Loading Recovery ---
  message("--- Step 5: Per-Omic Correlations ---")
  per_block_res <- per_omic_loading_recovery(
    eval_res = eval_loadings, 
    sim_object = sim_object,
    plot_heatmaps = TRUE, 
    plot_scatter = FALSE
  )
}


# ==============================================================================
# 6. FEATURE SELECTION OVERLAP (JACCARD)
# ==============================================================================

evaluate_feature_jaccard <- function(sim_object, res, top_n = 300) {
  eval_load <- evaluate_loading_recovery(sim_object, res$fabia, verbose=FALSE)
  true_L <- eval_load$true_all
  est_L  <- eval_load$est_all_flipped
  matches <- eval_load$matched
  
  jaccards <- numeric(nrow(matches))
  
  for(i in 1:nrow(matches)) {
    t_idx <- matches$true_comp[i]
    e_idx <- matches$est_comp[i]
    
    top_true <- order(abs(true_L[, t_idx]), decreasing = TRUE)[1:top_n]
    top_est  <- order(abs(est_L[, e_idx]), decreasing = TRUE)[1:top_n]
    
    intersection <- length(intersect(top_true, top_est))
    union_set    <- length(union(top_true, top_est))
    jaccards[i] <- intersection / union_set
  }
  names(jaccards) <- paste0("Factor_", matches$true_comp)
  return(jaccards)
}

if(interactive()) {
  j_scores <- evaluate_feature_jaccard(sim_object, res_kurtosis, top_n = 300)
  barplot(j_scores, main="Feature Selection Overlap (Jaccard)", ylim=c(0,1))
}

# ==============================================================================
# 7. CLASSIFICATION ACCURACY (AUC-ROC)
# ==============================================================================

evaluate_sample_auc <- function(sim_object, res, factor_idx_matched) {
  true_factors <- get_true_factors_from_sim(sim_object)
  est_scores <- res$fabia$factor_scores
  aucs <- numeric()
  
  for(i in seq_along(factor_idx_matched)) {
    est_col <- factor_idx_matched[i]
    true_binary <- abs(true_factors[, i]) > 1e-3 
    
    if(length(unique(true_binary)) > 1) {
      roc_obj <- pROC::roc(true_binary, est_scores[, est_col], quiet=TRUE)
      aucs[i] <- roc_obj$auc
    } else {
      aucs[i] <- NA 
    }
  }
  return(aucs)
}

if(interactive()) {
  matched_indices <- eval_factors$summary$matched_idx 
  auc_scores <- evaluate_sample_auc(sim_object, res_kurtosis, matched_indices)
  print(data.frame(Factor = 1:length(auc_scores), AUC = auc_scores))
}


# ==============================================================================
# 8. PER-FACTOR SNR TREND ANALYSIS
# ==============================================================================

if(interactive()) {
  
  # 1. Define range: 0.01 to 2.0 with steps
  snr_values <- seq(from = 0.01, to = 2.0, by = 0.4)
  
  factor_results <- data.frame(
    SNR = numeric(), Omic = character(), Factor = character(),
    Best_Correlation = numeric(), stringsAsFactors = FALSE
  )
  
  message(sprintf("Analyzing Factor-Specific Recovery across %d SNR steps...", length(snr_values)))
  pb <- txtProgressBar(min = 0, max = length(snr_values), style = 3)
  
  for (i in seq_along(snr_values)) {
    s <- snr_values[i]
    setTxtProgressBar(pb, i)
    
    # A. Re-simulate
    sim1_s <- simulateMultiOmics(
      vector_features = c(4000, 2500), n_samples = 100, n_factors = 3, snr = s, 
      signal.samples = c(3, 1), signal.features = list(c(4.5, 0.5), c(4.5, 0.5)),
      factor_structure = "mixed", num.factor = "multiple", seed = 123
    )
    sim2_s <- simulateMultiOmics(
      vector_features = c(4000, 3500), n_samples = 100, n_factors = 3, snr = s,
      signal.samples = c(3, 1), signal.features = list(c(2.5, 0.5), c(3, 2.5)),
      factor_structure = "mixed", num.factor = "multiple", seed = 123
    )
    sim_obj_s <- sim1_s
    sim_obj_s$omics$omic2 <- sim2_s$omics$omic2
    sim_obj_s[["list_betas"]][[2]] <- sim2_s[["list_betas"]][[2]]
    simX_s <- sim_obj_s$omics
    names(simX_s) <- paste0("omic", seq_along(simX_s))
    
    # B. Run Method (Corrected Logic)
    res_s <- suppressMessages(
      mfa_weighted_fabia(simX_s, gamma = 0.5, p = 3, alpha = 0.01, seed = 1)
    )
    
    # C. Evaluate
    eval_load_s <- evaluate_loading_recovery(sim_obj_s, res_s$fabia, verbose = FALSE)
    per_omic_s <- per_omic_loading_recovery(eval_res = eval_load_s, sim_object = sim_obj_s, 
                                            plot_heatmaps = FALSE, plot_scatter = FALSE)
    
    # D. Extract Max
    extract_best_matches <- function(omic_res, omic_name, snr_val) {
      cormat <- omic_res$cormat
      if (!is.null(cormat) && nrow(cormat) > 0) {
        for (f_idx in 1:nrow(cormat)) {
          best_val <- max(abs(cormat[f_idx, ]), na.rm = TRUE)
          factor_results <<- rbind(factor_results, data.frame(
            SNR = snr_val, Omic = omic_name, 
            Factor = rownames(cormat)[f_idx], Best_Correlation = best_val
          ))
        }
      }
    }
    extract_best_matches(per_omic_s$omic1, "Omic1", s)
    extract_best_matches(per_omic_s$omic2, "Omic2", s)
  }
  close(pb)
  
  # Plot
  factor_results$Factor <- factor(factor_results$Factor, levels = c("factor1", "factor2", "factor3"))
  p_snr <- ggplot(factor_results, aes(x = SNR, y = Best_Correlation, color = Factor)) +
    facet_wrap(~Omic) + 
    geom_smooth(method = "loess", se = FALSE, alpha = 0.2, size = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("factor1" = "#E41A1C", "factor2" = "#377EB8", "factor3" = "#4DAF4A")) +
    scale_y_continuous(limits = c(0, 1.05)) +
    labs(title = "Factor Recovery vs SNR", y = "Max Correlation") +
    theme_minimal()
  
  print(p_snr)
}


# ==============================================================================
# 9. JACCARD INDEX TREND ANALYSIS 
# ==============================================================================

if(interactive()) {
  
  snr_values <- seq(from = 0.01, to = 2.0, by = 0.4)
  top_n_features <- 400 
  
  jaccard_results <- data.frame(
    SNR = numeric(), Omic = character(), Factor = character(),
    Jaccard_Index = numeric(), stringsAsFactors = FALSE
  )
  
  message(sprintf("Analyzing Jaccard Overlap (Top %d) across %d SNR steps...", 
                  top_n_features, length(snr_values)))
  pb <- txtProgressBar(min = 0, max = length(snr_values), style = 3)
  
  for (i in seq_along(snr_values)) {
    s <- snr_values[i]
    setTxtProgressBar(pb, i)
    
    # A. Re-simulate
    sim1_s <- simulateMultiOmics(vector_features = c(4000, 2500), n_samples = 100, n_factors = 3, snr = s, signal.samples = c(3, 1), signal.features = list(c(4.5, 0.5), c(4.5, 0.5)), factor_structure = "mixed", num.factor = "multiple", seed = 123)
    sim2_s <- simulateMultiOmics(vector_features = c(4000, 3500), n_samples = 100, n_factors = 3, snr = s, signal.samples = c(3, 1), signal.features = list(c(2.5, 0.5), c(3, 2.5)), factor_structure = "mixed", num.factor = "multiple", seed = 123)
    sim_obj_s <- sim1_s; sim_obj_s$omics$omic2 <- sim2_s$omics$omic2; sim_obj_s[["list_betas"]][[2]] <- sim2_s[["list_betas"]][[2]]
    simX_s <- sim_obj_s$omics; names(simX_s) <- paste0("omic", seq_along(simX_s))
    
    # B. Run Method
    res_s <- suppressMessages(mfa_weighted_fabia(simX_s, gamma = 0.5, p = 3, alpha = 0.01, seed = 1))
    
    # C. Jaccard Helper
    get_true_block_robust <- function(sim_obj, block_idx, n_global_factors = 3) {
      betas_list <- sim_obj$list_betas[[block_idx]]
      if (length(betas_list) == 0) return(NULL)
      n_features <- length(betas_list[[1]])
      mat_out <- matrix(0, nrow = n_features, ncol = n_global_factors)
      colnames(mat_out) <- paste0("factor", 1:n_global_factors)
      
      if (is.list(betas_list)) {
        for (nm in names(betas_list)) {
          idx <- if (grepl("beta", nm)) as.integer(sub("beta", "", nm)) else as.integer(nm)
          if (!is.na(idx) && idx <= n_global_factors) mat_out[, idx] <- as.numeric(betas_list[[nm]])
        }
      }
      return(mat_out)
    }
    
    process_jaccard <- function(omic_name, block_idx) {
      true_mat <- get_true_block_robust(sim_obj_s, block_idx, n_global_factors = 3)
      if (is.null(true_mat) || !omic_name %in% names(res_s$fabia$loadings_per_block)) return(NULL)
      est_mat <- res_s$fabia$loadings_per_block[[omic_name]]
      
      if (ncol(est_mat) > 0) {
        cormat <- cor(true_mat, est_mat, use="pairwise.complete.obs")
        cormat[is.na(cormat)] <- 0
        
        for (f in 1:ncol(true_mat)) {
          if (sd(true_mat[, f]) < 1e-6) {
            jaccard_results <<- rbind(jaccard_results, data.frame(SNR = s, Omic = omic_name, Factor = colnames(true_mat)[f], Jaccard_Index = 0))
            next
          }
          best_match_idx <- which.max(abs(cormat[f, ]))
          if (length(best_match_idx) > 0) {
            top_true <- order(abs(true_mat[, f]), decreasing = TRUE)[1:top_n_features]
            top_est  <- order(abs(est_mat[, best_match_idx]), decreasing = TRUE)[1:top_n_features]
            j_val <- length(intersect(top_true, top_est)) / length(union(top_true, top_est))
            jaccard_results <<- rbind(jaccard_results, data.frame(SNR = s, Omic = omic_name, Factor = colnames(true_mat)[f], Jaccard_Index = j_val))
          }
        }
      }
    }
    
    process_jaccard("omic1", 1)
    process_jaccard("omic2", 2)
  }
  close(pb)
  
  jaccard_results$Factor <- factor(jaccard_results$Factor, levels = c("factor1", "factor2", "factor3"))
  
  p_jac <- ggplot(jaccard_results, aes(x = SNR, y = Jaccard_Index, color = Factor)) +
    facet_wrap(~Omic) + 
    geom_smooth(method = "loess", se = FALSE, alpha = 0.2, size = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("factor1" = "#E41A1C", "factor2" = "#377EB8", "factor3" = "#4DAF4A")) +
    scale_y_continuous(limits = c(0, 1.0)) +
    labs(title = "Jaccard Index vs SNR", y = "Jaccard Index") +
    theme_minimal()
  
  print(p_jac)
}

