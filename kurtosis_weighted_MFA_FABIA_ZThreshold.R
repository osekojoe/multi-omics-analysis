# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================

# List of required packages
required_pkgs <- c("SUMO", "FactoMineR", "fabia", "clue", "pheatmap", 
                   "ggplot2", "matrixStats", "reshape2", "e1071", "pROC")

# Install missing packages
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install)) install.packages(to_install)

# Load libraries
library(SUMO)         # For simulateMultiOmics
library(FactoMineR)   # For MFA logic
library(fabia)        # Core factorization algorithm
library(clue)         # Hungarian algorithm (matching factors)
library(pheatmap)     # Visualization
library(ggplot2)      # Plotting
library(matrixStats)  # Efficient matrix operations
library(reshape2)     # Data reshaping
library(e1071)        # Kurtosis calculation
library(pROC)         # Classification metrics

# ==============================================================================
# 2. DATA SIMULATION
# ==============================================================================

message("--- Step 1: Simulating Data ---")

# Simulation 1: Structure A
sim1 <- simulateMultiOmics(
  vector_features = c(4000, 2500),
  n_samples = 100,
  n_factors = 3,
  snr = 0.5,
  signal.samples = c(3, 1),
  signal.features = list(c(4.5, 0.5), c(4.5, 0.5)),
  factor_structure = "mixed",
  num.factor = "multiple",
  seed = 123
)

# Simulation 2: Structure B (Different feature counts/signals)
sim2 <- simulateMultiOmics(
  vector_features = c(4000, 3500),
  n_samples = 100,
  n_factors = 3,
  snr = 0.5,
  signal.samples = c(3, 1),
  signal.features = list(c(2.5, 0.5), c(3, 2.5)),
  factor_structure = "mixed",
  num.factor = "multiple",
  seed = 123
)

# Combine omics from different simulations into one master object
sim_object <- sim1
sim_object$omics$omic2 <- sim2$omics$omic2
sim_object[["list_betas"]][[2]] <- sim2[["list_betas"]][[2]]

# Prepare final input list
simX_list <- sim_object$omics
names(simX_list) <- paste0("omic", seq_along(simX_list))

# --- PLOT: Raw Simulation Heatmap ---
plot_simData(sim_object = sim_object, type = "heatmap")

# ==============================================================================
# 3. CORE ALGORITHM: Kurtosis-Weighted MFA-FABIA
# ==============================================================================

#' Run FABIA with feature weighting based on Kurtosis and MFA Normalization
#'
#' @param Xlist List of matrices (omics data).
#' @param gamma Tempering parameter for weights (0 = no weighting).
#' @param p FABIA parameter (number of biclusters).
#' @param alpha FABIA parameter (sparseness). Lower values = denser loadings.
#' @param scale_vars Boolean. Scale variables?
#' @param block_equalize Boolean. Equalize block variance?
#' @param block_equalize_method Method for equalization ("trace", "frobenius", "none").
#' @param cap_weights_quantile Quantile to cap extreme weights (default 0.99).
#' @param weight_floor Minimum weight to prevent zeroing out features.
#' @param verbose Print diagnostic messages?
#' @param seed Random seed.
#' @return List containing weighted data, FABIA results, and weights.
mfa_weighted_fabia <- function(Xlist,
                               gamma = 0.5,
                               p = 4,
                               alpha = 0.5,
                               scale_vars = TRUE,
                               block_equalize = TRUE,
                               block_equalize_method = c("trace", "frobenius", "none"),
                               cap_weights_quantile = 0.99,
                               weight_floor = 1e-3,     # NEW: prevents zeroing out
                               verbose = TRUE,          # NEW: diagnostics
                               seed = 1) {
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  stopifnot(is.list(Xlist), length(Xlist) >= 2)
  set.seed(seed)
  
  block_equalize_method <- match.arg(block_equalize_method)
  
  clean_and_prep <- function(M) {
    M <- as.matrix(M)
    M[is.infinite(M)] <- NA_real_
    vars <- matrixStats::colVars(M, na.rm = TRUE)
    keep <- is.finite(vars) & vars > 1e-9
    M <- M[, keep, drop = FALSE]
    if (ncol(M) == 0) stop("A block has no non-constant finite-variance columns after cleaning.")
    M
  }
  
  Xlist_clean <- lapply(Xlist, clean_and_prep)
  
  # B) scaling + optional block equalization
  Zlist <- lapply(Xlist_clean, function(X) {
    
    if (scale_vars) {
      s <- scale(X, center = TRUE, scale = TRUE)
      s[is.na(s)] <- 0
    } else {
      s <- scale(X, center = TRUE, scale = FALSE)
      s[is.na(s)] <- 0
    }
    
    if (block_equalize && block_equalize_method != "none") {
      if (block_equalize_method == "trace") {
        tr <- sum(matrixStats::colVars(s), na.rm = TRUE)
        if (is.finite(tr) && tr > 1e-12) s <- s / sqrt(tr)
      } else if (block_equalize_method == "frobenius") {
        fn <- sqrt(sum(s^2, na.rm = TRUE))
        if (is.finite(fn) && fn > 1e-12) s <- s / fn
      }
    }
    s
  })
  
  names(Zlist) <- names(Zlist) %||% paste0("omic", seq_along(Zlist))
  for (n in names(Zlist)) colnames(Zlist[[n]]) <- paste0(n, "_", colnames(Zlist[[n]]))
  
  Z <- do.call(cbind, Zlist)
  grp <- vapply(Zlist, ncol, integer(1))
  
  # C) weights
  if (verbose) message("Calculating block-wise kurtosis for feature weighting...")
  
  safe_kurtosis <- function(M) {
    apply(M, 2, function(v) {
      tryCatch(e1071::kurtosis(v, type = 2, na.rm = TRUE), error = function(e) NA_real_)
    })
  }
  
  kurtosis_per_block <- list()
  weights_per_block  <- list()
  
  for (bn in names(Zlist)) {
    M <- Zlist[[bn]]
    k <- safe_kurtosis(M)
    
    # --- CRITICAL FIX: if gamma == 0 => true no-weighting (all ones)
    if (gamma == 0) {
      w <- rep(1, length(k))
      names(w) <- colnames(M)
      
      kurtosis_per_block[[bn]] <- k
      weights_per_block[[bn]]  <- w
      next
    }
    
    # More stable mapping: use positive kurtosis but with a floor so nothing becomes 0
    # (If you prefer abs(k), replace pmax(0,k) with abs(k))
    w_raw <- pmax(0, k)
    w_raw[!is.finite(w_raw)] <- 0
    w_raw <- pmax(weight_floor, w_raw)  # floor prevents annihilation
    
    mean_w <- mean(w_raw, na.rm = TRUE)
    if (!is.finite(mean_w) || mean_w <= 0) mean_w <- 1
    
    w <- (w_raw / mean_w)^gamma
    w[!is.finite(w)] <- 1
    
    if (is.finite(cap_weights_quantile) && cap_weights_quantile < 1) {
      cap <- suppressWarnings(stats::quantile(w, probs = cap_weights_quantile, na.rm = TRUE, names = FALSE))
      if (is.finite(cap)) w <- pmin(w, cap)
    }
    
    # renormalize to mean 1 within block
    mw2 <- mean(w, na.rm = TRUE)
    if (is.finite(mw2) && mw2 > 0) w <- w / mw2
    
    names(w) <- colnames(M)
    kurtosis_per_block[[bn]] <- k
    weights_per_block[[bn]]  <- w
  }
  
  w_all <- unlist(weights_per_block, use.names = TRUE)
  w_all <- w_all[colnames(Z)]
  w_all[is.na(w_all)] <- 1  # if anything is missing, do not zero it out
  
  Z_w <- sweep(Z, 2, w_all, `*`)
  
  # Diagnostics: detect the “empty/constant” case before FABIA
  if (verbose) {
    col_sds <- matrixStats::colSds(Z_w, na.rm = TRUE)
    msg <- sprintf("Pre-FABIA Z_w: dim=%dx%d; mean(colSD)=%.3g; min(colSD)=%.3g; frac(|x|<1e-12)=%.2f%%",
                   nrow(Z_w), ncol(Z_w),
                   mean(col_sds), min(col_sds),
                   100 * mean(abs(Z_w) < 1e-12))
    message(msg)
  }
  
  if (ncol(Z_w) == 0) stop("Z_w has 0 columns after preprocessing.")
  if (sd(as.vector(Z_w), na.rm = TRUE) < 1e-12) stop("Z_w is (near) constant after preprocessing/weighting.")
  
  cuts <- cumsum(grp)
  weighted_blocks <- lapply(seq_along(grp), function(i) {
    if (i == 1) Z_w[, 1:cuts[1], drop = FALSE]
    else Z_w[, (cuts[i - 1] + 1):cuts[i], drop = FALSE]
  })
  names(weighted_blocks) <- names(Zlist)
  
  # D) FABIA
  fab_fit <- fabia::fabia(as.matrix(t(Z_w)), p = p, alpha = alpha, random = 1)
  
  factor_scores <- t(fab_fit@Z)
  loadings_all  <- fab_fit@L
  
  loadings_per_block <- list()
  for (bn in names(Zlist)) {
    pattern <- paste0("^", bn, "_")
    idx <- grep(pattern, rownames(loadings_all))
    if (length(idx) > 0) loadings_per_block[[bn]] <- loadings_all[idx, , drop = FALSE]
  }
  
  fab_res <- list(
    FABIA              = fab_fit,
    factor_scores      = factor_scores,
    loadings_per_block = loadings_per_block,
    group_sizes        = grp
  )
  
  list(
    weighted              = Z_w,
    weighted_blocks       = weighted_blocks,
    fabia                 = fab_res,
    weights               = w_all,
    weights_per_block     = weights_per_block,
    kurtosis_per_block    = kurtosis_per_block,
    block_equalize        = block_equalize,
    block_equalize_method = block_equalize_method,
    gamma                 = gamma
  )
}
# ==============================================================================
# 4. EVALUATION HELPERS
# ==============================================================================

#' Retrieve True Factor Matrix from Simulation Object
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

#' Evaluate Factor Score Recovery (Hungarian Matching)
evaluate_factor_recovery <- function(sim_object, res, n_factors = NULL, plot = TRUE, return_all = FALSE) {
  
  if (is.list(res) && !is.null(res$fabia) && !is.null(res$fabia$factor_scores)) {
    est <- as.matrix(res$fabia$factor_scores)
  } else if (!is.null(res$factor_scores)) {
    est <- as.matrix(res$factor_scores)
  } else {
    stop("Cannot find estimated factor scores.")
  }
  
  true <- as.matrix(get_true_factors_from_sim(sim_object))
  
  K_est <- ncol(est); K_true <- ncol(true)
  if (is.null(n_factors)) n_factors <- min(K_est, K_true)
  
  estK <- est[, seq_len(min(ncol(est), n_factors)), drop = FALSE]
  trueK <- true[, seq_len(min(ncol(true), n_factors)), drop = FALSE]
  K <- ncol(trueK)
  
  # Compute pairwise correlation
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  cor_mat_na <- cor_mat; cor_mat_na[is.na(cor_mat_na)] <- 0
  
  # Hungarian Assignment (Maximize diagonal correlation)
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
  
  # Metrics (R2 and MSE)
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

#' Evaluate Loading Matrix Recovery (Global)
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
    
    # Logic to parse beta lists from SUMO simulation
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


#' Evaluate Factor Score Recovery using Greedy Matching
#' 
#' Matches factors by iteratively picking the highest correlation in the matrix,
#' assigning that pair, and removing them from further consideration.
evaluate_factor_recovery_greedy <- function(sim_object, res, n_factors = NULL, plot = TRUE) {
  
  # 1. Extract Scores
  if (is.list(res) && !is.null(res$fabia) && !is.null(res$fabia$factor_scores)) {
    est <- as.matrix(res$fabia$factor_scores)
  } else if (!is.null(res$factor_scores)) {
    est <- as.matrix(res$factor_scores)
  } else {
    stop("Cannot find estimated factor scores.")
  }
  
  true <- as.matrix(get_true_factors_from_sim(sim_object))
  
  # 2. Dimensions
  K_est <- ncol(est); K_true <- ncol(true)
  if (is.null(n_factors)) n_factors <- min(K_est, K_true)
  
  estK <- est[, seq_len(min(ncol(est), n_factors)), drop = FALSE]
  trueK <- true[, seq_len(min(ncol(true), n_factors)), drop = FALSE]
  K <- ncol(trueK)
  
  # 3. Compute Correlation Matrix
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  cor_mat_abs <- abs(cor_mat)
  cor_mat_abs[is.na(cor_mat_abs)] <- 0
  
  # 4. GREEDY MATCHING LOGIC
  assigned_cols <- rep(NA, K) # To store which Est column matches True i
  
  # We work on a copy of the matrix so we can "cross out" rows/cols
  temp_mat <- cor_mat_abs
  
  # Iterate K times to find K pairs
  for(step in 1:K) {
    # Find position of the global maximum in the remaining matrix
    max_idx <- which(temp_mat == max(temp_mat), arr.ind = TRUE)
    
    # If multiple maxes (rare), take the first one
    best_true_idx <- max_idx[1, 1]
    best_est_idx  <- max_idx[1, 2]
    
  }
  
  # --- Greedy Loop ---
  pairs_found <- list()
  temp_mat <- cor_mat_abs
  
  for(i in 1:K) {
    # Find max value coordinates
    loc <- which(temp_mat == max(temp_mat), arr.ind = TRUE)
    r <- loc[1, 1] # True Factor Index
    c <- loc[1, 2] # Est Factor Index
    
    pairs_found[[i]] <- c(true=r, est=c)
    
    # "Cross out" this row and column by setting to -1 (so max won't pick them again)
    temp_mat[r, ] <- -1
    temp_mat[, c] <- -1
  }
  
  # Convert pairs list to the `assigned_cols` vector format expected by plotting code
  # assigned_cols[i] should be the Estimated index that matches True Factor i
  assigned_cols <- integer(K)
  for(p in pairs_found) {
    assigned_cols[p["true"]] <- p["est"]
  }
  
  # 5. Metrics & Alignment (Same as original)
  matched_corrs <- numeric(K); matched_signs <- numeric(K)
  for (i in seq_len(K)) {
    j <- assigned_cols[i]
    cval <- if (is.na(cor_mat[i, j])) 0 else cor_mat[i, j]
    matched_corrs[i] <- cval
    matched_signs[i] <- if (sign(cval) == 0) 1 else sign(cval)
  }
  
  est_matched_signed <- estK[, assigned_cols, drop = FALSE] %*% diag(matched_signs)
  colnames(est_matched_signed) <- paste0("Est_matched_", seq_len(K))
  
  # R2 / MSE Calculation
  r2 <- numeric(K); mse <- numeric(K)
  for (i in seq_len(K)) {
    fit <- lm(trueK[, i] ~ est_matched_signed[, i])
    ssr <- sum((predict(fit) - mean(trueK[, i]))^2)
    sst <- sum((trueK[, i] - mean(trueK[, i]))^2)
    r2[i] <- if (sst == 0) NA else ssr / sst
    mse[i] <- mean((trueK[, i] - est_matched_signed[, i])^2)
  }
  
  if (plot) {
    pheatmap(cor_mat, main = "Correlation (Greedy Matching)",
             cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
  }
  
  summary_tbl <- data.frame(
    true_factor = colnames(trueK), matched_idx = assigned_cols,
    corr = matched_corrs, R2 = r2, MSE = mse, stringsAsFactors = FALSE
  )
  
  return(list(cor_mat = cor_mat, summary = summary_tbl, mean_abs_cor = mean(abs(matched_corrs))))
}

# ==============================================================================
# 5. JACCARD ANALYSIS FUNCTIONS (FEATURE OVERLAP)
# ==============================================================================

#' Helper: Calculate Z-score Boolean (Outlier Detection)
#' 
get_outliers_zscore <- function(x, z_thresh = 3.0) {
  if(sd(x, na.rm = TRUE) < 1e-9) return(rep(FALSE, length(x)))
  z_scores <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  return(abs(z_scores) > z_thresh)
}

#' Evaluate Feature Recovery (Jaccard) Per Omic Layer
evaluate_jaccard_per_omic <- function(sim_object, res, z_thresh = 3.0) {
  
  # 1. Global Factor Matching
  eval_global <- evaluate_loading_recovery(sim_object, res$fabia, verbose=FALSE)
  matches <- eval_global$matched
  true_global <- eval_global$true_all
  est_global  <- eval_global$est_all_flipped
  
  # 2. Define Block Indices
  block_sizes <- eval_global$block_sizes
  block_names <- names(block_sizes)
  cuts <- cumsum(block_sizes)
  
  idx_list <- lapply(seq_along(block_sizes), function(i) {
    if (i == 1) seq_len(cuts[1]) else (cuts[i - 1] + 1):cuts[i]
  })
  names(idx_list) <- block_names
  
  # 3. Iterate Per Omic -> Per Factor
  results_df <- data.frame()
  
  for(bname in block_names) {
    rows <- idx_list[[bname]]
    true_block <- true_global[rows, , drop=FALSE]
    est_block  <- est_global[rows, , drop=FALSE]
    
    for(i in 1:nrow(matches)) {
      t_col <- matches$true_comp[i]
      e_col <- matches$est_comp[i]
      
      v_true <- true_block[, t_col]
      v_est  <- est_block[, e_col]
      
      bool_true <- get_outliers_zscore(v_true, z_thresh)
      bool_est  <- get_outliers_zscore(v_est,  z_thresh)
      
      # 
      if(sum(bool_true) == 0) {
        j_val <- 0
      } else {
        intersect_len <- sum(bool_true & bool_est)
        union_len     <- sum(bool_true | bool_est)
        j_val <- if(union_len > 0) intersect_len / union_len else 0
      }
      
      results_df <- rbind(results_df, data.frame(
        Omic = bname,
        Factor = paste0("Factor ", t_col),
        Jaccard = j_val,
        Threshold = z_thresh,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(results_df)
}

# ==============================================================================
# 6. SENSITIVITY ANALYSIS MODULES (Refined for Production)
# ==============================================================================

#' Run SNR Sensitivity Analysis (Trend)
run_snr_sensitivity_analysis <- function() {
  snr_values <- seq(from = 0.01, to = 2.0, by = 0.1)
  factor_results <- data.frame()
  
  message(sprintf("Analyzing Factor Recovery across %d SNR steps...", length(snr_values)))
  pb <- txtProgressBar(min = 0, max = length(snr_values), style = 3)
  
  for (i in seq_along(snr_values)) {
    s <- snr_values[i]
    setTxtProgressBar(pb, i)
    
    # 1. Re-simulate
    sim1_s <- simulateMultiOmics(vector_features = c(4000, 2500), n_samples = 100, n_factors = 3, snr = s, signal.samples = c(3, 1), signal.features = list(c(4.5, 0.5), c(4.5, 0.5)), factor_structure = "mixed", num.factor = "multiple", seed = 123)
    sim2_s <- simulateMultiOmics(vector_features = c(4000, 3500), n_samples = 100, n_factors = 3, snr = s, signal.samples = c(3, 1), signal.features = list(c(2.5, 0.5), c(3, 2.5)), factor_structure = "mixed", num.factor = "multiple", seed = 123)
    sim_obj_s <- sim1_s
    sim_obj_s$omics$omic2 <- sim2_s$omics$omic2
    sim_obj_s[["list_betas"]][[2]] <- sim2_s[["list_betas"]][[2]]
    simX_s <- sim_obj_s$omics; names(simX_s) <- paste0("omic", seq_along(simX_s))
    
    # 2. Run Algorithm (Updated to use mfa_weighted_fabia)
    # Wrap in tryCatch to prevent loop crash on singular matrices at low SNR
    res_s <- tryCatch({
      suppressMessages(mfa_weighted_fabia(simX_s, gamma = 0.5, p = 3, alpha = 0.5, seed = 1, verbose = FALSE))
    }, error = function(e) return(NULL))
    
    if(is.null(res_s)) next
    
    # 3. Evaluate Per Omic
    eval_load_s <- evaluate_loading_recovery(sim_obj_s, res_s$fabia, verbose = FALSE)
    per_omic_s <- per_omic_loading_recovery(eval_res = eval_load_s, sim_object = sim_obj_s, plot_heatmaps = FALSE, plot_scatter = FALSE)
    
    # 4. Extract Max Correlation
    extract_best <- function(omic_res, omic_name) {
      if(!is.null(omic_res$cormat) && nrow(omic_res$cormat) > 0) {
        for(f_idx in 1:nrow(omic_res$cormat)) {
          best_val <- max(abs(omic_res$cormat[f_idx, ]), na.rm = TRUE)
          factor_results <<- rbind(factor_results, data.frame(SNR = s, Omic = omic_name, Factor = rownames(omic_res$cormat)[f_idx], Best_Correlation = best_val))
        }
      }
    }
    extract_best(per_omic_s$omic1, "Omic1")
    extract_best(per_omic_s$omic2, "Omic2")
  }
  close(pb)
  
  if(nrow(factor_results) == 0) {
    message("No results generated for SNR Sensitivity.")
    return(NULL)
  }
  
  # Plot
  factor_results$Factor <- factor(factor_results$Factor, levels = c("factor1", "factor2", "factor3"))
  p <- ggplot(factor_results, aes(x = SNR, y = Best_Correlation, color = Factor)) +
    facet_wrap(~Omic) + 
    geom_smooth(method = "loess", se = FALSE, alpha = 0.2, size = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("factor1" = "#E41A1C", "factor2" = "#377EB8", "factor3" = "#4DAF4A")) +
    scale_y_continuous(limits = c(0, 1.05)) +
    labs(title = "Factor Recovery vs SNR", y = "Max Correlation") + theme_minimal()
  
  print(p)
}

#' Run Jaccard Trend Analysis (Z-Score Method)

run_jaccard_trend_analysis <- function(z_threshold = 3.0) {
  snr_values <- seq(from = 0.01, to = 2.0, by = 0.1)
  
  # Initialize with columns to prevent ggplot error if empty
  jaccard_results <- data.frame(
    SNR = numeric(), 
    Omic = character(), 
    Factor = character(), 
    Jaccard_Index = numeric(), 
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Analyzing Jaccard Overlap (Z > %.1f) across %d SNR steps...", z_threshold, length(snr_values)))
  pb <- txtProgressBar(min = 0, max = length(snr_values), style = 3)
  
  for (i in seq_along(snr_values)) {
    s <- snr_values[i]
    setTxtProgressBar(pb, i)
    
    # 1. Re-simulate
    sim1_s <- simulateMultiOmics(vector_features = c(4000, 2500), n_samples = 100, n_factors = 3, snr = s, signal.samples = c(3, 1), signal.features = list(c(4.5, 0.5), c(4.5, 0.5)), factor_structure = "mixed", num.factor = "multiple", seed = 123)
    sim2_s <- simulateMultiOmics(vector_features = c(4000, 3500), n_samples = 100, n_factors = 3, snr = s, signal.samples = c(3, 1), signal.features = list(c(2.5, 0.5), c(3, 2.5)), factor_structure = "mixed", num.factor = "multiple", seed = 123)
    sim_obj_s <- sim1_s; sim_obj_s$omics$omic2 <- sim2_s$omics$omic2; sim_obj_s[["list_betas"]][[2]] <- sim2_s[["list_betas"]][[2]]
    simX_s <- sim_obj_s$omics; names(simX_s) <- paste0("omic", seq_along(simX_s))
    
    # 2. Run Method (using mfa_weighted_fabia)
    res_s <- tryCatch({
      suppressMessages(mfa_weighted_fabia(simX_s, gamma = 0.5, p = 3, alpha = 0.5, seed = 1, verbose = FALSE))
    }, error = function(e) return(NULL))
    
    if(is.null(res_s)) next
    
    # 3. Helpers for this specific loop
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
    
    # 4. Calculate Jaccard per Omic
    process_jaccard <- function(omic_name, block_idx) {
      true_mat <- get_true_block_robust(sim_obj_s, block_idx, n_global_factors = 3)
      if (is.null(true_mat) || !omic_name %in% names(res_s$fabia$loadings_per_block)) return(NULL)
      est_mat <- res_s$fabia$loadings_per_block[[omic_name]]
      
      if (ncol(est_mat) > 0) {
        
        # --- FIX: Handle Zero Variance Before Correlation ---
        # Identify columns with non-zero variance
        valid_true <- apply(true_mat, 2, sd, na.rm=TRUE) > 1e-12
        valid_est  <- apply(est_mat, 2, sd, na.rm=TRUE) > 1e-12
        
        # Initialize correlation matrix with zeros
        cormat <- matrix(0, nrow = ncol(true_mat), ncol = ncol(est_mat))
        
        # Only compute correlation for valid pairs
        if (any(valid_true) && any(valid_est)) {
          # subset matrices to valid columns only
          t_sub <- true_mat[, valid_true, drop=FALSE]
          e_sub <- est_mat[, valid_est, drop=FALSE]
          
          # compute
          sub_cor <- cor(t_sub, e_sub, use="pairwise.complete.obs")
          sub_cor[is.na(sub_cor)] <- 0
          
          # map back to full matrix
          cormat[valid_true, valid_est] <- sub_cor
        }
        # ----------------------------------------------------
        
        for (f in 1:ncol(true_mat)) {
          fname <- colnames(true_mat)[f]
          v_true <- true_mat[, f]
          bool_true <- get_outliers_zscore(v_true, z_threshold)
          
          if (sum(bool_true) == 0) {
            jaccard_results <<- rbind(jaccard_results, data.frame(SNR = s, Omic = omic_name, Factor = fname, Jaccard_Index = 0))
            next
          }
          
          best_match_idx <- which.max(abs(cormat[f, ]))
          
          # Ensure we actually have a match (and the correlation isn't just 0 everywhere)
          if (length(best_match_idx) > 0 && max(abs(cormat[f, ])) > 0) {
            v_est <- est_mat[, best_match_idx]
            bool_est <- get_outliers_zscore(v_est, z_threshold)
            intersect_len <- sum(bool_true & bool_est)
            union_len     <- sum(bool_true | bool_est)
            j_val <- if(union_len > 0) intersect_len / union_len else 0
            
            jaccard_results <<- rbind(jaccard_results, data.frame(SNR = s, Omic = omic_name, Factor = fname, Jaccard_Index = j_val))
          } else {
            # If no correlation found (e.g. estimate was constant), Jaccard is 0
            jaccard_results <<- rbind(jaccard_results, data.frame(SNR = s, Omic = omic_name, Factor = fname, Jaccard_Index = 0))
          }
        }
      }
    }
    
    process_jaccard("omic1", 1)
    process_jaccard("omic2", 2)
  }
  
  close(pb)
  
  # Check if dataframe is populated
  if(nrow(jaccard_results) == 0) {
    message("\nWarning: No Jaccard results generated. Check method convergence.")
    return(NULL)
  }
  
  # Plot
  jaccard_results$Factor <- factor(jaccard_results$Factor, levels = c("factor1", "factor2", "factor3"))
  p <- ggplot(jaccard_results, aes(x = SNR, y = Jaccard_Index, color = Factor)) +
    facet_wrap(~Omic) + 
    geom_smooth(method = "loess", se = FALSE, alpha = 0.2, size = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("factor1" = "#E41A1C", "factor2" = "#377EB8", "factor3" = "#4DAF4A")) +
    scale_y_continuous(limits = c(0, 1.05)) +
    labs(title = paste0("Feature Selection Overlap (Z > ", z_threshold, ")"), x = "SNR", y = "Jaccard Index") +
    theme_minimal()
  
  print(p)
}

# ==============================================================================
# 7. MAIN EXECUTION BLOCK
# ==============================================================================

if (interactive()) {
  message("\n--- Step 1: Running Core Method ---")
  # Use the fixed function
  res_kurtosis <- mfa_weighted_fabia(simX_list, gamma = 0.5, p = 3, alpha = 0.5, scale_vars=TRUE, block_equalize = TRUE, block_equalize_method="trace", seed=1)
  
  # --- Step 2: Compare Matching Algorithms ---
  
  message("\n--- Evaluation A: Optimal Matching (Hungarian) ---")
  eval_opt <- evaluate_factor_recovery(sim_object, res_kurtosis, n_factors = 3, plot = FALSE)
  print(eval_opt$summary)
  
  message("\n--- Evaluation B: Greedy Matching ---")
  eval_greedy <- evaluate_factor_recovery_greedy(sim_object, res_kurtosis, n_factors = 3, plot = TRUE)
  print(eval_greedy$summary)
  
  # Check if they differ
  if(identical(eval_opt$summary$matched_idx, eval_greedy$summary$matched_idx)) {
    message("\nResult: Both algorithms found the same matching.")
  } else {
    message("\nResult: Algorithms found DIFFERENT matchings.")
  }
  
  message("\n--- Step 3: Evaluating Per-Omic Loadings ---")
  eval_loadings <- evaluate_loading_recovery(sim_object, res_kurtosis$fabia, verbose = TRUE)
  per_block_res <- per_omic_loading_recovery(eval_res = eval_loadings, sim_object = sim_object, plot_heatmaps = TRUE)
  
  message("\n--- Step 4: Evaluating Feature Selection (Jaccard) ---")
  jaccard_df <- evaluate_jaccard_per_omic(sim_object, res_kurtosis, z_thresh = 2.0)
  print(jaccard_df)
  
  message("\n--- Step 5: Running SNR Sensitivity Analysis ---")
  run_snr_sensitivity_analysis()
  
  message("\n--- Step 6: Running Jaccard Trend Analysis ---")
  run_jaccard_trend_analysis(z_threshold = 1.5)
}

