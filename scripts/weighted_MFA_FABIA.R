# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================

# ==============================================================================

# Ensure required packages are installed
required_pkgs <- c("SUMO", "FactoMineR", "fabia", "clue", "pheatmap", 
                   "ggplot2", "matrixStats", "reshape2")
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


# ==============================================================================
# 2. CORE ALGORITHM: MFA-Weighted FABIA
# ==============================================================================

#' Run FABIA with feature weighting derived from MFA
#'
#' @param Xlist List of matrices (omics data).
#' @param L Number of MFA components to use for weights (0 = no weighting).
#' @param gamma Tempering parameter for weights.
#' @param p FABIA parameter (p-percentile).
#' @param alpha FABIA parameter (sparseness).
#' @param block_equalize Boolean. Equalize block variance before MFA?
#' @param scale_vars Boolean. Scale variables?
#' @param seed Random seed.
#' @return List containing MFA results, weighted data, FABIA results, and weights.
#' Run FABIA with feature weighting derived from MFA
mfa_weighted_fabia <- function(Xlist,
                               L = 2,
                               gamma = 0.5,
                               p = 3, 
                               alpha = 0.2,
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
  
  # --- 2) MFA ---
  df_all <- as.data.frame(Z)
  mfa <- FactoMineR::MFA(
    df_all,
    group      = grp,
    type       = rep("s", length(grp)),
    ncp        = max(L, 5),
    name.group = names(Zlist),
    graph      = FALSE
  )
  
  # --- 3) Build feature weights ---
  if (L > 0) {
    cos2 <- try(mfa$quanti.var$cos2, silent = TRUE)
    if (inherits(cos2, "try-error") || is.null(cos2)) {
      coords <- mfa$quanti.var$coord[, seq_len(L), drop = FALSE]
      w_raw  <- rowMeans(coords^2)
    } else {
      w_raw  <- rowMeans(cos2[, seq_len(L), drop = FALSE])
    }
    w <- w_raw[colnames(Z)]
    w <- w / mean(w, na.rm = TRUE)
    w[!is.finite(w)] <- 1
    w <- w^gamma
  } else {
    w <- rep(1, ncol(Z))
    names(w) <- colnames(Z)
  }
  
  # --- 4) Apply weights ---
  Z_w <- sweep(Z, 2, w, `*`)    
  
  cuts <- cumsum(grp)
  weighted_blocks <- lapply(seq_along(grp), function(i) {
    if (i == 1) Z_w[, 1:cuts[1], drop = FALSE]
    else Z_w[, (cuts[i - 1] + 1):cuts[i], drop = FALSE]
  })
  names(weighted_blocks) <- names(Zlist)
  
  # --- 5) FABIA on weighted data ---
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
    weights         = w
  )
}

# ==============================================================================
# 3. EVALUATION HELPERS
# ==============================================================================

# --- Helper to get true factors ---
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

# --- 3a. Factor Score Recovery ---
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
    # 1. Heatmap
    pheatmap(cor_mat, main = "Corr. score (true x estimated factors)",
             cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
    
    # 2. Scatter plots
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

# --- 3b. Loading Recovery (Global) ---
evaluate_loading_recovery <- function(sim_object, fab_res, verbose = TRUE) {
  K <- length(sim_object$list_alphas)
  block_names <- names(sim_object$omics)
  block_sizes <- vapply(sim_object$omics, ncol, integer(1))
  p_total <- sum(block_sizes)
  
  # Build true_all loadings matrix
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
  
  # Match and flip
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

# --- 3c. Per-Omic Loading Recovery (Visualizations) ---
#' Per-Omic Loading Recovery with Zero-SD Protection
per_omic_loading_recovery <- function(eval_res, sim_object, plot_heatmaps = TRUE, plot_scatter = TRUE, scatter_n = NULL) {
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
    
    # --- ROBUST CORRELATION CALCULATION ---
    # 1. Identify columns with Zero SD (Constant vectors)
    #    We use a small epsilon for floating point safety
    sd_true <- apply(true_block, 2, sd, na.rm = TRUE)
    sd_est  <- apply(est_block, 2, sd, na.rm = TRUE)
    
    valid_true <- sd_true > 1e-12
    valid_est  <- sd_est > 1e-12
    
    # 2. Initialize Correlation Matrix with Zeros
    #    (If SD is zero, correlation is mathematically 0 for our purposes)
    cormat_block <- matrix(0, nrow = ncol(true_block), ncol = ncol(est_block))
    rownames(cormat_block) <- colnames(true_block)
    colnames(cormat_block) <- colnames(est_block)
    
    # 3. Only calculate correlation for valid vectors
    if (any(valid_true) && any(valid_est)) {
      # Subset only valid columns
      sub_cor <- stats::cor(
        true_block[, valid_true, drop = FALSE], 
        est_block[, valid_est, drop = FALSE], 
        use = "pairwise.complete.obs"
      )
      # Fill the valid slots in the full matrix
      cormat_block[valid_true, valid_est] <- sub_cor
    }
    
    # 4. Clean up any remaining NAs (rare, but good practice)
    cormat_block[!is.finite(cormat_block)] <- 0
    # ---------------------------------------
    
    # 1. Plot Heatmap
    if (plot_heatmaps) {
      cm <- cormat_block
      # Clean labels for plot
      rownames(cm) <- paste0("T", seq_len(nrow(cm)))
      colnames(cm) <- paste0("E", seq_len(ncol(cm)))
      
      op <- par(no.readonly = TRUE)
      par(mfrow = c(1,1), mar = c(5,5,4,2))
      
      # Heatmap with conditional color scaling
      image(1:ncol(cm), 1:nrow(cm), t(cm[nrow(cm):1, , drop = FALSE]),
            axes = FALSE, xlab = "Estimated", ylab = "True",
            main = paste0("Block: ", bname, " (Correlation)"))
      axis(1, at = 1:ncol(cm), labels = colnames(cm), las = 2)
      axis(2, at = 1:nrow(cm), labels = rev(rownames(cm)), las = 2)
      
      # Text labels
      for (i in 1:nrow(cm)) {
        for (j in 1:ncol(cm)) {
          text(j, nrow(cm)-i+1, sprintf("%.2f", cm[i,j]), cex = 0.7)
        }
      }
      par(op)
    }
    
    # 2. Plot Scatter
    matched_in_block <- do.call(rbind, lapply(seq_len(nrow(matched)), function(i) {
      tcomp <- matched$true_comp[i]; ecomp <- matched$est_comp[i]
      
      # Use our pre-calculated robust matrix to get the value
      cor_val <- cormat_block[tcomp, ecomp]
      
      data.frame(true_comp = tcomp, est_comp = ecomp, cor_block = cor_val)
    }))
    
    if (plot_scatter) {
      to_plot <- matched_in_block
      if (!is.null(scatter_n)) to_plot <- to_plot[seq_len(min(nrow(to_plot), scatter_n)), ]
      
      for (ri in seq_len(nrow(to_plot))) {
        tc <- to_plot$true_comp[ri]; ec <- to_plot$est_comp[ri]
        
        # Check if vectors are constant (all zeros) to add a warning in the plot title
        is_const_true <- !valid_true[tc]
        is_const_est  <- !valid_est[ec]
        
        df <- data.frame(true = true_block[, tc], est = est_block[, ec])
        
        # Determine subtitle based on zero-variance
        sub_txt <- sprintf("r = %.2f", to_plot$cor_block[ri])
        if (is_const_true) sub_txt <- paste(sub_txt, "(True vector is constant/zero)")
        if (is_const_est) sub_txt <- paste(sub_txt, "(Est vector is constant/zero)")
        
        p <- ggplot(df, aes(x = true, y = est)) +
          geom_point(alpha = 0.4) + geom_smooth(method = "lm", se = FALSE) +
          labs(title = sprintf("%s: True %d vs Est %d", bname, tc, ec),
               subtitle = sub_txt)
        print(p)
      }
    }
    per_block_summary[[bname]] <- list(cormat = cormat_block, matched = matched_in_block)
  }
  return(per_block_summary)
}
# ==============================================================================
# 4. EXECUTION
# ==============================================================================

if (interactive()) {
  # --- A. Data Simulation ---
  sim1 <- simulateMultiOmics(
    vector_features = c(4000, 2500), n_samples = 100, n_factors = 3, snr = 2.0,
    signal.samples = c(3, 1), signal.features = list(c(4.5, 0.5), c(4.5, 0.5)),
    factor_structure = "mixed", num.factor = "multiple", seed = 123
  )
  sim2 <- simulateMultiOmics(
    vector_features = c(4000, 3500), n_samples = 100, n_factors = 3, snr = 1.0,
    signal.samples = c(3, 1), signal.features = list(c(2.5, 0.5), c(3, 2.5)),
    factor_structure = "mixed", num.factor = "multiple", seed = 123
  )
  sim_object <- sim1
  sim_object$omics$omic2 <- sim2$omics$omic2
  sim_object[["list_betas"]][[2]] <- sim2[["list_betas"]][[2]]
  simX_list <- sim_object$omics
  
  # PLOT 1: Data Simulation Heatmap
  plot_simData(sim_object = sim_object, type = "heatmap")
  
  # --- B. Run Method ---
  res_rad2 <- mfa_weighted_fabia(simX_list, L = 2, gamma = 0.5, p = 3, alpha = 0.7)
  
  # PLOT 2: Basic FABIA plots
  # (Standard plot method for fabia object)
  plot(res_rad2$fabia$FABIA, main = "FABIA Results")
  
  # PLOT 3: Factor Scores & Loadings (Base plots)
  factor_scores <- res_rad2$fabia$factor_scores
  factor_loadings <- res_rad2$fabia$loadings_per_block
  plot(factor_scores, main = "Factor Scores (Base Plot)")
  plot(factor_loadings$omic1, main = "Loadings Omic1 (Base Plot)")
  plot(factor_loadings$omic2, main = "Loadings Omic2 (Base Plot)")
  
  # --- C. Factor Recovery Evaluation ---
  message("Evaluating Factors...")
  # PLOT 4 & 5: Factor Correlation Heatmap & Scatter Plots (inside function)
  eval_factors <- evaluate_factor_recovery(sim_object, res_rad2, n_factors = 3, plot = TRUE, return_all = TRUE)
  print(eval_factors$matched_summary)
  
  # --- D. Loading Recovery Evaluation (Global) ---
  message("Evaluating Loadings...")
  eval_loadings <- evaluate_loading_recovery(sim_object, res_rad2$fabia, verbose = TRUE)
  
  # PLOT 6: Global Loading Correlation Heatmap
  pheatmap(eval_loadings$cor_mat, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           display_numbers = TRUE, 
           main = "Global Loading Correlation")
  
  # --- E. Per-Omic Loading Recovery (Specific Request) ---
  message("Plotting Per-Omic Correlations...")
  # PLOT 7 & 8: Per-Omic Correlation Heatmaps (e.g., Omic1 Correlation) and Scatters
  per_block_res <- per_omic_loading_recovery(
    eval_res = eval_loadings, 
    sim_object = sim_object,
    plot_heatmaps = TRUE, 
    plot_scatter = TRUE, 
    scatter_n = 5
  )
}
