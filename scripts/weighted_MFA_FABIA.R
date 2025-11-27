
library(SUMO) 
library(FactoMineR)
library(fabia)


# ---------------- SIMULATED DATA -----------------
sim1 <- simulateMultiOmics(
  vector_features = c(4000, 2500),
  n_samples = 100,
  n_factors = 3,
  snr = 2.0,
  signal.samples = c(3, 1),
  signal.features = list(
    c(4.5, 0.5),
    c(4.5, 0.5)
  ),
  factor_structure = "mixed",
  num.factor = "multiple",
  seed = 123
)
#plot_simData(sim_object = sim1, type = "heatmap")

# ---------------- SIMULATED DATA -----------------
sim2 <- simulateMultiOmics(
  vector_features = c(4000, 3500),
  n_samples = 100,
  n_factors = 3,
  snr = 1.0,
  signal.samples = c(3, 1),
  signal.features = list(
    c(2.5, 0.5),
    c(3, 2.5)
  ),
  factor_structure = "mixed",
  num.factor = "multiple",
  seed = 123
)
#plot_simData(sim_object = sim2, type = "heatmap")

# Create dataset from different simulation
omic1 = sim1$omics$omic1
omic2 = sim2$omics$omic2

sim_object = sim1
sim_object$omics$omic2 = sim2$omics$omic2
sim_object[["list_betas"]][[2]] <- sim2[["list_betas"]][[2]]

plot_simData(sim_object = sim_object, type = "heatmap")

# Final data
simX_list = sim_object$omics


# ------------------------------------------------------------
# 2) Option — MFA-derived feature weighting → FABIA
#    L = 0 → no MFA weighting (all weights = 1)
#    output structure parallel to above
# ------------------------------------------------------------
mfa_weighted_fabia <- function(Xlist,
                               L = 2,          # use first L MFA comps for weights; 0 = no weighting
                               gamma = 0.5,    # tempering
                               p = 3, alpha = 0.2,
                               block_equalize = TRUE,
                               scale_vars = TRUE,
                               seed = 1) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  stopifnot(is.list(Xlist), length(Xlist) >= 2)
  set.seed(seed)
  
  # helper to clean NA/Inf per block
  clean_block <- function(M) {
    M <- as.matrix(M)
    # drop columns that are all non-finite
    all_bad <- apply(M, 2, function(z) all(!is.finite(z)))
    if (any(all_bad)) {
      M <- M[, !all_bad, drop = FALSE]
    }
    # impute remaining NA/Inf by column mean
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
  
  # 1) preprocess ----------------------------------------------------
  Zlist <- lapply(Xlist, function(X) {
    X <- as.matrix(X)
    if (scale_vars) X <- scale(X) else X <- scale(X, scale = FALSE)
    X
  })
  
  if (block_equalize) {
    Zlist <- Map(function(Z, k) Z / sqrt(k),
                 Zlist,
                 lapply(Zlist, ncol))
  }
  
  # CLEAN HERE (this was missing!)
  Zlist <- lapply(Zlist, clean_block)
  
  names(Zlist) <- names(Zlist) %||% paste0("omic", seq_along(Zlist))
  for (g in seq_along(Zlist)) {
    if (is.null(colnames(Zlist[[g]]))) {
      colnames(Zlist[[g]]) <- paste0("V", seq_len(ncol(Zlist[[g]])))
    }
    colnames(Zlist[[g]]) <- paste0(names(Zlist)[g], "_", colnames(Zlist[[g]]))
  }
  
  Z   <- do.call(cbind, Zlist)             # n × p
  grp <- vapply(Zlist, ncol, integer(1))   # group sizes
  
  # 2) MFA -----------------------------------------------------------
  df_all <- as.data.frame(Z)
  mfa <- FactoMineR::MFA(
    df_all,
    group      = grp,
    type       = rep("s", length(grp)),
    ncp        = max(L, 5),
    name.group = names(Zlist),
    graph      = FALSE
  )
  
  # 3) build feature weights -----------------------------------------
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
  
  # 4) apply weights -------------------------------------------------
  Z_w <- sweep(Z, 2, w, `*`)   # n × p
  
  # split weighted matrix back into blocks
  cuts <- cumsum(grp)
  weighted_blocks <- lapply(seq_along(grp), function(i) {
    if (i == 1) {
      Z_w[, 1:cuts[1], drop = FALSE]
    } else {
      Z_w[, (cuts[i-1] + 1):cuts[i], drop = FALSE]
    }
  })
  names(weighted_blocks) <- names(Zlist)
  
  # 5) FABIA on weighted data ----------------------------------------
  fab_fit <- fabia::fabia(as.matrix(t(Z_w)), p = p, alpha = alpha)
  
  factor_scores <- t(fab_fit@Z)   # n × K
  loadings_all  <- fab_fit@L      # p × K
  
  idx <- split(seq_len(nrow(loadings_all)),
               rep(seq_along(grp), times = grp))
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
    weighted       = Z_w,             # weighted data
    weighted_blocks = weighted_blocks, # weighted blocks (per omic)
    fabia           = fab_res,
    weights         = w
  )
}


#### Test on SUMO simulated data -----------------------------------------------

## 2 MFA components used to compute feature weights
res_rad2 <- mfa_weighted_fabia(simX_list, L = 2, gamma = 0.5, p = 3, alpha = 0.7)
res_rad2$fabia # get plots


## plots for factor scores and loadings

factor_scores <- res_rad2[["fabia"]][["factor_scores"]]
plot(factor_scores)

factor_loadings <- res_rad2[["fabia"]][["loadings_per_block"]]
plot(factor_loadings$omic1)
plot(factor_loadings$omic2)

#### -- PERFORMANCE EVALUATION -------------------------------------------------------------------------

###########################################################################
# factor recovery for hybrid MFA→FABIA pipeline. It:
# Locates simulated true factors in sim_object,
# Computes Pearson correlation matrix between true and estimated factor scores,
# Finds the best one-to-one matching using the Hungarian algorithm 
# (maximize absolute correlation),
###########################################################################


# Required packages (install if missing)
pkgs <- c("clue", "pheatmap", "ggplot2")
to_install <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install)
library(clue); library(pheatmap); library(ggplot2)

# Use sim_object$list_alphas as true factors
get_true_factors_from_sim <- function(sim_object) {
  if (!is.null(sim_object$list_alphas) && is.list(sim_object$list_alphas)) {
    # convert list of vectors (alpha1, alpha2, ...) to matrix (n_samples x K)
    mats <- sim_object$list_alphas
    # ensure they are column vectors of same length
    mat <- do.call(cbind, lapply(mats, as.numeric))
    colnames(mat) <- names(mats)
    rownames(mat) <- if (!is.null(rownames(sim_object$omics[[1]]))) rownames(sim_object$omics[[1]]) else seq_len(nrow(mat))
    return(mat)
  }
  stop("Could not find true factors in sim_object$list_alphas.")
}

# Core evaluation function (uses Hungarian assignment on abs(correlation))
# Replace your previous evaluate_factor_recovery() with this corrected version
evaluate_factor_recovery <- function(sim_object, res, n_factors = NULL, plot = TRUE, return_all = FALSE) {
  # get estimated factors
  if (is.list(res) && !is.null(res$fabia) && !is.null(res$fabia$factor_scores)) {
    est <- as.matrix(res$fabia$factor_scores)
  } else if (!is.null(res$factor_scores)) {
    est <- as.matrix(res$factor_scores)
  } else {
    stop("Cannot find estimated factor scores in 'res'. Expected res$fabia$factor_scores or res$factor_scores.")
  }
  n <- nrow(est)
  
  # true factors from sim_object$list_alphas
  true <- get_true_factors_from_sim(sim_object)
  true <- as.matrix(true)
  if (nrow(true) != n) stop("Sample count mismatch between estimated and true factors.")
  
  K_est <- ncol(est); K_true <- ncol(true)
  if (is.null(n_factors)) n_factors <- min(K_est, K_true)
  
  estK <- est[, seq_len(min(ncol(est), n_factors)), drop = FALSE]
  trueK <- true[, seq_len(min(ncol(true), n_factors)), drop = FALSE]
  K <- ncol(trueK)
  
  # correlation matrix (true x estimated)
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  
  # handle NA correlations by setting them to 0 (no association) before building cost
  cor_mat_na <- cor_mat
  cor_mat_na[is.na(cor_mat_na)] <- 0
  
  # Build nonnegative cost matrix for solve_LSAP: use (1 - |correlation|)
  abs_cor <- abs(cor_mat_na)
  cost <- 1 - abs_cor
  # Ensure non-negative numeric matrix
  cost <- as.matrix(cost)
  if (any(cost < 0)) cost <- cost - min(cost)
  
  # Solve assignment (Hungarian) minimizing cost -> maximizes abs(correlation)
  assignment <- solve_LSAP(cost)   # now valid: non-negative entries
  assigned_cols <- as.integer(assignment)
  
  # compute matched correlations and signs (use original cor_mat, but replace NA with 0 for numeric)
  matched_corrs <- numeric(K); matched_signs <- numeric(K)
  for (i in seq_len(K)) {
    j <- assigned_cols[i]
    cval <- cor_mat[i, j]
    if (is.na(cval)) cval <- 0
    matched_corrs[i] <- cval
    matched_signs[i] <- sign(cval)
    if (matched_signs[i] == 0) matched_signs[i] <- 1  # if correlation == 0 choose sign +1
  }
  
  # sign-adjust matched estimated factors
  est_matched_signed <- estK[, assigned_cols, drop = FALSE] %*% diag(matched_signs)
  colnames(est_matched_signed) <- paste0("Est_matched_", seq_len(K))
  
  # metrics
  mean_abs_cor <- mean(abs(matched_corrs), na.rm = TRUE)
  prop_above_0.7 <- mean(abs(matched_corrs) >= 0.7, na.rm = TRUE)
  prop_above_0.5 <- mean(abs(matched_corrs) >= 0.5, na.rm = TRUE)
  
  r2 <- numeric(K); mse <- numeric(K)
  for (i in seq_len(K)) {
    fit <- lm(trueK[, i] ~ est_matched_signed[, i])
    ssr <- sum((predict(fit) - mean(trueK[, i]))^2)
    sst <- sum((trueK[, i] - mean(trueK[, i]))^2)
    r2[i] <- if (sst == 0) NA else ssr / sst
    mse[i] <- mean((trueK[, i] - est_matched_signed[, i])^2)
  }
  
  if (plot) {
    pheatmap(cor_mat,
             main = "Correlation (true factors × estimated factors)",
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE)
    df_list <- lapply(seq_len(K), function(i) {
      data.frame(sample = seq_len(n),
                 true = trueK[, i],
                 est = est_matched_signed[, i],
                 factor = paste0("True_", i))
    })
    df <- do.call(rbind, df_list)
    g <- ggplot(df, aes(x = true, y = est)) +
      geom_point(alpha = 0.6) +
      facet_wrap(~factor, scales = "free") +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = "True factor score", y = "Estimated (signed) factor score",
           title = "True vs Estimated (matched & sign-adjusted)")
    print(g)
  }
  
  summary_tbl <- data.frame(
    true_factor = colnames(trueK),
    matched_estimated_index = assigned_cols,
    pearson_corr = matched_corrs,
    abs_corr = abs(matched_corrs),
    sign = matched_signs,
    R2 = r2,
    MSE = mse,
    stringsAsFactors = FALSE
  )
  
  overall <- list(
    correlation_matrix = cor_mat,
    assignment = assigned_cols,
    matched_summary = summary_tbl,
    mean_abs_cor = mean_abs_cor,
    prop_abs_cor_ge_0.7 = prop_above_0.7,
    prop_abs_cor_ge_0.5 = prop_above_0.5
  )
  if (return_all) {
    overall$est_signed_matched <- est_matched_signed
    overall$true_matched <- trueK
  }
  return(overall)
}

# -------------------------
# Example usage (run after res_rad2 exists)
# -------------------------
res_rad3 <- mfa_weighted_fabia(simX_list, L = 3, gamma = 0.5, p = 3, alpha = 0.2)
evaluation <- evaluate_factor_recovery(sim_object = sim_object, res = res_rad3, n_factors = 3, plot = TRUE, return_all = TRUE)
print(evaluation$matched_summary)
cat(sprintf("Mean |corr| = %.3f  |  prop |corr|>=0.7 = %.2f\n", evaluation$mean_abs_cor, evaluation$prop_abs_cor_ge_0.7))



####-- Diagonize why alpha3 wasnt recovered
# packages
pkgs <- c("pheatmap", "ggplot2", "matrixStats")
to_install <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(to_install)) install.packages(to_install)
library(pheatmap); library(ggplot2); library(matrixStats)

# 1) Basic properties of true factors (list_alphas)
true_factors <- do.call(cbind, lapply(sim_object$list_alphas, as.numeric))
colnames(true_factors) <- names(sim_object$list_alphas)
sapply(as.data.frame(true_factors), function(x) c(mean = mean(x), sd = sd(x), var = var(x)))

# Show number of signal samples for each factor (from signal_annotation)
sig_samples <- sim_object$signal_annotation$samples
sapply(sig_samples, length)   # e.g. factor1=35, factor2=18, factor3=13

# Plot the true factor time/score distributions
df_tf <- data.frame(sample = 1:nrow(true_factors), true_factors)
library(reshape2)
dfm <- melt(df_tf, id.vars = "sample")
ggplot(dfm, aes(x = value)) + geom_histogram(bins = 30) + facet_wrap(~variable, scales = "free") +
  ggtitle("Distribution of true factor scores (α1..α3)")

# 2) Check SNR / per-omic rough signal strength
# Compute for each omic the variance explained by each true factor via simple regression:
omics <- sim_object$omics
for (g in names(omics)) {
  cat("\n=== Omic:", g, "===\n")
  X <- as.matrix(omics[[g]])   # samples x features
  # project each feature on true factors and compute average |beta| or R2
  betas <- apply(X, 2, function(f) {
    # regress f ~ true_factors (use all K)
    fit <- lm(f ~ true_factors)
    # return adjusted R^2 (or coefficients)
    summary(fit)$adj.r.squared
  })
  cat("Mean adj. R2 across features:", mean(betas, na.rm=TRUE), "\n")
  cat("Median adj. R2 across features:", median(betas, na.rm=TRUE), "\n")
  # Show fraction of features with adj.R2 > small thresholds
  cat("Frac features adj.R2 > 0.01:", mean(betas > 0.01, na.rm=TRUE), "\n")
  cat("Frac features adj.R2 > 0.05:", mean(betas > 0.05, na.rm=TRUE), "\n")
}

# 3) Inspect true loadings (list_betas). Identify factor3's loadings and magnitude
str(sim_object$list_betas)
# list_betas is list of lists per-omic -> find beta vectors for factor3
for (g in seq_along(sim_object$list_betas)) {
  cat("\nOmic index", g, ":\n")
  betalist <- sim_object$list_betas[[g]]
  print(names(betalist))
  # if there is a beta3 print summaries
  if ("beta3" %in% names(betalist)) {
    b3 <- as.numeric(betalist[["beta3"]])
    cat("beta3: mean abs =", mean(abs(b3)), "sd abs =", sd(abs(b3)), "nonzero count =", sum(b3 != 0), "\n")
    # optionally show top absolute loadings
    print(head(sort(abs(b3), decreasing = TRUE), 10))
  } else {
    cat("No beta3 in this omic.\n")
  }
}

# 4) Check where factor3 signals live: samples and features
cat("Signal sample indices for factor3:\n")
print(sim_object$signal_annotation$samples$factor3)   # small set (13)
cat("Number of factor3 features per omic (from signal_annotation):\n")
lapply(sim_object$signal_annotation$features, function(x) {
  if (!is.null(x$factor3)) length(x$factor3) else 0
})

# 5) Visualize factor3-related features from omic2 (heatmap)
# if factor3 features exist in omic2:
if (!is.null(sim_object$signal_annotation$features$omic2$factor3)) {
  idx_feat <- sim_object$signal_annotation$features$omic2$factor3
  X2 <- as.matrix(sim_object$omics$omic2)
  # subset: factor3 features (rows = samples, cols = features)
  Xsub <- X2[, idx_feat, drop = FALSE]
  # scale per feature for visualization
  Xsub_sc <- scale(Xsub)
  # reorder samples by true factor3 score
  ord <- order(true_factors[, "alpha3"])
  pheatmap(t(Xsub_sc[ord, , drop = FALSE]), show_rownames = FALSE, show_colnames = FALSE,
           main = paste("Heatmap of omic2 features for factor3 (nfeat=", ncol(Xsub), ")"))
} else {
  cat("No factor3 features recorded in omic2 signal_annotation.\n")
}

# 6) Check estimated loadings for the matched estimated factor (from evaluation result)
# Assuming you have 'evaluation' object from previous run, with assignment
if (exists("evaluation")) {
  print(evaluation$matched_summary)
  matched_idx <- evaluation$matched_summary$matched_estimated_index
  # matched_idx[i] is the column index in est factors; find which matched estimated factor corresponded to true alpha3
  cat("Matched estimated index for true alpha3:", matched_idx[3], "\n")
  # look at fabia loadings_per_block for that estimated factor
  if (!is.null(res_rad2$fabia$loadings_per_block)) {
    lpb <- res_rad2$fabia$loadings_per_block
    for (g in names(lpb)) {
      Lg <- lpb[[g]]
      if (ncol(Lg) >= matched_idx[3]) {
        ld <- Lg[, matched_idx[3]]
        cat("Omic", g, ": mean abs loading for matched factor:", mean(abs(ld)), "max abs:", max(abs(ld)), "\n")
      } else cat("Omic", g, ": matched factor index out of range\n")
    }
  }
}

#######################################################################################
# Evaluating loading recovery
#####################################################################################
# ---------------------------
# Complete evaluation script
# ---------------------------
# Required packages
if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(clue)
library(ggplot2)

# ---------------------------
# 1) evaluate_loading_recovery()
# ---------------------------
evaluate_loading_recovery <- function(sim_object, fab_res, verbose = TRUE) {
  # Determine true K from list_alphas
  K <- length(sim_object$list_alphas)
  if (K < 1) stop("sim_object$list_alphas must contain at least one factor.")
  # Block names and sizes
  block_names <- names(sim_object$omics)
  block_sizes <- vapply(sim_object$omics, ncol, integer(1))
  p_total <- sum(block_sizes)
  # Build true_all (p_total x K) placing betas into correct columns (zeros where missing)
  true_all <- matrix(0, nrow = p_total, ncol = K)
  colnames(true_all) <- paste0("factor", seq_len(K))
  row_start <- 1
  for (b in seq_along(block_names)) {
    bname <- block_names[b]
    ncols <- block_sizes[b]
    row_end <- row_start + ncols - 1
    betas_block <- sim_object$list_betas[[b]]
    if (!is.null(betas_block) && is.list(betas_block)) {
      beta_names <- names(betas_block)
      if (is.null(beta_names)) {
        # unnamed list: map sequentially to factors
        for (j in seq_along(betas_block)) {
          target_col <- j
          if (target_col <= K) {
            vec <- as.numeric(betas_block[[j]])
            if (length(vec) == ncols) true_all[row_start:row_end, target_col] <- vec
            else warning(sprintf("beta length (%d) != block size (%d) for block %s, col %d — skipping.",
                                 length(vec), ncols, bname, target_col))
          }
        }
      } else {
        # named list, expect names like "beta1","beta2",...
        for (nm in beta_names) {
          m <- regmatches(nm, regexec("beta(\\d+)$", nm))
          if (length(m) && length(m[[1]]) >= 2) {
            target_col <- as.integer(m[[1]][2])
          } else {
            target_col <- suppressWarnings(as.integer(nm))
            if (is.na(target_col)) next
          }
          if (target_col >= 1 && target_col <= K) {
            vec <- as.numeric(betas_block[[nm]])
            if (length(vec) == ncols) true_all[row_start:row_end, target_col] <- vec
            else warning(sprintf("beta '%s' length (%d) != block size (%d) for block %s — skipping.",
                                 nm, length(vec), ncols, bname))
          } else {
            warning(sprintf("beta '%s' refers to factor %s out of range 1..%d; ignoring.", nm, target_col, K))
          }
        }
      }
    } else if (!is.null(betas_block)) {
      # single numeric vector -> assume factor1
      vec <- as.numeric(betas_block)
      if (length(vec) == ncols) true_all[row_start:row_end, 1] <- vec
      else warning(sprintf("Unexpected beta format for block %s; skipping.", bname))
    }
    row_start <- row_end + 1
  }
  
  # Extract estimated loadings from FABIA result (p_total x K_est)
  est_all <- fab_res$FABIA@L
  if (!is.matrix(est_all)) est_all <- as.matrix(est_all)
  if (nrow(est_all) != nrow(true_all)) stop(sprintf(
    "Estimated loadings rows (%d) != true loadings rows (%d).",
    nrow(est_all), nrow(true_all)
  ))
  K_est <- ncol(est_all)
  K_use <- min(K, K_est)
  
  # Correlation matrix (true cols vs estimated cols)
  cor_mat <- stats::cor(true_all[, seq_len(K_use), drop = FALSE],
                        est_all[, seq_len(K_use), drop = FALSE],
                        use = "pairwise.complete.obs")
  abs_cor <- abs(cor_mat)
  abs_cor[!is.finite(abs_cor)] <- 0
  
  # Build non-negative cost for solve_LSAP: cost = max_val - abs_cor
  max_val <- max(abs_cor, na.rm = TRUE)
  if (!is.finite(max_val)) max_val <- 0
  cost_block <- matrix(max_val, nrow = nrow(abs_cor), ncol = ncol(abs_cor)) - abs_cor
  nr <- nrow(cost_block); nc <- ncol(cost_block); n <- max(nr, nc)
  cost <- matrix(max_val, n, n)
  cost[1:nr, 1:nc] <- cost_block
  
  # Solve assignment
  assignment <- clue::solve_LSAP(cost)
  assigned_pairs <- cbind(row = seq_len(nr), col = as.integer(assignment[seq_len(nr)]))
  assigned_pairs <- assigned_pairs[assigned_pairs[,2] <= nc, , drop = FALSE]
  
  # Flip signs of estimated components to maximize positive correlation, compute RMSE
  matched <- data.frame(true_comp = integer(0), est_comp = integer(0),
                        cor = numeric(0), abs_cor = numeric(0), rmse = numeric(0),
                        stringsAsFactors = FALSE)
  est_all_flipped <- est_all
  for (i in seq_len(nrow(assigned_pairs))) {
    tr <- assigned_pairs[i, "row"]; es <- assigned_pairs[i, "col"]
    cc <- suppressWarnings(cor(true_all[, tr], est_all[, es], use = "pairwise.complete.obs"))
    if (!is.finite(cc)) cc <- 0
    if (cc < 0) {
      est_all_flipped[, es] <- -est_all_flipped[, es]
      cc <- -cc
    }
    dif <- true_all[, tr] - est_all_flipped[, es]
    rmse <- sqrt(mean(dif^2, na.rm = TRUE))
    matched <- rbind(matched, data.frame(true_comp = tr, est_comp = es, cor = cc, abs_cor = abs(cc), rmse = rmse))
  }
  
  # Summary metrics
  mean_abs_cor <- mean(matched$abs_cor, na.rm = TRUE)
  median_abs_cor <- median(matched$abs_cor, na.rm = TRUE)
  mean_rmse <- mean(matched$rmse, na.rm = TRUE)
  
  # Per-block matched RMSE & cor
  cuts <- cumsum(block_sizes)
  idx_list <- lapply(seq_along(block_sizes), function(i) {
    if (i == 1) seq_len(cuts[1]) else (cuts[i-1] + 1):cuts[i]
  })
  names(idx_list) <- block_names
  
  per_block <- lapply(names(idx_list), function(block_name) {
    idx <- idx_list[[block_name]]
    true_block <- true_all[idx, , drop = FALSE]
    est_block <- est_all_flipped[idx, , drop = FALSE]
    matched_block <- do.call(rbind, lapply(seq_len(nrow(matched)), function(i) {
      data.frame(block = block_name,
                 true_comp = matched$true_comp[i],
                 est_comp = matched$est_comp[i],
                 cor = matched$cor[i],
                 rmse_block = sqrt(mean((true_block[, matched$true_comp[i]] - est_block[, matched$est_comp[i]])^2, na.rm = TRUE)),
                 stringsAsFactors = FALSE)
    }))
    list(block = block_name,
         cormat = tryCatch(stats::cor(true_block[, seq_len(K_use), drop = FALSE],
                                      est_block[, seq_len(K_use), drop = FALSE],
                                      use = "pairwise.complete.obs"),
                           error = function(e) matrix(NA, K_use, K_use)),
         matched = matched_block)
  })
  names(per_block) <- names(idx_list)
  
  result <- list(
    true_all = true_all,
    est_all = est_all,
    est_all_flipped = est_all_flipped,
    cor_mat = cor_mat,
    abs_cor = abs_cor,
    assignment = assigned_pairs,
    matched = matched,
    summary = list(mean_abs_cor = mean_abs_cor, median_abs_cor = median_abs_cor, mean_rmse = mean_rmse),
    per_block = per_block,
    block_sizes = block_sizes
  )
  
  if (verbose) {
    cat("=== Loading recovery summary ===\n")
    cat(sprintf("Components compared (K): %d (true) vs %d (estimated) -> using %d\n", K, K_est, K_use))
    cat(sprintf("Mean absolute matched correlation: %.4f\n", mean_abs_cor))
    cat(sprintf("Median absolute matched correlation: %.4f\n", median_abs_cor))
    cat(sprintf("Mean matched RMSE: %.4f\n\n", mean_rmse))
    cat("Matched pairs (true -> est) with correlations and RMSE:\n"); print(matched)
    cat("\nPer-block matched RMSE (summary):\n")
    block_summ <- do.call(rbind, lapply(per_block, function(pb) {
      data.frame(block = pb$block,
                 mean_rmse = mean(pb$matched$rmse_block, na.rm = TRUE),
                 median_rmse = median(pb$matched$rmse_block, na.rm = TRUE),
                 mean_cor = mean(pb$matched$cor, na.rm = TRUE),
                 median_cor = median(pb$matched$cor, na.rm = TRUE))
    }))
    print(block_summ)
  }
  
  return(result)
}

# ---------------------------
# 2) plotting helpers
# ---------------------------
plot_correlation_heatmap <- function(cormat, main = "Correlation (true vs estimated)") {
  if (is.null(rownames(cormat))) rownames(cormat) <- paste0("T", seq_len(nrow(cormat)))
  if (is.null(colnames(cormat))) colnames(cormat) <- paste0("E", seq_len(ncol(cormat)))
  # base image plot
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mar = c(5,5,4,2))
  image(1:ncol(cormat), 1:nrow(cormat), t(cormat[nrow(cormat):1, , drop = FALSE]),
        axes = FALSE, xlab = "Estimated components", ylab = "True components", main = main)
  axis(1, at = 1:ncol(cormat), labels = colnames(cormat), las = 2)
  axis(2, at = 1:nrow(cormat), labels = rev(rownames(cormat)), las = 2)
  for (i in 1:nrow(cormat)) for (j in 1:ncol(cormat)) text(j, nrow(cormat)-i+1, sprintf("%.2f", cormat[i,j]), cex = 0.8)
}

plot_matched_scatter <- function(true_all, est_all_flipped, matched_df, n_show = nrow(matched_df)) {
  for (i in seq_len(min(n_show, nrow(matched_df)))) {
    tcol <- matched_df$true_comp[i]
    ecol <- matched_df$est_comp[i]
    dat <- data.frame(true = true_all[, tcol], est = est_all_flipped[, ecol])
    p <- ggplot(dat, aes(x = true, y = est)) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "lm", se = FALSE) +
      labs(title = sprintf("True comp %d vs Est comp %d (r=%.3f)", tcol, ecol, matched_df$cor[i]),
           x = "True loading", y = "Estimated loading")
    print(p)
  }
}

# ---------------------------
# 3) support recovery helper
#     - two modes for reference truth:
#        a) use non-zero entries in true_all (default)
#        b) use sim_object$signal_annotation$features if present and requested
# ---------------------------
support_recovery <- function(eval_res, threshold = NULL, top_n = NULL,
                             use_signal_annotation = FALSE, verbose = TRUE) {
  # eval_res: result object from evaluate_loading_recovery()
  true_all <- eval_res$true_all
  est_all <- eval_res$est_all_flipped  # already sign-flipped to match true
  block_sizes <- eval_res$block_sizes
  block_names <- names(block_sizes)
  K <- ncol(true_all)
  p_total <- nrow(true_all)
  
  # Build true support sets (list of integer vectors of feature indices) per factor
  if (use_signal_annotation && !is.null(sim_object$signal_annotation$features)) {
    # try to extract from sim_object$signal_annotation$features; handle nested structure
    true_support <- vector("list", K)
    names(true_support) <- paste0("factor", seq_len(K))
    # signal_annotation$features may have per-block names (omic1, omic2) and factor lists inside
    pos <- 1
    for (b in seq_along(sim_object$signal_annotation$features)) {
      bname <- names(sim_object$signal_annotation$features)[b]
      block_feats <- sim_object$signal_annotation$features[[bname]]
      # block_feats may be list of factorX -> integer vector of indices (global or local)
      if (is.list(block_feats) && length(block_feats) > 0) {
        for (fnm in names(block_feats)) {
          # determine factor index from name like "factor2"
          m <- regmatches(fnm, regexec("factor(\\d+)$", fnm))
          if (length(m) && length(m[[1]]) >= 2) {
            fidx <- as.integer(m[[1]][2])
            if (fidx >= 1 && fidx <= K) {
              # add these indices to true_support[fidx]
              true_support[[fidx]] <- unique(c(true_support[[fidx]], as.integer(block_feats[[fnm]])))
            }
          }
        }
      }
    }
    # If annotation uses local (per-block) indices, you may need to convert to global indices.
    # Here we assume stored indices are global (consistent with sim_object concatenated_datasets).
  } else {
    # Default: true support = non-zero entries in true_all (absolute > 0)
    true_support <- lapply(seq_len(K), function(k) which(abs(true_all[, k]) > 0))
    names(true_support) <- paste0("factor", seq_len(K))
  }
  
  # Build estimated support sets: thresholding absolute loading values
  est_support <- vector("list", ncol(est_all))
  for (k in seq_len(ncol(est_all))) {
    vec <- abs(est_all[, k])
    if (!is.null(threshold)) {
      est_idx <- which(vec >= threshold)
    } else if (!is.null(top_n)) {
      # choose top_n features by absolute loading
      ord <- order(vec, decreasing = TRUE)
      est_idx <- head(ord, top_n)
    } else {
      # default threshold: choose features with abs loading > 99th percentile (sparse)
      tval <- quantile(vec, 0.99, na.rm = TRUE)
      est_idx <- which(vec >= tval)
    }
    est_support[[k]] <- est_idx
  }
  names(est_support) <- paste0("est", seq_len(ncol(est_all)))
  
  # Map estimated components to true components using matched pairs from eval_res
  mapping <- eval_res$matched
  # For each matched pair compute precision/recall/F1
  pr_table <- data.frame(true_comp = integer(0), est_comp = integer(0),
                         TP = integer(0), FP = integer(0), FN = integer(0),
                         precision = numeric(0), recall = numeric(0), F1 = numeric(0),
                         stringsAsFactors = FALSE)
  for (i in seq_len(nrow(mapping))) {
    tcomp <- mapping$true_comp[i]
    ecomp <- mapping$est_comp[i]
    true_set <- true_support[[tcomp]]
    est_set <- est_support[[ecomp]]
    TP <- length(intersect(true_set, est_set))
    FP <- length(setdiff(est_set, true_set))
    FN <- length(setdiff(true_set, est_set))
    prec <- if ((TP + FP) == 0) NA else TP / (TP + FP)
    rec <- if ((TP + FN) == 0) NA else TP / (TP + FN)
    F1 <- if (is.na(prec) || is.na(rec) || (prec + rec) == 0) NA else 2 * prec * rec / (prec + rec)
    pr_table <- rbind(pr_table, data.frame(true_comp = tcomp, est_comp = ecomp,
                                           TP = TP, FP = FP, FN = FN,
                                           precision = prec, recall = rec, F1 = F1,
                                           stringsAsFactors = FALSE))
  }
  
  if (verbose) {
    cat("=== Support recovery (precision/recall) ===\n")
    print(pr_table)
    cat("\nAverage F1 (across matched comps):", mean(pr_table$F1, na.rm = TRUE), "\n")
  }
  
  return(list(true_support = true_support, est_support = est_support, pr_table = pr_table))
}

# ---------------------------
# 4) Example run (use your objects)
# ---------------------------
# Make sure 'sim_object' and 'res_rad2' exist in your environment and res_rad2$fabia is present.
# Run evaluation:
eval_res <- evaluate_loading_recovery(sim_object = sim_object, fab_res = res_rad2$fabia, verbose = TRUE)

# Plot correlation heatmap
plot_correlation_heatmap(eval_res$cor_mat)

# Plot matched scatterplots (true vs estimated loadings)
plot_matched_scatter(eval_res$true_all, eval_res$est_all_flipped, eval_res$matched)

# Compute support recovery using default thresholding (99th percentile per estimated component)
supp_res <- support_recovery(eval_res, threshold = NULL, top_n = NULL, use_signal_annotation = FALSE, verbose = TRUE)

# If you'd prefer to use signal_annotation$features as the ground truth (if indices are global), run:
if (!is.null(sim_object$signal_annotation$features)) {
  supp_res2 <- support_recovery(eval_res, threshold = NULL, top_n = NULL, use_signal_annotation = TRUE, verbose = TRUE)
}

# Save evaluation objects if desired
# save(eval_res, supp_res, file = "fabia_evaluation_results.RData")


####
# ------------------------------
# Per-omic loading-recovery
# ------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

per_omic_loading_recovery <- function(eval_res, sim_object, plot_heatmaps = TRUE, plot_scatter = TRUE, scatter_n = NULL) {
  # eval_res: result of evaluate_loading_recovery()
  # sim_object: simulation object used to build true loadings
  # scatter_n: number of matched pairs to show per block (NULL -> show all)
  true_all <- eval_res$true_all
  est_all_flipped <- eval_res$est_all_flipped
  matched <- eval_res$matched       # data.frame(true_comp, est_comp, cor, rmse)
  block_sizes <- eval_res$block_sizes
  block_names <- names(block_sizes)
  cuts <- cumsum(block_sizes)
  idx_list <- lapply(seq_along(block_sizes), function(i) {
    if (i == 1) seq_len(cuts[1]) else (cuts[i-1] + 1):cuts[i]
  })
  names(idx_list) <- block_names
  
  per_block_summary <- list()
  
  for (bname in block_names) {
    idx <- idx_list[[bname]]
    true_block <- true_all[idx, , drop = FALSE]           # rows = features in this omic
    est_block <- est_all_flipped[idx, , drop = FALSE]
    # correlation matrix for block (true factors vs estimated factors)
    # keep same number of columns as available
    Kt <- ncol(true_block); Ke <- ncol(est_block)
    # if a column is all zeros in both, cor returns NA — replace with 0 after
    cormat_block <- tryCatch({
      stats::cor(true_block, est_block, use = "pairwise.complete.obs")
    }, error = function(e) {
      matrix(NA, nrow = Kt, ncol = Ke)
    })
    cormat_block[!is.finite(cormat_block)] <- 0
    
    # For summary: use global matched pairs (which were determined on full p)
    # Find matched rows that map to valid estimated components (in range)
    matched_in_block <- do.call(rbind, lapply(seq_len(nrow(matched)), function(i) {
      tcomp <- matched$true_comp[i]; ecomp <- matched$est_comp[i]
      # compute block-specific correlation and block-specific RMSE
      cor_block_val <- 0
      rmse_block_val <- NA
      if (tcomp <= ncol(true_block) && ecomp <= ncol(est_block)) {
        cor_block_val <- suppressWarnings(cor(true_block[, tcomp], est_block[, ecomp], use = "pairwise.complete.obs"))
        if (!is.finite(cor_block_val)) cor_block_val <- 0
        dif <- true_block[, tcomp] - est_block[, ecomp]
        rmse_block_val <- sqrt(mean(dif^2, na.rm = TRUE))
      } else {
        # If the matched pair refers to a component not present in block (e.g., true factor has no effect here),
        # then correlation will be NA/0 and RMSE computed across zeros/est vector.
        cor_block_val <- NA
        rmse_block_val <- NA
      }
      data.frame(true_comp = tcomp, est_comp = ecomp,
                 cor_block = cor_block_val, rmse_block = rmse_block_val,
                 stringsAsFactors = FALSE)
    }))
    
    # Summaries: mean absolute correlation (ignoring NA), mean RMSE (ignoring NA)
    mean_abs_cor_block <- mean(abs(matched_in_block$cor_block), na.rm = TRUE)
    median_abs_cor_block <- median(abs(matched_in_block$cor_block), na.rm = TRUE)
    mean_rmse_block <- mean(matched_in_block$rmse_block, na.rm = TRUE)
    
    per_block_summary[[bname]] <- list(
      cormat = cormat_block,
      matched_block = matched_in_block,
      summary = list(mean_abs_cor = mean_abs_cor_block,
                     median_abs_cor = median_abs_cor_block,
                     mean_rmse = mean_rmse_block)
    )
    
    # print summary for this block
    cat("\n------ Block:", bname, "------\n")
    cat(sprintf("Features: %d | Matched pairs reported: %d\n", length(idx), nrow(matched_in_block)))
    cat(sprintf("Mean abs matched correlation (block): %.4f\n", mean_abs_cor_block))
    cat(sprintf("Median abs matched correlation (block): %.4f\n", median_abs_cor_block))
    cat(sprintf("Mean matched RMSE (block): %.4f\n\n", mean_rmse_block))
    cat("Matched pairs (block-level cor / rmse) -- NA means component not present in this block:\n")
    print(matched_in_block)
    
    # plotting: correlation heatmap for this block
    if (plot_heatmaps) {
      cm <- cormat_block
      if (all(is.na(cm))) {
        message("Block ", bname, ": correlation matrix is all NA/zero, skipping heatmap.")
      } else {
        # tidy labels
        rownames(cm) <- paste0("T", seq_len(nrow(cm)))
        colnames(cm) <- paste0("E", seq_len(ncol(cm)))
        op <- par(no.readonly = TRUE)
        on.exit(par(op), add = TRUE)
        par(mfrow = c(1,1), mar = c(5,5,4,2))
        image(1:ncol(cm), 1:nrow(cm), t(cm[nrow(cm):1, , drop = FALSE]),
              axes = FALSE, xlab = "Estimated components", ylab = "True components",
              main = paste0("Block: ", bname, " (corr)"))
        axis(1, at = 1:ncol(cm), labels = colnames(cm), las = 2)
        axis(2, at = 1:nrow(cm), labels = rev(rownames(cm)), las = 2)
        for (i in 1:nrow(cm)) for (j in 1:ncol(cm)) text(j, nrow(cm)-i+1, sprintf("%.2f", cm[i,j]), cex = 0.7)
      }
    }
    
    # plotting: scatterplots for matched pairs restricted to block features
    if (plot_scatter) {
      to_plot <- matched_in_block
      if (!is.null(scatter_n) && is.numeric(scatter_n)) to_plot <- to_plot[seq_len(min(nrow(to_plot), scatter_n)), , drop = FALSE]
      for (ri in seq_len(nrow(to_plot))) {
        tc <- to_plot$true_comp[ri]; ec <- to_plot$est_comp[ri]
        if (is.na(tc) || is.na(ec)) next
        # ensure tc/ec exist in block columns (they are indices relative to full K, but true_block and est_block have same K cols)
        if (tc > ncol(true_block) || ec > ncol(est_block)) next
        df <- data.frame(true = true_block[, tc], est = est_block[, ec])
        p <- ggplot(df, aes(x = true, y = est)) +
          geom_point(alpha = 0.4) +
          geom_smooth(method = "lm", se = FALSE) +
          labs(title = sprintf("%s: True comp %d vs Est comp %d (block only)", bname, tc, ec),
               subtitle = sprintf("block cor=%.3f | block RMSE=%.4f", to_plot$cor_block[ri], to_plot$rmse_block[ri]),
               x = "True loading (block features)", y = "Estimated loading (block features)")
        print(p)
      }
    }
  }
  
  return(per_block_summary)
}

# ------------------------------
# Run per-omic recovery (example)
# ------------------------------
# Ensure eval_res exists (from previous evaluate_loading_recovery call)
if (!exists("eval_res")) stop("eval_res not found. Run evaluate_loading_recovery() first and assign to eval_res.")
per_block_res <- per_omic_loading_recovery(eval_res = eval_res, sim_object = sim_object,
                                           plot_heatmaps = TRUE, plot_scatter = TRUE, scatter_n = 5)

# per_block_res contains per-block cormats, matched-block tables, and summaries.
# Example: show summary for omic1
if ("omic1" %in% names(per_block_res)) {
  print(per_block_res$omic1$summary)
  print(per_block_res$omic1$matched_block)
}

