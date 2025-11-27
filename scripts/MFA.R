
# Cleaned and formatted version of the Robust MFA block (drop-in)

# -------------------- Helpers & Packages --------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y

library(FactoMineR)   # for MFA
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clue)
library(fabia)


# -------------------- SIMULATED DATA --------------------
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
omic1 <- sim1$omics$omic1
omic2 <- sim2$omics$omic2

sim_object <- sim1
sim_object$omics$omic2 <- sim2$omics$omic2
sim_object[["list_betas"]][[2]] <- sim2[["list_betas"]][[2]]

# show heatmap of the combined simulated data (optional)
plot_simData(sim_object = sim_object, type = "heatmap")

# Final data list (named)
simX_list <- sim_object$omics
names(simX_list) <- paste0("omic", seq_along(simX_list))


# -------------------- 1) True Factor Scores --------------------
# Extract true scores from sim_object$list_alphas
true_scores <- do.call(cbind, lapply(sim_object$list_alphas, as.numeric))
colnames(true_scores) <- names(sim_object$list_alphas)
rownames(true_scores) <- rownames(sim_object$omics[[1]]) %||% seq_len(nrow(true_scores))

n_factors <- ncol(true_scores)

df_scores_true <- as.data.frame(true_scores) %>%
  mutate(Sample = rownames(true_scores)) %>%
  pivot_longer(cols = -Sample, names_to = "Factor", values_to = "Score")

ggplot(df_scores_true, aes(x = Sample, y = Score, color = Factor, group = Factor)) +
  geom_point(alpha = 0.6) +
  #geom_line(alpha = 0.3) +
  facet_wrap(~Factor, scales = "free_y") +
  labs(title = "True Factor Scores (Simulated Data)", x = "Sample", y = "Score") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# -------------------- 2) True Loadings per Block --------------------
get_true_loadings_per_block <- function(sim_object) {
  true_loadings_list <- list()
  for (i in seq_along(sim_object$list_betas)) {
    block_betas <- sim_object$list_betas[[i]]
    block_mats <- lapply(block_betas, function(b) matrix(as.numeric(b), ncol = 1))
    block_combined <- do.call(cbind, block_mats)
    rownames(block_combined) <- paste0("omic", i, "_feature_", seq_len(nrow(block_combined)))
    true_loadings_list[[i]] <- block_combined
  }
  names(true_loadings_list) <- paste0("omic", seq_along(sim_object$list_betas))
  true_loadings_list
}

true_loadings <- get_true_loadings_per_block(sim_object)

# Convert to long format for plotting
loadings_long <- lapply(names(true_loadings), function(block) {
  mat <- true_loadings[[block]]
  df <- as.data.frame(mat)
  df$Feature <- rownames(mat)
  df$Block <- block
  pivot_longer(df, cols = -c(Feature, Block), names_to = "Factor", values_to = "Loading")
})

df_loadings_true <- bind_rows(loadings_long)

ggplot(df_loadings_true, aes(x = Feature, y = Loading, color = Factor)) +
  geom_point(alpha = 0.4, size = 0.5) +
  facet_wrap(~Block + Factor, scales = "free_y") +
  labs(title = "True Loadings per Block (Simulated Data)", y = "Loading") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Quick base-R plots (optional)
par(mfrow = c(2, 2))
plot(sim_object[["list_betas"]][[1]][["beta2"]], main = "beta2 omic1", ylab = "loadings", cex = 0.4)
plot(sim_object[["list_betas"]][[2]][["beta1"]], main = "beta1 omic2", ylab = "loadings", cex = 0.4)
plot(sim_object[["list_betas"]][[2]][["beta2"]], main = "beta2 omic2", ylab = "loadings", cex = 0.4)
plot(sim_object[["list_betas"]][[2]][["beta3"]], main = "beta3 omic2", ylab = "loadings", cex = 0.4)
par(mfrow = c(1, 1))


# -------------------- Evaluation Helpers --------------------
evaluate_scores <- function(true_mat, est_mat, n_factors_eval = NULL, plots = TRUE, title_prefix = "") {
  true_mat <- as.matrix(true_mat)
  est_mat  <- as.matrix(est_mat)
  
  if (nrow(true_mat) != nrow(est_mat)) stop("sample rows mismatch between true and estimated score matrices")
  
  K_true <- ncol(true_mat)
  K_est  <- ncol(est_mat)
  n_f <- n_factors_eval %||% min(K_true, K_est)
  
  trueK <- true_mat[, seq_len(n_f), drop = FALSE]
  estK  <- est_mat[, seq_len(min(K_est, n_f)), drop = FALSE]
  
  # correlation matrix (true factors as rows, est factors as cols)
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  cor_mat0 <- cor_mat
  cor_mat0[is.na(cor_mat0)] <- 0
  abs_cor <- abs(cor_mat0)
  cost <- 1 - abs_cor
  
  assign <- clue::solve_LSAP(as.matrix(cost))   # returns column index for each row
  cols <- as.integer(assign)
  
  # get matched corrs + signs
  matched_corrs <- numeric(n_f)
  matched_signs <- numeric(n_f)
  for (i in seq_len(n_f)) {
    j <- cols[i]
    cval <- cor_mat[i, j]
    if (is.na(cval)) cval <- 0
    matched_corrs[i] <- cval
    sgn <- sign(cval)
    if (sgn == 0) sgn <- 1
    matched_signs[i] <- sgn
  }
  
  # construct matched estimated (signed)
  est_matched_signed <- estK[, cols, drop = FALSE] %*% diag(matched_signs)
  colnames(est_matched_signed) <- paste0("Est_matched_", seq_len(n_f))
  
  # metrics
  mean_abs_cor <- mean(abs(matched_corrs), na.rm = TRUE)
  prop_ge_0.7 <- mean(abs(matched_corrs) >= 0.7, na.rm = TRUE)
  prop_ge_0.5 <- mean(abs(matched_corrs) >= 0.5, na.rm = TRUE)
  
  R2  <- numeric(n_f)
  MSE <- numeric(n_f)
  for (i in seq_len(n_f)) {
    fit <- lm(trueK[, i] ~ est_matched_signed[, i])
    ssr <- sum((predict(fit) - mean(trueK[, i]))^2)
    sst <- sum((trueK[, i] - mean(trueK[, i]))^2)
    R2[i]  <- if (sst == 0) NA else ssr / sst
    MSE[i] <- mean((trueK[, i] - est_matched_signed[, i])^2)
  }
  
  df_summary <- data.frame(
    true_factor = colnames(trueK) %||% paste0("True", seq_len(n_f)),
    matched_estimated_index = cols,
    pearson_corr = matched_corrs,
    abs_corr = abs(matched_corrs),
    sign = matched_signs,
    R2 = R2,
    MSE = MSE,
    stringsAsFactors = FALSE
  )
  
  if (plots) {
    pheatmap::pheatmap(cor_mat, main = paste0(title_prefix, " Score Correlation (true × est)"),
                       cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
    
    df <- data.frame(
      true = as.vector(trueK),
      est  = as.vector(est_matched_signed),
      factor = factor(rep(seq_len(n_f), each = nrow(trueK)))
    )
    
    print(
      ggplot(df, aes(x = true, y = est)) +
        geom_point(alpha = 0.6) +
        facet_wrap(~factor, scales = "free") +
        geom_smooth(method = "lm", se = FALSE) +
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


evaluate_loadings_per_block <- function(true_loadings_list, est_loadings_list, factor_assignment, matched_signs, plots = TRUE, title_prefix = "") {
  blocks <- names(true_loadings_list)
  res_block <- list()
  
  for (b in blocks) {
    true_mat <- as.matrix(true_loadings_list[[b]])
    est_mat  <- as.matrix(est_loadings_list[[b]])
    
    message("Block: ", b, " | true dims: ", paste(dim(true_mat), collapse = "x"),
            " | est dims: ", paste(dim(est_mat), collapse = "x"))
    
    # align rows if possible (by name); if not, fall back to index-based alignment
    if (!is.null(rownames(true_mat)) && !is.null(rownames(est_mat)) && all(rownames(true_mat) %in% rownames(est_mat))) {
      est_mat <- est_mat[rownames(true_mat), , drop = FALSE]
    } else if (nrow(true_mat) == nrow(est_mat)) {
      warning("Row names differ but feature counts match; using index-based alignment for block: ", b)
    } else {
      stop("Feature names/count mismatch for block ", b, " (true vs estimated).")
    }
    
    # trim factor_assignment if it references est factors that don't exist
    if (max(factor_assignment, na.rm = TRUE) > ncol(est_mat)) {
      warning("factor_assignment refers to factor indices > available estimated factors in block ", b,
              ". Will trim assignment to available estimated factors.")
      valid_idx <- which(factor_assignment <= ncol(est_mat))
      factor_assignment_use <- factor_assignment[valid_idx]
      matched_signs_use <- matched_signs[valid_idx]
    } else {
      factor_assignment_use <- factor_assignment
      matched_signs_use <- matched_signs
    }
    
    # build matched estimated loadings for the assignment
    est_matched <- est_mat[, factor_assignment_use, drop = FALSE] %*% diag(matched_signs_use)
    colnames(est_matched) <- paste0("Est_matched_", seq_len(ncol(est_matched)))
    
    # comparable number of factors
    K_true_block <- if (!is.null(ncol(true_mat))) ncol(true_mat) else 1
    K_est_block  <- if (!is.null(ncol(est_matched))) ncol(est_matched) else 1
    K <- min(K_true_block, K_est_block)
    
    if (K == 0) {
      warning("No comparable factor columns for block ", b, "; skipping.")
      res_block[[b]] <- list(summary = data.frame(), est_matched = est_matched)
      next
    }
    
    # compute per-factor correlations and RMSEs (up to K)
    corrs <- numeric(K)
    rmse  <- numeric(K)
    for (k in seq_len(K)) {
      a <- true_mat[, k]
      bvec <- est_matched[, k]
      if (all(is.na(a)) || all(is.na(bvec))) {
        corrs[k] <- NA
        rmse[k] <- NA
      } else {
        corrs[k] <- suppressWarnings(cor(a, bvec, use = "pairwise.complete.obs"))
        rmse[k]  <- sqrt(mean((a - bvec)^2, na.rm = TRUE))
      }
    }
    
    df <- data.frame(
      block = b,
      factor = paste0("Factor", seq_len(K)),
      pearson_corr = corrs,
      RMSE = rmse,
      stringsAsFactors = FALSE
    )
    
    if (plots) {
      cm <- tryCatch(cor(true_mat, est_matched, use = "pairwise.complete.obs"),
                     error = function(e) matrix(NA, nrow = ncol(true_mat), ncol = ncol(est_matched)))
      
      pheatmap(cm, main = paste0(title_prefix, " Loadings corr (true × est) block=", b),
               cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = FALSE)
      
      long_df <- bind_rows(lapply(seq_len(K), function(k) {
        data.frame(
          feature = rownames(true_mat),
          true_loading = true_mat[, k],
          est_loading  = est_matched[, k],
          factor = paste0("Factor", k),
          stringsAsFactors = FALSE
        )
      }), .id = NULL)
      
      print(
        ggplot(long_df, aes(x = true_loading, y = est_loading)) +
          geom_point(alpha = 0.4, size = 0.6) +
          facet_wrap(~factor, scales = "free") +
          geom_smooth(method = "lm", se = FALSE) +
          labs(title = paste0(title_prefix, " True vs Est Loadings (", b, ")"))
      )
    }
    
    res_block[[b]] <- list(summary = df, est_matched = est_matched)
  }
  
  res_block
}


# -------------------- Run MFA and Evaluation --------------------
# Preprocess simX_list to scaled matrices and name blocks
simX_list_clean <- lapply(simX_list, function(X) scale(as.matrix(X)))
block_names <- names(simX_list_clean) %||% paste0("omic", seq_along(simX_list_clean))

for (i in seq_along(simX_list_clean)) {
  colnames(simX_list_clean[[i]]) <- paste0(block_names[i], "_", seq_len(ncol(simX_list_clean[[i]])))
}

combined_df <- do.call(cbind, simX_list_clean)
group_sizes <- sapply(simX_list_clean, ncol)

# Fit MFA: set number of components to compute
ncp_fit <- 5
mfa_res <- MFA(
  combined_df,
  group = group_sizes,
  type = rep("s", length(group_sizes)),
  ncp = ncp_fit,
  name.group = block_names,
  graph = FALSE
)

# Extract MFA scores and loadings
mfa_scores_mat <- as.matrix(mfa_res$ind$coord)           # samples x components
mfa_loadings_full <- as.matrix(mfa_res$quanti.var$coord) # features x components


# Re-split full loadings per block
cuts <- cumsum(group_sizes)
mfa_loadings_per_block <- lapply(seq_along(group_sizes), function(i) {
  start <- if (i == 1) 1 else (cuts[i - 1] + 1)
  end   <- cuts[i]
  mfa_loadings_full[start:end, , drop = FALSE]
})
names(mfa_loadings_per_block) <- block_names

par(mfrow = c(2, 2))
plot(mfa_res$ind$coord, main="MFA scores", cex=0.3)
plot(mfa_loadings_per_block$omic1, main="MFA loadings omic1", cex=0.3)
plot(mfa_loadings_per_block$omic2, main="MFA loadings omic2", cex=0.3)
par(mfrow = c(1, 1))

# Determine safe number of factors to evaluate
n_available_mfa_comp <- ncol(mfa_scores_mat)
n_available_mfa_load_comp <- ncol(mfa_loadings_full)
n_true_factors <- ncol(true_scores)

n_f_mfa <- min(n_true_factors, n_available_mfa_comp, n_available_mfa_load_comp)
message("MFA: evaluating ", n_f_mfa, " factors (min of true factors, MFA scores, MFA loadings)")

# Evaluate MFA scores
mfa_score_eval <- evaluate_scores(true_scores, mfa_scores_mat,
                                  n_factors_eval = n_f_mfa,
                                  plots = TRUE, title_prefix = "MFA")

cat("\nMFA score summary:\n")
print(mfa_score_eval$summary)
cat("MFA mean abs corr:", mfa_score_eval$mean_abs_cor, "\n")
cat("MFA prop |corr|>=0.7:", mfa_score_eval$prop_abs_cor_ge_0.7, "\n")

# Prepare per-block loadings for evaluation (subset to n_f_mfa columns if available)
mfa_loadings_for_eval <- lapply(mfa_loadings_per_block, function(mat) {
  ncols <- ncol(mat)
  if (ncols < n_f_mfa) {
    warning("Block has fewer MFA components (", ncols, ") than n_f_mfa (", n_f_mfa, "). Will evaluate up to available columns.")
    return(mat[, seq_len(ncols), drop = FALSE])
  } else {
    return(mat[, seq_len(n_f_mfa), drop = FALSE])
  }
})

# -------------------- Align true_loadings names with MFA variable names --------------------
for (i in seq_along(true_loadings)) {
  block <- names(true_loadings)[i]
  nfeat <- nrow(true_loadings[[i]])
  new_names <- paste0(block, "_", seq_len(nfeat))
  message("Setting rownames for true_loadings[[", block, "]] to first/last: ", new_names[1], " ... ", new_names[length(new_names)])
  rownames(true_loadings[[i]]) <- new_names
}

# Diagnostics: check first few names
for (i in seq_along(mfa_loadings_per_block)) {
  block <- names(mfa_loadings_per_block)[i]
  rn_true <- rownames(true_loadings[[block]])
  rn_est  <- rownames(mfa_loadings_per_block[[block]])
  message("Block ", block, ": true first name = ", if (length(rn_true) > 0) rn_true[1] else "NULL",
          " | est first name = ", if (length(rn_est) > 0) rn_est[1] else "NULL")
}

# -------------------- Evaluate MFA loadings (per-block) --------------------
mfa_load_eval <- evaluate_loadings_per_block(
  true_loadings,
  mfa_loadings_for_eval,
  factor_assignment = mfa_score_eval$assignment,
  matched_signs = mfa_score_eval$matched_signs,
  plots = TRUE,
  title_prefix = "MFA"
)

# Print concise summaries
cat("\nMFA loadings (per block) summary after name-alignment:\n")
for (b in names(mfa_load_eval)) {
  print(head(mfa_load_eval[[b]]$summary, 10))
}


# -------------------- Heatmaps --------------------
# MFA score correlation heatmap
cor_mat_mfa <- mfa_score_eval$correlation_matrix
rownames(cor_mat_mfa) <- paste0("True_F", seq_len(nrow(cor_mat_mfa)))
colnames(cor_mat_mfa) <- paste0("MFA_F", seq_len(ncol(cor_mat_mfa)))

pheatmap::pheatmap(
  cor_mat_mfa,
  main = "MFA Score Correlation (True × Estimated)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 10
)

# MFA loadings correlation heatmap per block
for (block_name in names(mfa_load_eval)) {
  block_res <- mfa_load_eval[[block_name]]
  if (is.null(block_res$est_matched) || ncol(block_res$est_matched) == 0) {
    warning("No matched factors for block: ", block_name, "; skipping heatmap.")
    next
  }
  
  true_mat <- as.matrix(true_loadings[[block_name]])
  est_mat  <- as.matrix(block_res$est_matched)
  
  # Align rownames if needed
  if (!is.null(rownames(true_mat)) && !is.null(rownames(est_mat)) && all(rownames(true_mat) %in% rownames(est_mat))) {
    est_mat <- est_mat[rownames(true_mat), , drop = FALSE]
  } else if (nrow(true_mat) != nrow(est_mat)) {
    warning("Row mismatch for block ", block_name, "; skipping heatmap.")
    next
  }
  
  cor_block <- cor(true_mat, est_mat, use = "pairwise.complete.obs")
  rownames(cor_block) <- paste0("True_F", seq_len(nrow(cor_block)))
  colnames(cor_block) <- paste0("Est_F", seq_len(ncol(cor_block)))
  
  pheatmap::pheatmap(
    cor_block,
    main = paste0("MFA Loadings Correlation (Block=", block_name, ")"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 8
  )
}
