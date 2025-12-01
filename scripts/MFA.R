# ==============================================================================
# Robust MFA Block (with Diagnostic Plots)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Helpers & Packages
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y

library(FactoMineR) # for MFA
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clue)
library(fabia)

# ------------------------------------------------------------------------------
# 1. Simulated Data Generation
# ------------------------------------------------------------------------------
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

# Combine simulations
sim_object <- sim1
sim_object$omics$omic2 <- sim2$omics$omic2
sim_object[["list_betas"]][[2]] <- sim2[["list_betas"]][[2]]

# Final data list
simX_list <- sim_object$omics
names(simX_list) <- paste0("omic", seq_along(simX_list))

# --- PLOT: Raw Simulation Heatmap ---
plot_simData(sim_object = sim_object, type = "heatmap")


# ------------------------------------------------------------------------------
# 2. Extract Truth & Visualize True Structure
# ------------------------------------------------------------------------------

# --- True Scores ---
true_scores <- do.call(cbind, lapply(sim_object$list_alphas, as.numeric))
colnames(true_scores) <- names(sim_object$list_alphas)
rownames(true_scores) <- rownames(sim_object$omics[[1]]) %||% seq_len(nrow(true_scores))
n_factors <- ncol(true_scores)

# --- PLOT: True Factor Scores ---
df_scores_true <- as.data.frame(true_scores) %>%
  mutate(Sample = rownames(true_scores)) %>%
  pivot_longer(cols = -Sample, names_to = "Factor", values_to = "Score")

print(ggplot(df_scores_true, aes(x = Sample, y = Score, color = Factor, group = Factor)) +
        geom_point(alpha = 0.6) +
        facet_wrap(~Factor, scales = "free_y") +
        labs(title = "True Factor Scores (Simulated Data)", x = "Sample", y = "Score") +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))

# --- True Loadings ---
get_true_loadings_per_block <- function(sim_object) {
  true_loadings_list <- list()
  for (i in seq_along(sim_object$list_betas)) {
    block_betas <- sim_object$list_betas[[i]]
    block_mats <- lapply(block_betas, function(b) matrix(as.numeric(b), ncol = 1))
    block_combined <- do.call(cbind, block_mats)
    
    # Naming convention: omicX_feature_Y
    rownames(block_combined) <- paste0("omic", i, "_feature_", seq_len(nrow(block_combined)))
    true_loadings_list[[i]] <- block_combined
  }
  names(true_loadings_list) <- paste0("omic", seq_along(sim_object$list_betas))
  true_loadings_list
}

true_loadings <- get_true_loadings_per_block(sim_object)

# --- PLOT: True Loadings ---
loadings_long <- lapply(names(true_loadings), function(block) {
  mat <- true_loadings[[block]]
  df <- as.data.frame(mat)
  df$Feature <- rownames(mat)
  df$Block <- block
  pivot_longer(df, cols = -c(Feature, Block), names_to = "Factor", values_to = "Loading")
})
df_loadings_true <- bind_rows(loadings_long)

print(ggplot(df_loadings_true, aes(x = Feature, y = Loading, color = Factor)) +
        geom_point(alpha = 0.4, size = 0.5) +
        facet_wrap(~Block + Factor, scales = "free_y") +
        labs(title = "True Loadings per Block (Simulated Data)", y = "Loading") +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))

# --- PLOT: Base R Quick Check ---
par(mfrow = c(2, 2))
plot(sim_object[["list_betas"]][[1]][["beta2"]], main = "True Beta2 (Omic1)", ylab = "Loadings", cex = 0.4)
plot(sim_object[["list_betas"]][[2]][["beta1"]], main = "True Beta1 (Omic2)", ylab = "Loadings", cex = 0.4)
plot(sim_object[["list_betas"]][[2]][["beta2"]], main = "True Beta2 (Omic2)", ylab = "Loadings", cex = 0.4)
plot(sim_object[["list_betas"]][[2]][["beta3"]], main = "True Beta3 (Omic2)", ylab = "Loadings", cex = 0.4)
par(mfrow = c(1, 1))


# ------------------------------------------------------------------------------
# 3. Evaluation Functions (Logic + Integrated Plotting)
# ------------------------------------------------------------------------------

evaluate_scores <- function(true_mat, est_mat, n_factors_eval = NULL, plots = TRUE, title_prefix = "") {
  true_mat <- as.matrix(true_mat)
  est_mat  <- as.matrix(est_mat)
  
  if (nrow(true_mat) != nrow(est_mat)) stop("Sample rows mismatch.")
  
  K_true <- ncol(true_mat)
  K_est  <- ncol(est_mat)
  n_f <- n_factors_eval %||% min(K_true, K_est)
  
  trueK <- true_mat[, seq_len(n_f), drop = FALSE]
  estK  <- est_mat[, seq_len(min(K_est, n_f)), drop = FALSE]
  
  # 1. Match Factors (Hungarian Algorithm on Correlation)
  cor_mat <- cor(trueK, estK, use = "pairwise.complete.obs")
  cost <- 1 - abs(ifelse(is.na(cor_mat), 0, cor_mat))
  cols <- as.integer(clue::solve_LSAP(as.matrix(cost)))
  
  # 2. Compute Metrics & Fix Signs
  matched_corrs <- numeric(n_f)
  matched_signs <- numeric(n_f)
  
  for (i in seq_len(n_f)) {
    j <- cols[i]
    cval <- if(is.na(cor_mat[i, j])) 0 else cor_mat[i, j]
    matched_corrs[i] <- cval
    matched_signs[i] <- if(sign(cval) == 0) 1 else sign(cval)
  }
  
  est_matched_signed <- estK[, cols, drop = FALSE] %*% diag(matched_signs)
  colnames(est_matched_signed) <- paste0("Est_matched_", seq_len(n_f))
  
  # 3. Visualizations (Integrated)
  if (plots) {
    pheatmap::pheatmap(cor_mat,
                       main = paste0(title_prefix, " Score Correlation (True vs Est)"),
                       cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
                       number_color = "black", fontsize_number = 10)
    
    df_plot <- data.frame(
      true = as.vector(trueK),
      est  = as.vector(est_matched_signed),
      factor = factor(rep(seq_len(n_f), each = nrow(trueK)))
    )
    print(ggplot(df_plot, aes(x = true, y = est)) +
            geom_point(alpha = 0.6) +
            facet_wrap(~factor, scales = "free") +
            geom_smooth(method = "lm", se = FALSE) +
            labs(title = paste0(title_prefix, " True vs Estimated Scores")))
  }
  
  list(
    correlation_matrix = cor_mat,
    assignment = cols,
    matched_signs = matched_signs,
    summary = data.frame(
      true_factor = colnames(trueK), matched_idx = cols,
      pearson = matched_corrs, abs_corr = abs(matched_corrs)
    ),
    mean_abs_cor = mean(abs(matched_corrs), na.rm = TRUE),
    est_signed_matched = est_matched_signed
  )
}

evaluate_loadings_per_block <- function(true_loadings_list, est_loadings_list, factor_assignment, matched_signs, plots = TRUE, title_prefix = "") {
  blocks <- names(true_loadings_list)
  res_block <- list()
  
  for (b in blocks) {
    true_mat <- as.matrix(true_loadings_list[[b]])
    est_mat  <- as.matrix(est_loadings_list[[b]])
    
    # Align rows
    if (!is.null(rownames(true_mat)) && !is.null(rownames(est_mat)) && all(rownames(true_mat) %in% rownames(est_mat))) {
      est_mat <- est_mat[rownames(true_mat), , drop = FALSE]
    } else if(nrow(true_mat) != nrow(est_mat)) {
      stop("Feature count mismatch in block ", b)
    }
    
    # Apply matching from Score Evaluation
    valid_idx <- which(factor_assignment <= ncol(est_mat))
    factor_assignment_use <- factor_assignment[valid_idx]
    matched_signs_use <- matched_signs[valid_idx]
    
    est_matched <- est_mat[, factor_assignment_use, drop = FALSE] %*% diag(matched_signs_use)
    colnames(est_matched) <- paste0("Est_matched_", seq_len(ncol(est_matched)))
    
    K <- min(ncol(true_mat), ncol(est_matched))
    
    if (plots && K > 0) {
      cm <- cor(true_mat, est_matched, use = "pairwise.complete.obs")
      pheatmap(cm,
               main = paste0(title_prefix, " Loadings Corr (Block=", b, ")"),
               cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
               number_color = "black", fontsize_number = 8)
      
      long_df <- bind_rows(lapply(seq_len(K), function(k) {
        data.frame(
          true_loading = true_mat[, k],
          est_loading  = est_matched[, k],
          factor = paste0("Factor", k)
        )
      }))
      print(ggplot(long_df, aes(x = true_loading, y = est_loading)) +
              geom_point(alpha = 0.4, size = 0.6) +
              facet_wrap(~factor, scales = "free") +
              geom_smooth(method = "lm", se = FALSE) +
              labs(title = paste0(title_prefix, " Loadings Scatter (", b, ")")))
    }
    
    corrs <- sapply(seq_len(K), function(k) cor(true_mat[,k], est_matched[,k], use="pair"))
    res_block[[b]] <- list(summary = data.frame(block=b, factor=paste0("F",1:K), pearson=corrs), est_matched = est_matched)
  }
  res_block
}

# ------------------------------------------------------------------------------
# 4. Execution: Run MFA
# ------------------------------------------------------------------------------

# Scale and naming
simX_list_clean <- lapply(simX_list, function(X) scale(as.matrix(X)))
block_names <- names(simX_list_clean)
for (i in seq_along(simX_list_clean)) {
  colnames(simX_list_clean[[i]]) <- paste0(block_names[i], "_", seq_len(ncol(simX_list_clean[[i]])))
}

combined_df <- do.call(cbind, simX_list_clean)
group_sizes <- sapply(simX_list_clean, ncol)

# Fit MFA
mfa_res <- MFA(combined_df, group = group_sizes, type = rep("s", length(group_sizes)),
               ncp = 5, name.group = block_names, graph = FALSE)

# Extract matrices
mfa_scores_mat <- as.matrix(mfa_res$ind$coord)
mfa_loadings_full <- as.matrix(mfa_res$quanti.var$coord)

# Split loadings back to blocks
cuts <- cumsum(group_sizes)
mfa_loadings_per_block <- lapply(seq_along(group_sizes), function(i) {
  start <- if (i == 1) 1 else (cuts[i - 1] + 1)
  mfa_loadings_full[start:cuts[i], , drop = FALSE]
})
names(mfa_loadings_per_block) <- block_names

# --- PLOT: MFA Diagnostics ---
par(mfrow = c(2, 2))
plot(mfa_res$ind$coord, main = "MFA Scores (Samples)", cex = 0.3)
plot(mfa_loadings_per_block$omic1, main = "MFA Loadings (Omic1)", cex = 0.3)
plot(mfa_loadings_per_block$omic2, main = "MFA Loadings (Omic2)", cex = 0.3)
par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# 5. Execution: Evaluation (Math & Plots)
# ------------------------------------------------------------------------------

n_f_mfa <- min(ncol(true_scores), ncol(mfa_scores_mat))

message("--- Evaluating Scores ---")
mfa_score_eval <- evaluate_scores(
  true_scores, 
  mfa_scores_mat,
  n_factors_eval = n_f_mfa,
  plots = TRUE, 
  title_prefix = "MFA"
)

cat("MFA Score Summary:\n")
print(mfa_score_eval$summary)

# Update True Loadings Names to match MFA output names (required for alignment)
for (i in seq_along(true_loadings)) {
  rn <- paste0(names(true_loadings)[i], "_", seq_len(nrow(true_loadings[[i]])))
  rownames(true_loadings[[i]]) <- rn
}

mfa_loadings_eval_input <- lapply(mfa_loadings_per_block, function(x) x[, 1:n_f_mfa, drop=FALSE])

message("\n--- Evaluating Loadings ---")
mfa_load_eval <- evaluate_loadings_per_block(
  true_loadings,
  mfa_loadings_eval_input,
  factor_assignment = mfa_score_eval$assignment,
  matched_signs = mfa_score_eval$matched_signs,
  plots = TRUE,
  title_prefix = "MFA"
)

cat("MFA Loadings Summary (Block 1):\n")
print(head(mfa_load_eval[[1]]$summary))
