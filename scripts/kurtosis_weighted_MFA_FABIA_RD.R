# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================

# Combined list of all required packages from both sections
required_pkgs <- c("FactoMineR", "fabia", "clue", "pheatmap", 
                   "ggplot2", "matrixStats", "reshape2", "e1071", "gtools", 
                   "gridExtra", "dplyr")

to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install)) install.packages(to_install)

library(FactoMineR)
library(fabia)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(gtools)      # For natural sorting (mixedsort)
library(matrixStats) # Efficient matrix operations
library(reshape2)
library(e1071)       # For Kurtosis calculation
library(gridExtra)


# ==============================================================================
# 2. DATA LOADING & PREPROCESSING
# ==============================================================================

message("--- Step 1: Loading & Aligning Data ---")

# Load Metadata
D17_metadata <- read.csv(file = "../multi-omics radiation_data/D17_metadata.csv",
                         header = TRUE, row.names = 1, check.names = TRUE)

# Load Omics Data (Input: Features as Rows, Samples as Columns)
D17_mRNA_data <- read.csv(file = "../multi-omics radiation_data/D17_mRNA_cor_normalized_tmm.csv",
                          header = TRUE, row.names = 1, check.names = TRUE)
D17_protein_data <- read.csv(file = "../multi-omics radiation_data/D17_Protein_cor_normalized_quant.csv",
                             header = TRUE, row.names = 1, check.names = TRUE)

# Transpose to Standard format (Samples as Rows, Features as Columns)
mRNA_mat <- t(as.matrix(D17_mRNA_data))
prot_mat <- t(as.matrix(D17_protein_data))

# Align Samples
common_samples <- intersect(rownames(mRNA_mat), rownames(prot_mat))
if(length(common_samples) == 0) stop("No matching samples found!")

# Subset to common samples
X_list <- list(
  mRNA = mRNA_mat[common_samples, ],
  Protein = prot_mat[common_samples, ]
)

message(sprintf("Aligned %d Samples across %d mRNA and %d Protein features.", 
                length(common_samples), ncol(X_list$mRNA), ncol(X_list$Protein)))

# ==============================================================================
# 3. VISUALIZATION: Global Heatmap
# ==============================================================================

message("--- Step 2: Visualizing Global Structure ---")

# Helper to clean and scale
clean_and_scale <- function(M) {
  vars <- matrixStats::colVars(M)
  M <- M[, vars > 1e-9] # Remove constant features
  s <- scale(M)
  s[is.na(s)] <- 0
  return(s)
}

# Process for Heatmap
Z_vis <- cbind(clean_and_scale(X_list$mRNA), clean_and_scale(X_list$Protein))

# Annotation Track
annotation_df <- data.frame(
  Type = factor(rep(c("mRNA", "Protein"),
                    times = c(ncol(clean_and_scale(X_list$mRNA)),
                              ncol(clean_and_scale(X_list$Protein))))))
rownames(annotation_df) <- colnames(Z_vis)

# Order Samples naturally (X13, X15...)
sorted_samples <- gtools::mixedsort(rownames(Z_vis))
Z_vis_ordered <- t(Z_vis[sorted_samples, ])

# Render Global Heatmap
pheatmap(Z_vis_ordered,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = seq(-3, 3, length.out = 100),
         show_rownames = FALSE, show_colnames = TRUE,
         cluster_cols = FALSE, cluster_rows = FALSE, # Manual Order
         annotation_row = annotation_df,
         main = "Global Multi-Omics Data Structure")


# ==============================================================================
# 4. CORE ALGORITHM: Kurtosis-Weighted MFA-FABIA (Updated Robust Version)
# ==============================================================================

message("--- Step 3: Running Hybrid MFA-FABIA ---")

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
                               seed = 123) {
  
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
  
  return(list(
    model              = fab_fit,        # Kept for backward compat with your visualization code
    fabia              = fab_fit,        # Kept for consistency with newer code
    factor_scores      = factor_scores,
    loadings           = loadings_all,
    loadings_per_block = loadings_per_block,
    weights            = w_all,
    kurtosis           = unlist(kurtosis_per_block)
  ))
}

# Run the Algorithm (Updated Parameters)
# Note: Using block_equalize_method="trace" for better stability than simple eigenvalue scaling
res_hybrid <- mfa_weighted_fabia(X_list, 
                                 gamma = 0.5, 
                                 p = 4, 
                                 alpha = 0.05, 
                                 scale_vars = TRUE,
                                 block_equalize = TRUE,
                                 block_equalize_method = "trace", 
                                 seed = 123)

# ==============================================================================
# VISUALIZATION: Annotated Factor Scores
# ==============================================================================

message("Generating annotated factor score heatmap")

# 1. load metadata
if(!exists("D17_metadata")) {
  D17_metadata <- read.csv(file = "../multi-omics radiation_data/D17_metadata.csv",
                           header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

# 2. Prepare Annotation Dataframe
heatmap_data <- res_hybrid$factor_scores # ((Samples x Factors))

# Clean "X" prefix to match metadata IDs (e.g. "X13" -> 13)
clean_ids <- as.integer(sub("^X", "", rownames(heatmap_data)))

# Match IDs to find the correct metadata row for each sample
meta_idx <- match(clean_ids, D17_metadata$Internal.ID)

# Create the annotation table
annot_df <- data.frame(
  Description = as.factor(D17_metadata$Description[meta_idx]),
  Treatment   = as.factor(D17_metadata$Treatment[meta_idx]),
  Control     = as.factor(D17_metadata$Control[meta_idx])
)

# Row names of annotation must match heatmap row names
rownames(annot_df) <- rownames(heatmap_data)

# 3. Define Custom Colors 
ann_colors <- list(
  Treatment = c("TRUE" = "#E41A1C", "FALSE" = "#F0F0F0"), # Red for Treated
  Control   = c("TRUE" = "#4DAF4A", "FALSE" = "#F0F0F0")  # Green for Control
)

pheatmap(heatmap_data,
         cluster_cols = FALSE,      # Don't reorder Factors (keep 1, 2, 3, 4)
         cluster_rows = TRUE,       # Cluster Samples to see which group together
         display_numbers = TRUE,    # Show score values
         annotation_row = annot_df,      
         annotation_colors = ann_colors, 
         main = "Factor Scores (sample activity)")



# ==============================================================================
# 5. SCATTER PLOTS (for SAMPLE SCORES)
#           # 👉 How strongly a sample expresses a latent factor)
# ==============================================================================

message("--- Step 4: Generating Scatter Plots ---")

scores_df <- as.data.frame(res_hybrid$factor_scores)
colnames(scores_df) <- paste0("Factor", 1:ncol(scores_df))

# Clean ID: Remove "X" prefix and ensure integer match
scores_df$ID_clean <- as.integer(sub("^X", "", rownames(scores_df)))

# Merge
scores_merged <- merge(
  x = scores_df, 
  y = D17_metadata, 
  by.x = "ID_clean", 
  by.y = "Internal.ID", 
  all.x = TRUE
)

# Add Label back
scores_merged$SampleLabel <- rownames(scores_df)[match(scores_merged$ID_clean, scores_df$ID_clean)]

# convert groups to factors
scores_merged$Description <- as.factor(scores_merged$Description)
scores_merged$Control     <- as.factor(scores_merged$Control)
scores_merged$Treatment   <- as.factor(scores_merged$Treatment)

# --- B. Define Theme ---
my_theme <- theme_minimal() +
  theme(plot.title = element_text(face="bold", size=12),
        legend.position = "right")

# Check if merge worked
if(nrow(scores_merged) == 0) {
  stop("Merge failed! No common IDs")
} else {
  message(sprintf("Successfully merged metadata for %d samples.", nrow(scores_merged)))
}

# ==============================================================================
# PART 1: Standard Scatter Plots (Factor 1 vs Factor 2)
# ==============================================================================

# 1) Scatter Plot: NO GROUPING
p1 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.7, color = "black") + 
  geom_text(vjust = -0.8, size = 3) +
  labs(title = "1. Scores (ungrouped)", x="Factor 1", y="Factor 2") +
  my_theme

# 2) Scatter Plot: Grouped by 'Description'
p2 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, color = Description, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "2. Scores by 'Description'", x="Factor 1", y="Factor 2") +
  my_theme

# 3) Scatter Plot: Grouped by 'Control'
p3 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, color = Control, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "3. Scores by 'Control'", x="Factor 1", y="Factor 2") +
  my_theme

# 4) Scatter Plot: Grouped by 'Treatment'
p4 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, color = Treatment, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "4. Scores by 'Treatment'", x="Factor 1", y="Factor 2") +
  my_theme

# Display Part 1
grid.arrange(p1, p2, p3, p4, nrow = 2)


# 1) Scatter Plot: NO GROUPING
p11 <- ggplot(scores_merged, aes(x = Factor3, y = Factor4, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.7, color = "black") + 
  geom_text(vjust = -0.8, size = 3) +
  labs(title = "1. Scores (ungrouped)", x="Factor 3", y="Factor 4") +
  my_theme

# 2) Scatter Plot: Grouped by 'Description'
p21 <- ggplot(scores_merged, aes(x = Factor3, y = Factor4, color = Description, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "2. Scores by 'Description'", x="Factor 3", y="Factor 4") +
  my_theme

# 3) Scatter Plot: Grouped by 'Control'
p31 <- ggplot(scores_merged, aes(x = Factor3, y = Factor4, color = Control, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "3. Scores by 'Control'", x="Factor 3", y="Factor 4") +
  my_theme

# 4) Scatter Plot: Grouped by 'Treatment'
p41 <- ggplot(scores_merged, aes(x = Factor3, y = Factor4, color = Treatment, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "4. Scores by 'Treatment'", x="Factor 3", y="Factor 4") +
  my_theme

# Display Part 1
grid.arrange(p11, p21, p31, p41, nrow = 2)


# ==============================================================================
# PART 2: Factor Scan Plots (All Factors vs Samples)
# Plots for factor scores v samples; for each factor
# ==============================================================================

long_scores <- melt(scores_merged,
                    id.vars = c("SampleLabel", "Description", "Control", "Treatment"),
                    measure.vars = colnames(scores_df)[grep("Factor", colnames(scores_df))],
                    variable.name = "Factor", 
                    value.name = "Score")

# Define Plotting Function
plot_factor_scan <- function(data, group_col, title_text) {
  
  # 1. Determine Sample Order based on UNIQUE samples
  # (Sort by Group first, then by Sample Name)
  unique_meta <- unique(data[, c("SampleLabel", group_col)])
  unique_meta <- unique_meta[order(unique_meta[[group_col]], unique_meta$SampleLabel), ]
  
  # 2. Apply this unique order
  data$OrderedSample <- factor(data$SampleLabel, levels = unique_meta$SampleLabel)
  
  # 3. Plot
  ggplot(data, aes(x = OrderedSample, y = Score, color = .data[[group_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_segment(aes(x = OrderedSample, xend = OrderedSample, y = 0, yend = Score), alpha=0.5) +
    geom_point(size = 3, alpha = 0.9) +
    
    # Facet by Factor (Show F1, F2, F3, F4)      
    facet_wrap(~Factor, ncol = 2, scales = "free_y") +
    
    labs(title = title_text,
         subtitle = "Ordered by Group",
         x = "Samples", 
         y = "Factor Score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8),
          strip.text = element_text(face = "bold", size = 10))
}

# 5) Scan by 'Control'
p5 <- plot_factor_scan(long_scores, "Control", "5. Factor Scan: Control")

# 6) Scan by 'Treatment'
p6 <- plot_factor_scan(long_scores, "Treatment", "6. Factor Scan: Treatment")

# 7) Scan by 'Description'
p7 <- plot_factor_scan(long_scores, "Description", "7. Factor Scan: Description")

# Display Part 2
print(p5)
print(p6)
print(p7)


# ==============================================================================
# 6. EXPLORATORY VISUALIZATIONS FOR LOADINGS
# ==============================================================================


message("Generating Exploratory Loading Plots...")

# --- PREP: Create a master dataframe of all loadings ---
all_loadings <- as.data.frame(res_hybrid$loadings)
colnames(all_loadings) <- paste0("Factor", 1:ncol(all_loadings))
all_loadings$FeatureID <- rownames(all_loadings)
# Identify Omic Type
all_loadings$Type <- ifelse(grepl("^mRNA_", all_loadings$FeatureID), "mRNA", "Protein")
# Clean Feature Name (remove prefix for display)
all_loadings$GeneName <- sub("^mRNA_|^Protein_", "", all_loadings$FeatureID)

# ==============================================================================
# VISUALIZATION A: Ranked Lollipop Plot (Top 20 Features per Factor)
# ==============================================================================
# Goal: Explore the top biological drivers for each factor.

plot_list <- list()

for(i in 1:ncol(res_hybrid$factor_scores)) {
  factor_col <- paste0("Factor", i)
  
  # 1. Select Top 20 features by absolute loading weight
  top_data <- all_loadings %>%
    arrange(desc(abs(get(factor_col)))) %>%
    head(20)
  
  # 2. Create Lollipop Plot
  p <- ggplot(top_data, aes(x = reorder(GeneName, abs(get(factor_col))), 
                            y = get(factor_col), 
                            color = Type)) +
    geom_segment(aes(x = GeneName, xend = GeneName, y = 0, yend = get(factor_col)), color="grey") +
    geom_point(size = 3) +
    coord_flip() + # Flip to make labels readable
    scale_color_manual(values = c("mRNA" = "#1f77b4", "Protein" = "#ff7f0e")) +
    labs(title = paste("Top Drivers: Factor", i),
         y = "Loading Weight", x = NULL) +
    theme_minimal() +
    theme(legend.position = "none") # Hide legend to save space (add back if needed)
  
  plot_list[[i]] <- p
}

# Arrange in a grid
grid.arrange(grobs = plot_list, ncol = 2, top = "Top 20 Features per Factor (Ranked)")


# ==============================================================================
# VISUALIZATION B: Clustered Heatmap of Top Features
### LOADING HEATMAPS (Global & Per-Omic)
# ==============================================================================
# Goal: Check specificity. Do Factor 1 genes also light up in Factor 2?

message("Generating Unified Loading Heatmaps...")

#' Plot a Clustered Heatmap of Top Loadings
#' 
#' @param loading_matrix Matrix of feature loadings (Features x Factors)
#' @param title_text Main title for the plot
#' @param top_n Number of top features per factor to select
#' @param annotate_type Boolean. If TRUE, adds an "mRNA/Protein" annotation bar 
#'         and hides row names (best for Global plots). If FALSE, strips prefixes 
#'         and shows row names (best for Per-Omic plots).
plot_loading_heatmap <- function(loading_matrix, title_text, top_n = 50, annotate_type = FALSE) {
  
  # 1. Identify Top N features for EACH factor
  # Loop through columns (Factors) and find indices of highest absolute weights
  top_indices <- unique(unlist(lapply(1:ncol(loading_matrix), function(i) {
    order(abs(loading_matrix[, i]), decreasing = TRUE)[1:top_n]
  })))
  
  # 2. Subset the matrix
  subset_mat <- loading_matrix[top_indices, , drop = FALSE]
  
  # 3. Configure Annotation & Row Names
  annot_row <- NA
  annot_col <- NA
  show_rows <- TRUE
  
  if(annotate_type) {
    # --- Global Mode ---
    # Create "Type" annotation based on row prefixes
    types <- ifelse(grepl("^mRNA_", rownames(subset_mat)), "mRNA", "Protein")
    annot_row <- data.frame(Type = factor(types))
    rownames(annot_row) <- rownames(subset_mat)
    annot_col <- list(Type = c(mRNA = "#1f77b4", Protein = "#ff7f0e"))
    
    # Hide row names because 200+ mixed genes/proteins are too crowded
    show_rows <- FALSE 
  } else {
    # --- Per-Omic Mode ---
    # Clean prefixes (e.g., "mRNA_GeneA" -> "GeneA") for readability
    rownames(subset_mat) <- sub("^mRNA_|^Protein_", "", rownames(subset_mat))
    show_rows <- TRUE
  }
  
  # 4. Define Dynamic Color Limits
  # Scales the color palette to the data's specific range (crucial for Protein vs mRNA)
  limit <- max(abs(subset_mat))
  breaks_list <- seq(-limit, limit, length.out = 100)
  
  # 5. Render Heatmap
  pheatmap(subset_mat,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           breaks = breaks_list,
           annotation_row = annot_row,
           annotation_colors = annot_col,
           cluster_cols = FALSE,        # Keep Factors ordered (1, 2, 3...)
           cluster_rows = TRUE,         # Cluster features to show shared patterns
           show_rownames = show_rows,
           fontsize_row = 4,
           main = paste0(title_text, " (Union of Top ", top_n, " per Factor)"))
}

# --- EXECUTE PLOTS ---

# 1. Global Heatmap (Combined mRNA & Protein)
# shows broad patterns and omic-specificity (via annotation bar)
plot_loading_heatmap(res_hybrid$loadings, 
                     "Global Top Loadings", 
                     top_n = 50, 
                     annotate_type = FALSE)



# 2. mRNA Heatmap (Specific Transcriptomic Drivers)
# shows specific gene names
plot_loading_heatmap(res_hybrid$loadings_per_block$mRNA, 
                     "mRNA Top Loadings", 
                     top_n = 50, 
                     annotate_type = FALSE)

# 3. Protein Heatmap (Specific Proteomic Drivers)
# re-scaled color bar reveals subtle protein patterns invisible in the global plot
plot_loading_heatmap(res_hybrid$loadings_per_block$Protein, 
                     "Protein Top Loadings", 
                     top_n = 50, 
                     annotate_type = FALSE)


# ==============================================================================
# VISUALIZATION C: Sparsity Diagnostic (Density Plot)
# ==============================================================================
# 👉 Goal: Confirm FABIA worked. You want a high peak at 0 (noise) and "Heavy Tails" (signal).

melted_loadings <- melt(all_loadings, id.vars = c("FeatureID", "Type", "GeneName"))

ggplot(melted_loadings, aes(x = value, fill = Type)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("mRNA" = "#1f77b4", "Protein" = "#ff7f0e")) +
  facet_wrap(~variable, scales = "free") + # Facet by Factor
  xlim(-0.5, 0.5) + # Zoom in on the center to see the "tails" better
  labs(title = "Loading Distribution (Sparsity Check)",
       subtitle = "Ideal FABIA Result: Sharp peak at 0, long visible tails",
       x = "Loading Value", y = "Density") +
  theme_minimal()



# ==============================================================================
# 7. PER-OMIC LOADING EXPLORATION
# ==============================================================================


message("Generating Per-Omic Loading Plots...")

# Function to generate plots for a single omic block
visualize_block_loadings <- function(block_name, loading_matrix, top_n = 15) {
  
  # Convert matrix to dataframe
  df <- as.data.frame(loading_matrix)
  colnames(df) <- paste0("Factor", 1:ncol(df))
  df$FeatureID <- rownames(df)
  
  # Remove the prefix (e.g., "mRNA_") for cleaner labels
  prefix_pattern <- paste0("^", block_name, "_")
  df$GeneName <- sub(prefix_pattern, "", df$FeatureID)
  
  # List to store plots for this block
  plot_list <- list()
  
  # Iterate through each factor to make a Lollipop Plot
  for(i in 1:ncol(loading_matrix)) {
    factor_col <- paste0("Factor", i)
    
    # Get Top N features for this specific omic/factor combo
    top_feats <- df %>%
      arrange(desc(abs(get(factor_col)))) %>%
      head(top_n)
    
    # Create Plot
    p <- ggplot(top_feats, aes(x = reorder(GeneName, abs(get(factor_col))), 
                               y = get(factor_col))) +
      # Lollipop stick
      geom_segment(aes(x = GeneName, xend = GeneName, y = 0, yend = get(factor_col)), 
                   color = ifelse(block_name == "mRNA", "#1f77b4", "#ff7f0e"), 
                   alpha = 0.6) +
      # Lollipop head
      geom_point(size = 3, 
                 color = ifelse(block_name == "mRNA", "#1f77b4", "#ff7f0e")) +
      coord_flip() +
      labs(title = paste0(block_name, ": Factor ", i),
           x = NULL, y = "Loading Weight") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            axis.text.y = element_text(size = 8))
    
    plot_list[[i]] <- p
  }
  
  return(plot_list)
}

# --- GENERATE & DISPLAY PLOTS ---

# 1. Generate mRNA Plots
mRNA_plots <- visualize_block_loadings("mRNA", res_hybrid$loadings_per_block$mRNA, top_n = 15)

# 2. Generate Protein Plots
prot_plots <- visualize_block_loadings("Protein", res_hybrid$loadings_per_block$Protein, top_n = 15)

# --- DISPLAY 

# Option A: View mRNA Factors Side-by-Side
grid.arrange(grobs = mRNA_plots, ncol = 2, top = "Top mRNA Drivers per Factor")



# Option B: View Protein Factors Side-by-Side
grid.arrange(grobs = prot_plots, ncol = 2, top = "Top Protein Drivers per Factor")

# Option C: Comparative View (Factor 1 mRNA vs Factor 1 Protein)
# This is useful to see if the SAME biological pathway is active in both layers for Factor 1
grid.arrange(mRNA_plots[[1]], prot_plots[[1]], 
             mRNA_plots[[2]], prot_plots[[2]], 
             ncol = 2, 
             top = "Comparison: mRNA (Left) vs Protein (Right) for Factors 1 & 2")



# ==============================================================================
# 8. VARIANCE EXPLAINED (Per Factor / Per Omic)
# ==============================================================================

message("Generating Variance Percentage Table...")

# 1. Calculate Sum of Squared Loadings (SSL) per block/factor
# (Re-using logic from previous step for clarity)
ssl_matrix <- sapply(res_hybrid$loadings_per_block, function(loadings) {
  colSums(loadings^2)
})

# ssl_matrix is now:
#         mRNA  Protein
# Factor1 120.5 45.2
# Factor2 80.3  90.1
# ...

# 2. Calculate Totals per Factor (Row Sums)
total_variance_per_factor <- rowSums(ssl_matrix)

# 3. Calculate Percentages
percent_matrix <- (ssl_matrix / total_variance_per_factor) * 100

# 4. Format the Table
# Round to 2 decimal places and add "%" sign
variance_table <- as.data.frame(round(percent_matrix, 2))
colnames(variance_table) <- c("mRNA (%)", "Protein (%)")

# Add a "Dominant Layer" column for quick interpretation
variance_table$Dominant_Layer <- ifelse(variance_table$`mRNA (%)` > 50, "mRNA", "Protein")

# 5. Print the Table
print("--- Variance Explained Composition per Factor ---")
print(variance_table)

