# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================

required_pkgs <- c("FactoMineR", "fabia", "clue", "pheatmap", 
                   "ggplot2", "matrixStats", "reshape2", "e1071", "gtools")
to_install <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(to_install)) install.packages(to_install)

library(FactoMineR)
library(fabia)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(gtools)      # For natural sorting (mixedsort)
library(matrixStats)
library(reshape2)
library(e1071)       # For Kurtosis calculation
library(gridExtra)
library(reshape2)


# ==============================================================================
# 2. DATA LOADING & PREPROCESSING
# ==============================================================================

message("--- Step 1: Loading & Aligning Data ---")

# Load metadata
D17_metadata <- read.csv(file = "../multi-omics radiation_data/D17_metadata.csv",
                          header = TRUE, row.names = 1, check.names = TRUE)

# 2.1 Load Datasets (Input: Features as Rows, Samples as Columns)
D17_mRNA_data <- read.csv(file = "../multi-omics radiation_data/D17_mRNA_cor_normalized_tmm.csv",
                          header = TRUE, row.names = 1, check.names = TRUE)
D17_protein_data <- read.csv(file = "../multi-omics radiation_data/D17_Protein_cor_normalized_quant.csv",
                             header = TRUE, row.names = 1, check.names = TRUE)

# 2.2 Transpose to Standard format (Samples as Rows, Features as Columns)
mRNA_mat <- t(as.matrix(D17_mRNA_data))
prot_mat <- t(as.matrix(D17_protein_data))

# 2.3 Align Samples
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

# Helper to clean and scale for visualization
clean_and_scale <- function(M) {
  vars <- colVars(M)
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

# Render
pheatmap(Z_vis_ordered,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = seq(-3, 3, length.out = 100),
         show_rownames = FALSE, show_colnames = TRUE,
         cluster_cols = FALSE, cluster_rows = FALSE, # Manual Order
         annotation_row = annotation_df,
         main = "Global Multi-Omics")


# ==============================================================================
# 4. CORE ALGORITHM: Kurtosis-Weighted MFA-FABIA
# ==============================================================================

message("--- Step 3: Running Hybrid MFA-FABIA ---")

#' 0 MFA-FABIA Hybrid model pipeline 
#' 1 MFA Stage: Used for preprocessing and weighting. 
#'   It equalizes the variance between omics layers (so Proteomics doesn't drown
#'   ... out mRNA) and assesses the global structure.
#' 2 Kurtosis Weighting: Instead of standard variance, we weight features by
#'   ... Kurtosis.
#'  High kurtosis indicates "heavy tails" (non-Gaussian expression), 
#'  which is the specific signal FABIA looks for.
#' 3 FABIA Stage: Performs sparse matrix factorization to find biclusters 
#'   ... (subsets of samples active on subsets of genes).
#' Run FABIA with feature weighting based on Kurtosis (Non-Normality)
#'
#' @param Xlist List of matrices (omics data).
#' @param gamma Tempering parameter for weights (0 = no weighting).
#' @param p FABIA parameter (number of biclusters).
#' @param alpha FABIA parameter (sparseness). Lower values = denser loadings.
#' @param block_equalize Boolean. Equalize block variance
#' @param scale_vars Boolean. Scale variables?
#' @param seed Random seed.
#' @return List containing model object, factor scores, loadings, weights, and kurtosis values.
mfa_weighted_fabia <- function(Xlist, 
                               gamma = 0.5, 
                               p = 4, 
                               alpha = 0.05, 
                               block_equalize = TRUE,
                               scale_vars = TRUE,
                               seed = 123) {
  
  set.seed(seed)
  
  # A. Robust Preprocessing (Remove Zero-Variance)
  Xlist_clean <- lapply(Xlist, function(M) {
    vars <- matrixStats::colVars(as.matrix(M))
    M[, vars > 1e-9, drop = FALSE]
  })
  
  # B. Block-size–normalized PCA 
  #Zlist <- lapply(Xlist_clean, function(X) {
  #  if(scale_vars) {
  #    s <- scale(X); s[is.na(s)] <- 0
  #  } else {
  #    s <- scale(X, scale = FALSE)
  #  }
  #  if(block_equalize) s <- s / sqrt(ncol(s))
  #  return(s)
  #})
  
  
  # B. MFA block equalizing 
  # 👉 (each block is weighted by the inverse of its first PCA eigenvalue:)
  # 👉 Ensures each block contributes equally to the first global dimension.
  # 👉 This equalizes structural variance, not just block size.
  
  Zlist <- lapply(Xlist_clean, function(X) {
    # scale variables
    s <- scale(X)
    s[is.na(s)] <- 0
    if(block_equalize) {
      pca_block <- prcomp(s, center = FALSE, scale. = FALSE)
      lambda1 <- pca_block$sdev[1]^2
      s <- s / sqrt(lambda1)
    }
    return(s)
  })
  
  
  # 👉 Rename features with prefix (Crucial for splitting later)
  for(n in names(Zlist)) colnames(Zlist[[n]]) <- paste0(n, "_", colnames(Zlist[[n]]))
  
  # Combine
  Z <- do.call(cbind, Zlist)
  
  # C. Feature Weighting (Kurtosis)
  # 👉 Features were weighted by excess kurtosis to emphasize sparse, heavy-tailed 
  # 👉 signals characteristic of bicluster structure, following the assumptions
  # 👉 of FABIA :: Laplace-like priors, Sparse loadings, Few active features per factor
  kurt_vals <- apply(Z, 2, e1071::kurtosis, type = 2, na.rm = TRUE)
  w_raw <- pmax(0, kurt_vals)
  mean_w <- mean(w_raw, na.rm = TRUE); if(mean_w==0) mean_w <- 1
  w <- (w_raw / mean_w)^gamma
  w[!is.finite(w)] <- 0
  
  # Apply Weights
  Z_w <- sweep(Z, 2, w, `*`)
  
  # D. Run FABIA
  fab_fit <- fabia::fabia(as.matrix(t(Z_w)), p = p, alpha = alpha, random = 1)
  
  # Format Outputs
  factor_scores <- t(fab_fit@Z) # Samples x Factors
  loadings_all  <- fab_fit@L    # Features x Factors
  
  # --- Split Loadings per Block ---
  # We use the prefixes we added earlier to reliably split the matrix back up
  loadings_per_block <- list()
  for(n in names(Zlist)) {
    # Find rows corresponding to this block (e.g., starts with "mRNA_")
    pattern <- paste0("^", n, "_")
    idx <- grep(pattern, rownames(loadings_all))
    
    if(length(idx) > 0) {
      loadings_per_block[[n]] <- loadings_all[idx, , drop=FALSE]
    }
  }
  
  return(list(
    model = fab_fit,
    factor_scores = factor_scores,
    loadings = loadings_all,         # Global matrix
    loadings_per_block = loadings_per_block, # Separated list
    weights = w,
    kurtosis = kurt_vals
  ))
}

# Run the Algorithm
res_hybrid <- mfa_weighted_fabia(X_list, 
                                 gamma = 0.5, 
                                 p = 4, 
                                 alpha = 0.05,
                                 block_equalize = TRUE,
                                 scale_vars = TRUE)

# Visualize Scores Heatmap
pheatmap(res_hybrid$factor_scores, cluster_cols = FALSE, display_numbers = TRUE,
         main = "Factor Scores (Sample Activity)")




# ==============================================================================
# 5. SCATTER PLOTS (for SAMPLE SCORES)
          # 👉 How strongly a sample expresses a latent factor)
# ==============================================================================

message("--- Step 4: Generating Scatter Plots ---")

# --- A. Load Metadata ---

D17_metadata <- read.csv(file = "../multi-omics radiation_data/D17_metadata.csv",
                         header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# --- B. Prepare Score Data ---
scores_df <- as.data.frame(res_hybrid$factor_scores)
colnames(scores_df) <- paste0("Factor", 1:ncol(scores_df))

# --- C. ALIGNMENT (Match "X13" to 13) ---

# 1. Clean the Sample IDs from the analysis (remove "X" prefix if present)
# "X13" -> "13"
scores_df$ID_clean <- sub("^X", "", rownames(scores_df))

# 2. Ensure IDs are the same data type (Integers) for merging
scores_df$ID_clean <- as.integer(scores_df$ID_clean)

# --- D. Merge Data ---
# Merge based on the cleaned ID from scores and 'Internal.ID' from metadata
scores_merged <- merge(
  x = scores_df, 
  y = D17_metadata, 
  by.x = "ID_clean",      # The cleaned ID (e.g., 13)
  by.y = "Internal.ID",   # The integer ID in metadata (e.g., 13)
  all.x = TRUE            # Keep all samples from analysis
)

# Check if merge worked
if(nrow(scores_merged) == 0) {
  stop("Merge failed! No common IDs")
} else {
  message(sprintf("Successfully merged metadata for %d samples.", nrow(scores_merged)))
}

# Add Sample Label (Original "X13" style name for plotting)
scores_merged$SampleLabel <- rownames(scores_df)[match(scores_merged$ID_clean, scores_df$ID_clean)]

# Convert grouping columns to factors
scores_merged$Description <- as.factor(scores_merged$Description)
scores_merged$Control     <- as.factor(scores_merged$Control)
scores_merged$Treatment   <- as.factor(scores_merged$Treatment)

# --- E. Generate Plots ---

my_theme <- theme_minimal() +
  theme(plot.title = element_text(face="bold", size=12),
        legend.position = "right")

# 1) Scatter Plot: NO GROUPING
p1 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 4, alpha = 0.7, color = "black") + 
  geom_text(vjust = -0.8, size = 3) +
  labs(title = "1. Scores (ungrouped)", x="Factor 1", y="Factor 2") +
  my_theme

# 2) Scatter Plot: Grouped by 'Description'
p2 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, color = Description, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "2. Scores by 'Description'", x="Factor 1", y="Factor 2") +
  my_theme

# 3) Scatter Plot: Grouped by 'Control'
p3 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, color = Control, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "3. Scores by 'Control'", x="Factor 1", y="Factor 2") +
  my_theme

# 4) Scatter Plot: Grouped by 'Treatment'
p4 <- ggplot(scores_merged, aes(x = Factor1, y = Factor2, color = Treatment, label = SampleLabel)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "4. Scores by 'Treatment'", x="Factor 1", y="Factor 2") +
  my_theme

# --- F. Print Plots ---
print(p1)
print(p2)
print(p3)
print(p4)

# 2x2 grid 
grid.arrange(p1, p2, p3, p4, nrow = 2)




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
# ==============================================================================
# Goal: Check specificity. Do Factor 1 genes also light up in Factor 2?

# 1. Extract Top 50 features from EACH factor and combine them
top_features_indices <- unique(unlist(lapply(1:ncol(res_hybrid$factor_scores), function(i) {
  order(abs(res_hybrid$loadings[,i]), decreasing = TRUE)[1:50]
})))

# 2. Subset the Loading Matrix
loading_subset <- res_hybrid$loadings[top_features_indices, ]

# 3. Create Annotation for Rows (Omic Type)
type_annot <- data.frame(Type = ifelse(grepl("^mRNA_", rownames(loading_subset)), "mRNA", "Protein"))
rownames(type_annot) <- rownames(loading_subset)

# 4. Render Heatmap
# We use a divergent palette: Blue (-) ... White (0) ... Red (+)
pheatmap(loading_subset,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = seq(-max(abs(loading_subset)), max(abs(loading_subset)), length.out = 100),
         annotation_row = type_annot,
         annotation_colors = list(Type = c(mRNA="#1f77b4", Protein="#ff7f0e")),
         show_rownames = FALSE, # Hide names if too crowded
         cluster_cols = FALSE,  # Keep Factors ordered 1..4
         main = "Heatmap of Top Loaded Features (Union of Top 50 per Factor)")


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



