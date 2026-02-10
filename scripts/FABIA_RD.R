# ==============================================================================
# SCRIPT: FABIA ON D17 DATASET (CLEANED)
# ==============================================================================

# 1. SETUP & LIBRARIES
# ==============================================================================
library(fabia)
library(pheatmap)
library(matrixStats)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

set.seed(123)

# 2. DATA LOADING & ALIGNMENT
# ==============================================================================
message("--- Step 1: Loading Data ---")

# Load Metadata
D17_metadata <- read.csv(file = "../multi-omics radiation_data/D17_metadata.csv",
                         header = TRUE, row.names = 1, check.names = FALSE)

# Load Omics Data
mRNA_raw <- read.csv(file = "../multi-omics radiation_data/D17_mRNA_cor_normalized_tmm.csv",
                     header = TRUE, row.names = 1)
prot_raw <- read.csv(file = "../multi-omics radiation_data/D17_Protein_cor_normalized_quant.csv",
                     header = TRUE, row.names = 1)

# Transpose to Samples x Features
mRNA_mat <- t(as.matrix(mRNA_raw))
prot_mat <- t(as.matrix(prot_raw))

# Align Samples
common <- intersect(rownames(mRNA_mat), rownames(prot_mat))
if(length(common) == 0) stop("No matching samples found.")

mRNA_common <- mRNA_mat[common, ]
prot_common <- prot_mat[common, ]

message(sprintf("Aligned %d samples. mRNA features: %d, Protein features: %d", 
                length(common), ncol(mRNA_common), ncol(prot_common)))

# 3. PREPROCESSING
# ==============================================================================
message("--- Step 2: Preprocessing ---")

# Combine matrices
X_combined <- cbind(mRNA_common, prot_common)

# Clean: Remove zero-variance columns
vars <- colVars(X_combined)
X_clean <- X_combined[, vars > 1e-9]

# Scale: Center and Scale (Critical for FABIA)
X_scaled <- scale(X_clean)
X_scaled[is.na(X_scaled)] <- 0

# Transpose for FABIA (Needs Features x Samples)
X_input <- t(X_scaled)

message(sprintf("Final Input Matrix: %d Features x %d Samples", nrow(X_input), ncol(X_input)))

# 4. RUN FABIA
# ==============================================================================
message("--- Step 3: Running FABIA ---")

# p=4: Extract 4 biclusters
# alpha=0.05: Sparseness level
res_fabia <- fabia(X_input, p = 2, alpha = 0.05, random = 1)

# 5. POST-PROCESSING (Extract & Organsize Results)
# ==============================================================================
message("--- Step 4: Extracting Results ---")

# A. Extract Factor Scores
scores <- t(res_fabia@Z) # Samples x Factors
colnames(scores) <- paste0("Factor", 1:4)

# B. Extract Global Loadings
global_loadings <- res_fabia@L

# C. Split Loadings back into mRNA and Protein
n_mrna_features <- ncol(mRNA_common)
loadings_mRNA <- global_loadings[1:n_mrna_features, , drop = FALSE]
loadings_Prot <- global_loadings[(n_mrna_features + 1):nrow(global_loadings), , drop = FALSE]

# D. Prepare Merged Dataframes for Plotting (Done once here)
scores_df <- as.data.frame(scores)
scores_df$SampleID <- rownames(scores_df)
scores_df$ID_clean <- as.integer(gsub("^X", "", scores_df$SampleID))

# Merge with Metadata
scores_merged <- merge(x = scores_df, y = D17_metadata, 
                       by.x = "ID_clean", by.y = "Internal.ID", all.x = TRUE)

# Convert metadata to factors
scores_merged$Description <- as.factor(scores_merged$Description)
scores_merged$Control     <- as.factor(scores_merged$Control)
scores_merged$Treatment   <- as.factor(scores_merged$Treatment)

# Create Long Format for Scan Plots
long_scores <- melt(scores_merged,
                    id.vars = c("SampleID", "Description", "Control", "Treatment"),
                    measure.vars = colnames(scores), 
                    variable.name = "Factor", value.name = "Score")

# 6. VISUALIZATION 1: FACTOR SCORES HEATMAP
# ==============================================================================
message("Generating Scores Heatmap...")

# Prepare Annotation
meta_idx <- match(scores_df$ID_clean, D17_metadata$Internal.ID)
annot_df <- data.frame(
  Description = as.factor(D17_metadata$Description[meta_idx]),
  Treatment   = as.factor(D17_metadata$Treatment[meta_idx]),
  Control     = as.factor(D17_metadata$Control[meta_idx])
)
rownames(annot_df) <- rownames(scores)

# Define Colors
ann_colors <- list(
  Treatment = c("TRUE" = "firebrick", "FALSE" = "grey90"),
  Control   = c("TRUE" = "forestgreen", "FALSE" = "grey90"),
  Description = c("CorP17_2Gy" = "#E41A1C", "CorP17_Ctrl" = "#4DAF4A",
                  "CorP17FA_2Gy" = "#377EB8", "CorP17FA_Ctrl" = "#984EA3")
)

pheatmap(scores, annotation_row = annot_df, annotation_colors = ann_colors,
         cluster_cols = FALSE, cluster_rows = TRUE, display_numbers = TRUE,
         main = "FABIA: Factor Scores")

# 7. VISUALIZATION 2: LOADING PLOTS (Lollipop & Heatmap)
# ==============================================================================
message("Generating Loading Plots...")

# Function: Lollipop Plot
plot_top_loadings <- function(loading_mat, omic_name, bar_color, top_n = 15) {
  df <- as.data.frame(loading_mat)
  colnames(df) <- paste0("Factor", 1:ncol(df))
  df$FeatureID <- rownames(df)
  plot_list <- list()
  
  for(i in 1:ncol(loading_mat)) {
    factor_col <- paste0("Factor", i)
    top_data <- df %>% arrange(desc(abs(get(factor_col)))) %>% head(top_n)
    
    p <- ggplot(top_data, aes(x = reorder(FeatureID, abs(get(factor_col))), y = get(factor_col))) +
      geom_segment(aes(x = FeatureID, xend = FeatureID, y = 0, yend = get(factor_col)), 
                   color = bar_color, alpha = 0.6) +
      geom_point(size = 3, color = bar_color) +
      coord_flip() +
      labs(title = paste0(omic_name, ": Factor ", i), x = NULL, y = "Loading Weight") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10, face = "bold"), axis.text.y = element_text(size = 8))
    plot_list[[i]] <- p
  }
  return(plot_list)
}

# Function: Loading Heatmap
plot_loading_heatmap <- function(loading_matrix, title_text, top_n = 40) {
  top_indices <- unique(unlist(lapply(1:ncol(loading_matrix), function(i) {
    order(abs(loading_matrix[, i]), decreasing = TRUE)[1:top_n]
  })))
  subset_mat <- loading_matrix[top_indices, , drop = FALSE]
  limit <- max(abs(subset_mat))
  
  pheatmap(subset_mat,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           breaks = seq(-limit, limit, length.out = 100),
           cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE, fontsize_row = 5,
           main = paste0(title_text, " (Top ", top_n, "/Factor)"))
}

# Execute Loading Plots
grid.arrange(grobs = plot_top_loadings(loadings_mRNA, "mRNA", "#1f77b4"), ncol = 2, top = "Top mRNA Drivers")
grid.arrange(grobs = plot_top_loadings(loadings_Prot, "Protein", "#ff7f0e"), ncol = 2, top = "Top Protein Drivers")

plot_loading_heatmap(loadings_mRNA, "FABIA: mRNA Top Loadings")
plot_loading_heatmap(loadings_Prot, "FABIA: Protein Top Loadings")

# 8. VISUALIZATION 3: FACTOR SCAN PLOTS
# ==============================================================================
message("Generating Factor Scan Plots...")

plot_factor_scan <- function(data, group_col, title_text) {
  # Sort data
  unique_meta <- unique(data[, c("SampleID", group_col)])
  unique_meta <- unique_meta[order(unique_meta[[group_col]], unique_meta$SampleID), ]
  data$OrderedSample <- factor(data$SampleID, levels = unique_meta$SampleID)
  
  ggplot(data, aes(x = OrderedSample, y = Score, color = .data[[group_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_segment(aes(x = OrderedSample, xend = OrderedSample, y = 0, yend = Score), alpha=0.5) +
    geom_point(size = 3, alpha = 0.9) +
    facet_wrap(~Factor, ncol = 2, scales = "free_y") +
    labs(title = title_text, subtitle = paste("Ordered by", group_col), x = "Samples", y = "Factor Score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8),
          strip.text = element_text(face = "bold", size = 10), legend.position = "right")
}

print(plot_factor_scan(long_scores, "Control", "FABIA: Factor Scan (Control)"))
print(plot_factor_scan(long_scores, "Treatment", "FABIA: Factor Scan (Treatment)"))
print(plot_factor_scan(long_scores, "Description", "FABIA: Factor Scan (Description)"))

# 9. VISUALIZATION 4: SCATTERPLOTS
# ==============================================================================
message("Generating Scatterplots...")

create_scatter_panel <- function(data, x_fac, y_fac, title_prefix) {
  my_theme <- theme_minimal() + theme(plot.title = element_text(face="bold", size=10), legend.position = "right")
  
  p1 <- ggplot(data, aes(x = .data[[x_fac]], y = .data[[y_fac]], label = SampleID)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_point(size = 2, alpha = 0.7, color = "black") + geom_text(vjust = -0.8, size = 3) +
    labs(title = paste0("1. ", title_prefix, " (Ungrouped)"), x=x_fac, y=y_fac) + my_theme
  
  p2 <- ggplot(data, aes(x = .data[[x_fac]], y = .data[[y_fac]], color = Description)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = paste0("2. ", title_prefix, " by Description"), x=x_fac, y=y_fac) + my_theme
  
  p3 <- ggplot(data, aes(x = .data[[x_fac]], y = .data[[y_fac]], color = Control)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("FALSE"="grey70", "TRUE"="forestgreen")) +
    labs(title = paste0("3. ", title_prefix, " by Control"), x=x_fac, y=y_fac) + my_theme
  
  p4 <- ggplot(data, aes(x = .data[[x_fac]], y = .data[[y_fac]], color = Treatment)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("FALSE"="grey70", "TRUE"="firebrick")) +
    labs(title = paste0("4. ", title_prefix, " by Treatment"), x=x_fac, y=y_fac) + my_theme
  
  grid.arrange(p1, p2, p3, p4, nrow = 2, top = paste("FABIA:", x_fac, "vs", y_fac))
}

create_scatter_panel(scores_merged, "Factor1", "Factor2", "Scores")
create_scatter_panel(scores_merged, "Factor3", "Factor4", "Scores")

