# ==============================================================================
# SCRIPT: STANDARD MFA ON D17 DATASET
# ==============================================================================

# 1. SETUP
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(matrixStats)
library(pheatmap)

# 2. DATA LOADING & ALIGNMENT
message("Loading Data...")
# Load Metadata
D17_metadata <- read.csv(file = "../multi-omics radiation_data/D17_metadata.csv",
                         header = TRUE, row.names = 1, check.names = FALSE)

# Load Omics
mRNA_raw <- read.csv(file = "../multi-omics radiation_data/D17_mRNA_cor_normalized_tmm.csv",
                     header = TRUE, row.names = 1)
prot_raw <- read.csv(file = "../multi-omics radiation_data/D17_Protein_cor_normalized_quant.csv",
                     header = TRUE, row.names = 1)

# Transpose to Samples x Features
mRNA_mat <- t(as.matrix(mRNA_raw))
prot_mat <- t(as.matrix(prot_raw))

# Align Samples
common <- intersect(rownames(mRNA_mat), rownames(prot_mat))
mRNA_common <- mRNA_mat[common, ]
prot_common <- prot_mat[common, ]

# 3. PREPARE FOR MFA
# ------------------------------------------------------------------------------
# FactoMineR requires a single dataframe and a vector defining group sizes
df_mfa <- cbind(as.data.frame(mRNA_common), as.data.frame(prot_common))

# Define structure
group_sizes <- c(ncol(mRNA_common), ncol(prot_common))
group_names <- c("mRNA", "Protein")
group_types <- c("c", "c") # "c" for continuous data

message(sprintf("Running MFA on %d samples. Groups: mRNA (%d), Protein (%d)", 
                length(common), group_sizes[1], group_sizes[2]))

# 4. RUN MFA
# ------------------------------------------------------------------------------
# ncp=5: Keep first 5 dimensions
res_mfa <- MFA(df_mfa, 
               group = group_sizes, 
               type = group_types, 
               name.group = group_names, 
               ncp = 5, 
               graph = FALSE)

# 5. VISUALIZATION
# ------------------------------------------------------------------------------

# A. Scree Plot (Variance Explained)
fviz_screeplot(res_mfa, addlabels = TRUE) +
  ggtitle("MFA: Scree Plot (Global Variance)")

# B. Group Plot (Inertia)
# Checks if mRNA dominates Protein or vice versa
# What it shows: The "center of gravity" for the mRNA vs. Protein datasets.
# Use: Diagnostics. If one point is far from the origin and the other is close, 
# the analysis is dominated by one omic layer (unbalanced).
fviz_mfa_var(res_mfa, "group") +
  ggtitle("MFA: Contribution of Omics Layers")

# C. Sample Plot (PCA-like Map)
# Align metadata to common samples for coloring
meta_aligned <- D17_metadata[match(gsub("^X", "", common), D17_metadata$Internal.ID), ]

# Individual Factor Map (Sample Space):
# What it shows: Samples projected into the latent space (e.g., Dim 1 vs. Dim 2).
# Use: Identifying clusters (e.g., Treatment vs. Control) and outliers
fviz_mfa_ind(res_mfa, 
             habillage = as.factor(meta_aligned$Treatment), # <--- Use habillage for grouping
             geom = "point",
             pointsize = 3,
             palette = c("grey", "red"),     # Colors map to the groups in habillage
             addEllipses = TRUE, 
             ellipse.type = "confidence",
             repel = TRUE) +
  ggtitle("MFA: Sample Space (Colored by Treatment)")

fviz_mfa_ind(res_mfa, 
             habillage = as.factor(meta_aligned$Control),
             geom = "point",
             pointsize = 3,
             palette = c("grey", "red"),
             addEllipses = TRUE, 
             ellipse.type = "confidence",
             repel = TRUE) +
  ggtitle("MFA: Sample Space (Colored by Control)")

fviz_mfa_ind(res_mfa, 
             habillage = as.factor(meta_aligned$Description),
             geom = "point",
             pointsize = 3,
             addEllipses = TRUE, 
             ellipse.type = "confidence",
             repel = TRUE) +
  ggtitle("MFA: Sample Space (Colored by Description)")



# Top 20 contributing features for Dimension 1
fviz_contrib(res_mfa, choice = "quanti.var", axes = 1, top = 20)

# Extract MFA coordinates (analogous to Factor Scores)
mfa_scores <- res_mfa$ind$coord[, 1:4] # First 4 dims

# Plot using same annotation as FABIA
pheatmap(mfa_scores, 
         annotation_row = annot_df,      # Use same annotation from FABIA script
         annotation_colors = ann_colors, # Use same colors
         main = "MFA: Dimension Scores (Global)")

# Variable Correlation Circle (Global Loadings):
# What it shows: How strongly features correlate with the global dimensions.
# Use: Identifying which genes/proteins drive the separation of samples.
fviz_mfa_var(res_mfa, "quanti.var")
