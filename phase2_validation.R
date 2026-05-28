# ============================================================
# LIVER DECONVOLUTION PROJECT — PHASE 2
# In silico validation using 100 synthetic liver tissues
# ============================================================

# ---- Step 2.1: Load libraries and Phase 1 objects ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)

setwd("/Users/mayankgandhi/Downloads/Projects/LiverDecon")

liver_ref     <- readRDS("data/processed/liver_ref.rds")
liver_sampled <- readRDS("data/processed/liver_sampled.rds")

cat("liver_ref cells:", ncol(liver_ref), "\n")
cat("liver_sampled cells:", ncol(liver_sampled), "\n")

# ---- Step 2.2: Identify held-out cells ----

# Cells used in reference
sampled_barcodes <- colnames(liver_sampled)

# Held-out cells = in liver_ref but NOT in liver_sampled
heldout_cells <- liver_ref[, !colnames(liver_ref) %in% sampled_barcodes]

cat("Held-out cells available:", ncol(heldout_cells), "\n")
table(heldout_cells@meta.data$Lineage)

# ---- Step 2.3: Generate 100 synthetic tissues ----

set.seed(42)

n_synthetic  <- 100
n_cells_each <- 1500
lineages     <- c("Hepatocyte", "Kupffer", "LSEC", "Lymphoid")

# Storage matrices
synthetic_counts <- matrix(0,
                           nrow = nrow(heldout_cells),
                           ncol = n_synthetic)
rownames(synthetic_counts) <- rownames(heldout_cells)
colnames(synthetic_counts) <- paste0("Synthetic_", 1:n_synthetic)

known_proportions <- matrix(0,
                            nrow = n_synthetic,
                            ncol = length(lineages))
rownames(known_proportions) <- paste0("Synthetic_", 1:n_synthetic)
colnames(known_proportions) <- lineages

# Generate each synthetic tissue
for (i in 1:n_synthetic) {
  
  props <- abs(rnorm(length(lineages), mean = 0.25, sd = 0.25))
  props <- props / sum(props)
  
  n_per_lineage <- round(props * n_cells_each)
  
  for (j in seq_along(lineages)) {
    available <- sum(heldout_cells@meta.data$Lineage == lineages[j])
    n_per_lineage[j] <- min(n_per_lineage[j], available)
  }
  
  true_props <- n_per_lineage / sum(n_per_lineage)
  known_proportions[i, ] <- true_props
  
  selected_cells <- c()
  for (j in seq_along(lineages)) {
    lineage_cells <- colnames(heldout_cells)[
      heldout_cells@meta.data$Lineage == lineages[j]
    ]
    selected <- sample(lineage_cells, n_per_lineage[j], replace = FALSE)
    selected_cells <- c(selected_cells, selected)
  }
  
  cell_matrix <- GetAssayData(heldout_cells[, selected_cells], layer = "counts")
  synthetic_counts[, i] <- rowSums(cell_matrix)
}

cat("Synthetic tissues generated:", ncol(synthetic_counts), "\n")
cat("Known proportions summary:\n")
print(round(colMeans(known_proportions), 3))

# Save synthetic tissues
saveRDS(list(counts            = synthetic_counts,
             known_proportions = known_proportions),
        file = "data/processed/synthetic_tissues.rds")

cat("Saved synthetic tissues to data/processed/\n")

# ---- Step 2.4: Export synthetic tissues for CIBERSORTx ----

synthetic_df <- as.data.frame(synthetic_counts)
synthetic_df <- tibble::rownames_to_column(synthetic_df, var = "GeneSymbol")

# Replace dashes with underscores for CIBERSORTx compatibility
synthetic_df$GeneSymbol <- gsub("-", "_", synthetic_df$GeneSymbol)

write.table(synthetic_df,
            file      = "data/processed/synthetic_mixture_clean.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("CIBERSORTx mixture file written\n")
cat("Genes:", nrow(synthetic_counts), "\n")
cat("Samples:", ncol(synthetic_counts), "\n")