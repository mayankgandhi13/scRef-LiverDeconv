# ============================================================
# LIVER DECONVOLUTION PROJECT — PHASE 1
# Building single-cell reference from MacParland 2018 (GSE115469)
# ============================================================

# ---- Step 1: Load libraries ----
library(GEOquery)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

# ---- Step 2: Set working directory and folder structure ----
setwd("/Users/mayankgandhi/Downloads/Projects/LiverDecon")

dir.create("data/raw",        showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed",  showWarnings = FALSE, recursive = TRUE)
dir.create("data/reference",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Step 3: Load metadata ----
metadata <- read.table("data/raw/GSE115469_CellClusterType.txt",
                       header = FALSE,
                       sep = "\t",
                       skip = 1,
                       col.names = c("CellName", "Sample",
                                     "CellBarcode", "Cluster", "CellType"))

# Quick checks
dim(metadata)
unique(metadata$CellType)
table(metadata$CellType)

# ---- Step 4: Load count matrix ----
counts <- read.csv("data/raw/GSE115469_Data.csv",
                   row.names = 1)

# Confirm dimensions: 20007 genes x 8444 cells
dim(counts)

# Confirm cell names match between counts and metadata
all(colnames(counts) == metadata$CellName)

# ---- Step 5: Quality control ----

# Create Seurat object
liver_sc <- CreateSeuratObject(
  counts  = counts,
  project = "LiverAtlas"
)

# Calculate mitochondrial percentage per cell
liver_sc[["percent.mt"]] <- PercentageFeatureSet(liver_sc, pattern = "^MT-")

# Look at QC metrics
head(liver_sc@meta.data)

# Visualise QC distributions across all 5 donors
VlnPlot(liver_sc,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3)

# Save the plot
ggsave("outputs/figures/QC_before_filtering.png", width = 12, height = 5)

# Apply QC filters (MacParland's own cutoffs)
liver_sc <- subset(liver_sc,
                   subset = nCount_RNA   >= 1500 &
                     percent.mt   <= 50   &
                     nFeature_RNA >= 200)

# How many cells remain?
dim(liver_sc)

# ---- Step 6: Add cell type labels and merge into lineages ----

# Match metadata to the cells remaining after QC
liver_sc$CellType <- metadata$CellType[match(colnames(liver_sc), 
                                             metadata$CellName)]

# Check it worked
head(liver_sc$CellType)

# Check updated counts per cell type after QC
table(liver_sc$CellType)

# Define lineage mapping
lineage_map <- c(
  "Hepatocyte_1"                = "Hepatocyte",
  "Hepatocyte_2"                = "Hepatocyte",
  "Hepatocyte_3"                = "Hepatocyte",
  "Hepatocyte_4"                = "Hepatocyte",
  "Hepatocyte_5"                = "Hepatocyte",
  "Hepatocyte_6"                = "Hepatocyte",
  "Central_venous_LSECs"        = "LSEC",
  "Periportal_LSECs"            = "LSEC",
  "Inflammatory_Macrophage"     = "Kupffer",
  "Non-inflammatory_Macrophage" = "Kupffer",
  "alpha-beta_T_Cells"          = "Lymphoid",
  "gamma-delta_T_Cells_1"       = "Lymphoid",
  "gamma-delta_T_Cells_2"       = "Lymphoid",
  "NK-like_Cells"               = "Lymphoid",
  "Plasma_Cells"                = "Lymphoid"
)

# Bypass Seurat's metadata assignment using direct slot access
cell_types <- as.character(liver_sc$CellType)
lineage_vector <- lineage_map[cell_types]
names(lineage_vector) <- colnames(liver_sc)
liver_sc@meta.data$Lineage <- lineage_vector

# Check results
table(liver_sc@meta.data$Lineage, useNA = "ifany")

# ---- Step 7: Sample 200 cells per lineage ----

# Remove excluded cell types (NA lineage)
liver_ref <- liver_sc[, !is.na(liver_sc@meta.data$Lineage)]

# Confirm only 4 lineages remain
table(liver_ref@meta.data$Lineage)

# Set seed for reproducibility
set.seed(42)

# Sample 200 cells per lineage
sampled_cells <- liver_ref@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  group_by(Lineage) %>%
  slice_sample(n = 200) %>%
  pull(cell_barcode)

# Subset Seurat object to sampled cells only
liver_sampled <- liver_ref[, sampled_cells]

# Confirm: should be exactly 800 cells (200 x 4 lineages)
table(liver_sampled@meta.data$Lineage)

# ---- Step 8: Export CIBERSORTx reference ----

# Extract raw count matrix
ref_matrix <- GetAssayData(liver_sampled, layer = "counts")

# Build column names: Lineage_CellNumber (e.g. Hepatocyte_1, Kupffer_1)
lineage_labels <- liver_sampled@meta.data$Lineage
cell_numbers   <- ave(seq_along(lineage_labels),
                      lineage_labels,
                      FUN = seq_along)
colnames(ref_matrix) <- paste0(lineage_labels, "_", cell_numbers)

# Convert to data frame and add GeneSymbol column
ref_df <- as.data.frame(as.matrix(ref_matrix))
ref_df <- tibble::rownames_to_column(ref_df, var = "GeneSymbol")

# Write to file
write.table(ref_df,
            file      = "data/reference/liver_CIBERSORTx_reference.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("File written:", nrow(ref_df), "genes x", ncol(ref_df)-1, "cells\n")

# ---- Step 9: Export BayesPrism reference ----

# Identify genes to exclude
ribo_genes <- rownames(liver_sampled)[grepl("^RP[SL]", rownames(liver_sampled))]
mito_genes <- rownames(liver_sampled)[grepl("^MT-",    rownames(liver_sampled))]
exclude     <- unique(c(ribo_genes, mito_genes))

cat("Genes excluded:", length(exclude), "\n")

# Filter matrix and transpose (BayesPrism needs cells x genes)
bp_matrix <- ref_matrix[!rownames(ref_matrix) %in% exclude, ]
bp_sc     <- t(as.matrix(bp_matrix))

# Cell type labels
bp_cell_type  <- liver_sampled@meta.data$Lineage
bp_cell_state <- liver_sampled@meta.data$Lineage

# Save as RDS
saveRDS(
  list(sc_matrix  = bp_sc,
       cell_type  = bp_cell_type,
       cell_state = bp_cell_state),
  file = "data/reference/liver_BayesPrism_reference.rds"
)

cat("BayesPrism reference saved:", nrow(bp_sc), "cells x", ncol(bp_sc), "genes\n")


