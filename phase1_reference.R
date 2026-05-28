# ============================================================
# LIVER DECONVOLUTION PROJECT — PHASE 1
# Building single-cell reference from MacParland 2018 (GSE115469)
# ============================================================

# ---- Step 1.1: Load libraries ----
library(GEOquery)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

# ---- Step 1.2: Set working directory and folder structure ----
setwd("/Users/mayankgandhi/Downloads/Projects/LiverDecon")

dir.create("data/raw",        showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed",  showWarnings = FALSE, recursive = TRUE)
dir.create("data/reference",  showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Step 1.3: Load metadata ----
metadata <- read.table("data/raw/GSE115469_CellClusterType.txt",
                       header = FALSE,
                       sep = "\t",
                       skip = 1,
                       col.names = c("CellName", "Sample",
                                     "CellBarcode", "Cluster", "CellType"))

dim(metadata)
unique(metadata$CellType)
table(metadata$CellType)

# ---- Step 1.4: Load count matrix ----
counts <- read.csv("data/raw/GSE115469_Data.csv",
                   row.names = 1)

dim(counts)
all(colnames(counts) == metadata$CellName)

# ---- Step 1.5: Quality control ----
liver_sc <- CreateSeuratObject(
  counts  = counts,
  project = "LiverAtlas"
)

liver_sc[["percent.mt"]] <- PercentageFeatureSet(liver_sc, pattern = "^MT-")

head(liver_sc@meta.data)

VlnPlot(liver_sc,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        ncol = 3)

ggsave("outputs/figures/QC_before_filtering.png", width = 12, height = 5)

liver_sc <- subset(liver_sc,
                   subset = nCount_RNA   >= 1500 &
                     percent.mt   <= 50   &
                     nFeature_RNA >= 200)

dim(liver_sc)

# ---- Step 1.6: Add cell type labels and merge into lineages ----
liver_sc$CellType <- metadata$CellType[match(colnames(liver_sc),
                                             metadata$CellName)]

head(liver_sc$CellType)
table(liver_sc$CellType)

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

cell_types <- as.character(liver_sc$CellType)
lineage_vector <- lineage_map[cell_types]
names(lineage_vector) <- colnames(liver_sc)
liver_sc@meta.data$Lineage <- lineage_vector

table(liver_sc@meta.data$Lineage, useNA = "ifany")

# ---- Step 1.7: Sample 200 cells per lineage ----
liver_ref <- liver_sc[, !is.na(liver_sc@meta.data$Lineage)]

table(liver_ref@meta.data$Lineage)

set.seed(42)

sampled_cells <- liver_ref@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  group_by(Lineage) %>%
  slice_sample(n = 200) %>%
  pull(cell_barcode)

liver_sampled <- liver_ref[, sampled_cells]

table(liver_sampled@meta.data$Lineage)

# ---- Step 1.8: Export CIBERSORTx reference ----
ref_matrix <- GetAssayData(liver_sampled, layer = "counts")

# Round to integers — CIBERSORTx requires whole number counts
ref_matrix <- round(ref_matrix)

lineage_labels <- liver_sampled@meta.data$Lineage
cell_numbers   <- ave(seq_along(lineage_labels),
                      lineage_labels,
                      FUN = seq_along)
colnames(ref_matrix) <- paste0(lineage_labels, "_", cell_numbers)

ref_df <- as.data.frame(as.matrix(ref_matrix))
ref_df <- tibble::rownames_to_column(ref_df, var = "GeneSymbol")

# Replace dashes with underscores for CIBERSORTx compatibility
ref_df$GeneSymbol <- gsub("-", "_", ref_df$GeneSymbol)

# Keep only genes expressed in at least 10 cells
gene_counts <- rowSums(ref_df[,-1] > 0)
ref_df      <- ref_df[gene_counts >= 10, ]

write.table(ref_df,
            file      = "data/reference/liver_MacParland_final.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("Final reference written:", nrow(ref_df), "genes x", ncol(ref_df)-1, "cells\n")

# ---- Step 1.9: Export BayesPrism reference ----
ribo_genes <- rownames(liver_sampled)[grepl("^RP[SL]", rownames(liver_sampled))]
mito_genes <- rownames(liver_sampled)[grepl("^MT-",    rownames(liver_sampled))]
exclude     <- unique(c(ribo_genes, mito_genes))

cat("Genes excluded:", length(exclude), "\n")

bp_matrix <- ref_matrix[!rownames(ref_matrix) %in% exclude, ]
bp_sc     <- t(as.matrix(bp_matrix))

bp_cell_type  <- liver_sampled@meta.data$Lineage
bp_cell_state <- liver_sampled@meta.data$Lineage

saveRDS(
  list(sc_matrix  = bp_sc,
       cell_type  = bp_cell_type,
       cell_state = bp_cell_state),
  file = "data/reference/liver_BayesPrism_reference.rds"
)

cat("BayesPrism reference saved:", nrow(bp_sc), "cells x", ncol(bp_sc), "genes\n")

# ---- Save R objects for Phase 2 ----
saveRDS(liver_ref,     file = "data/processed/liver_ref.rds")
saveRDS(liver_sampled, file = "data/processed/liver_sampled.rds")

cat("Objects saved to data/processed/\n")

# ---- Step 1.10: Export filtered CIBERSORTx reference ----

# Keep only genes expressed in at least 10 cells (reduces file size)
ref_check   <- as.data.frame(as.matrix(GetAssayData(liver_sampled, layer = "counts")))
gene_counts <- rowSums(ref_check > 0)
cat("Genes expressed in at least 10 cells:", sum(gene_counts >= 10), "\n")

ref_filtered    <- ref_check[gene_counts >= 10, ]
ref_filtered_df <- tibble::rownames_to_column(ref_filtered, var = "GeneSymbol")

# Apply dash to underscore fix
ref_filtered_df$GeneSymbol <- gsub("-", "_", ref_filtered_df$GeneSymbol)

# Rebuild column names
colnames(ref_filtered_df)[-1] <- paste0(lineage_labels, "_", cell_numbers)

write.table(ref_filtered_df,
            file      = "data/reference/liver_MacParland_filtered.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("Filtered reference written:", nrow(ref_filtered), "genes x", ncol(ref_filtered), "cells\n")
