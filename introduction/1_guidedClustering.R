#!/usr/bin/env Rscript
#-------------------------------------------------
# Title: Guided Clustering Tutorial
# Author: Amanda Zacharias
# Date: 2023-12-08
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
# module load StdEnv/2023 r/4.3.1
#
# Options -----------------------------------------

# Packages -----------------------------------------
library(dplyr) # 1.1.4; data wrangling
library(Seurat) # 5.0.1; scRNA data analysis
library(patchwork) # 1.1.3; combine ggplots
library(ggplot2) # 3.4.4

# Pathways -----------------------------------------
baseDir <- file.path("introduction", "1_guidedClustering")
# Input ===========
dataLink <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

# Output ===========
scDataDir <- file.path(dirname(baseDir), "data")
plotsDir <- file.path(baseDir, "plots")
rDataDir <- file.path(baseDir, "rData")
system(paste("mkdir", scDataDir, plotsDir, rDataDir))

# Load data -----------------------------------------
# Download data =========
# system(paste('wget -P', scDataDir, dataLink))
# system(paste('tar -xvzf', file.path(scDataDir, basename(dataLink)), "-C", scDataDir))

# Read in data as matrix ========
pbmc.data <- Read10X(file.path(scDataDir, "filtered_gene_bc_matrices", "hg19"))

# Create Seurat object ======
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# > pbmc
# An object of class Seurat
# 13714 features across 2700 samples within 1 assay
# Active assay: RNA (13714 features, 0 variable features)
# 1 layer present: counts

# Count matrix data structure -----------------------------------------
pbmc.data[1:3, 1:10] # 3 genes, 10 cells
# 3 x 10 sparse Matrix of class "dgCMatrix"
# [[ suppressing 10 column names ‘AAACATACAACCAC-1’, ‘AAACATTGAGCTAC-1’, ‘AAACATTGATCAGC-1’ ... ]]
# MIR1302-10 . . . . . . . . . .
# FAM138A    . . . . . . . . . .
# OR4F5      . . . . . . . . . .

# Seurat uses a sparse matrix to represent data and save memory.
# "." indicates that X feature is not present in Y cell

# Quality Control -----------------------------------------
# Basic metrics:
# 1. Number of unique genes in each cell
#   - catch empty droplets and doublets
# 2. Total number of molecules in each cell
#   - high association with metric 1.
# 3. Percentage of reads that map to mitochondrial genome
#   - dying cells have high mt
#   - use the MT- genes
# See automatically calculated metrics =======
pbmc@meta.data[1:5, ]
#                  orig.ident nCount_RNA nFeature_RNA
# AAACATACAACCAC-1     pbmc3k       2419          779
# AAACATTGAGCTAC-1     pbmc3k       4903         1352
# AAACATTGATCAGC-1     pbmc3k       3147         1129
# AAACCGTGCTTCCG-1     pbmc3k       2639          960
# AAACCGTGTATGCG-1     pbmc3k        980          521
# nCount_RNA = Total number of molecules in each cell
# nFeature_RNA = Number of unique features(genes) in each cell

# Calculate MT percentage =======
# [["newColumn]] Adds a column to the metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Violin plots =======
qcViolin <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
ggsave(plot = qcViolin, filename = "rawViolins.pdf", path = plotsDir,
       width = 370, height = 185, units = "mm")

nCountVsMT <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
nCountVsNFeature <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(plot = nCountVsMT + nCountVsNFeature,
       filename = "rawFeatureScatter.pdf", path = plotsDir,
       width = 370, height = 185, units = "mm")

rm(qcViolin)
rm(nCountVsMT)
rm(nCountVsNFeature)

# QC Filtering =======
# Keep cells with a unique number of features >200 and <2500.
# Keep cells if their percentage of MT is < 5%.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing -----------------------------------------
# Default parameters are being used.
# Feature counts for each cell are divided by the total counts for that cell
# and multiplied by the 10,000 scale factor.
# Then, they are natural-log transformed using `log1p()`, which is "log(1+x)".
# Stored in pbmc[["RNA"]]$data
pbmc <- NormalizeData(pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)
# Other methods exist!

# Feature selection -----------------------------------------
# Identify features with high variance across cells because they are
# likely to be biologically interesting.
# `FindVariableFeatures` models the mean-variance relationship inherent in scRNA.
# Using the vst method. Top 2,000 features are returnned.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Top 10 variable
top10 <- head(VariableFeatures(pbmc), 10)

# Plotting
variableNoLabels <- VariableFeaturePlot(pbmc)
variableLabels <- LabelPoints(plot = variableNoLabels, points = top10, repel = TRUE)
ggsave(plot = variableNoLabels + variableLabels,
       filename = "variableFeatures.pdf", path = plotsDir,
       width = 370, height = 185, units = "mm")
rm(variableNoLabels, variableLabels)

# Scaling for Dimensionality Reduction -----------------------------------------
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# Prevent highly expressed genes from dominating the downstream clustering.
# Only variable features are scaled by default.
# Results stored in pbmc[["RNA"]]$scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Regressing out unwanted variation (technical effects): use vars.to.regress parameter.
# Recommended to use the `SCTransform` function instead!

# Linear dimensional reduction (PCA) -----------------------------------------
# Calculate =======
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Inspect =======
# Access with pbmc[["pca"]]

# Seurat lists genes w/ positive and negative correlations across cells for 1st PCs
print(pbmc[["pca"]], dims = 1:2, nfeatures = 5)
# PC_ 1
# Positive:  CST3, TYROBP, LST1, AIF1, FTL
# Negative:  MALAT1, LTB, IL32, IL7R, CD2
# PC_ 2
# Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1
# Negative:  NKG7, PRF1, CST7, GZMB, GZMA

# Plots =======
vizDim <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
ggsave(plot = vizDim,
       filename = "vizDim.pdf", path = plotsDir,
       width = 185, height = 185, units = "mm")

dimPlot <- DimPlot(pbmc, reduction = "pca") + NoLegend()
ggsave(plot = dimPlot,
       filename = "dimPlot.pdf", path = plotsDir,
       width = 185, height = 185, units = "mm")

rm(vizDim, dimPlot)

pdf(file.path(plotsDir, "dimHeatmap.pdf"), width = 4, height = 4)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE, combine = TRUE)
dev.off()

pdf(file.path(plotsDir, "dimHeatmap1to10.pdf"), width = 7, height = 7)
DimHeatmap(pbmc, dims = 1:10, cells = 500, balanced = TRUE, combine = TRUE)
dev.off()

# Optimize number of PCs to include =======
# Elbow plot
elbow <- ElbowPlot(pbmc)
ggsave(plot = elbow,
       filename = "elbow.pdf", path = plotsDir,
       width = 90, height = 90, units = "mm")
# Recommendations:
# 1. Supervised
# 2. Elbow plot
#   - Lean on the higher side
#   - Try multiple different thresholds. results probs won't change much
# 3. Prior knowledge.
#   Ex. genes in PCs 12-13 seem important for rare DC and NK cells

# Cluster cells -----------------------------------------
# Seurat takes a graph-based approach. In brief:
#      1) group cells together based on similarity (KNN),
#      2) partition into groups
# Note for partitioning step:  0.4-1.2 'granularity' is usually
# good for data w/ ~3k cells. Resolution often increases w/ # of cells.
pbmc <- FindNeighbors(pbmc, dims = 1:10) # step 1
# dims = dimensions of reduction to use as input
pbmc <- FindClusters(pbmc, resolution = 0.5) # step 2

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
# AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1
#                2                3                2                1
# AAACCGTGTATGCG-1
#                6
# Levels: 0 1 2 3 4 5 6 7 8

# Non-linear dimensional reduction (UMAP/tSNE) -----------------------------------------
# Visualize groupings of cells. Should be similar to above clusters.
# Pros: learn about underlying structures in the dataset.
# Cons: Preserve local relationships (2 similar cells w/ be near each other),
#       but often lose global relationships.
#       Therefore, don't draw any biological interpretations from these plots!
pbmc <- RunUMAP(pbmc, dims = 1:10)
umapPlot <- DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave(plot = umapPlot,
       filename = "umapPlot.pdf", path = plotsDir,
       width = 90, height = 90, units = "mm")

# Save pbmc object -----------------------------------------
saveRDS(pbmc, file = file.path(rDataDir, "pbmc.rds"))

# Load pbmc object -----------------------------------------
pbmc <- ReadRDS(file = file.path(rDataDir, "pbmc.rds"))

# Finding differentially expressed features (markers) -----------------------------------------
# Define clusters by DEGs.
# Default: + and - markers for a single cluster (ident.1) compared to all other cells.
# Can also test groups of cells against all other cells.
# Uses the `presto` package to speed up DE testing.

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
# logfc.threshold = 0.1, min.pct = 0.01
# logfc.threshold: Limit testing to genes w/ ata least X-fold difference between two groups of cells, on average
# min.pct: Limit testing to genes in a min. fraction of cells in either of two populations.
head(cluster2.markers, n = 5)
#             p_val avg_log2FC pct.1 pct.2    p_val_adj
# IL32 2.892340e-90  1.3070772 0.947 0.465 3.966555e-86
# LTB  1.060121e-86  1.3312674 0.981 0.643 1.453850e-82
# CD3D 8.794641e-71  1.0597620 0.922 0.432 1.206097e-66
# IL7R 3.516098e-68  1.4377848 0.750 0.326 4.821977e-64
# LDHB 1.642480e-67  0.9911924 0.954 0.614 2.252497e-63
# ‘pct.1’: The percentage of cells where the gene is detected in the first group
# ‘pct.2’: The percentage of cells where the gene is detected in the second group
# 'p_val_adj': Adjusted p-value, based on bonferroni using all genes in the dataset

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
#                       p_val avg_log2FC pct.1 pct.2     p_val_adj
# FCGR3A        8.246578e-205   6.794969 0.975 0.040 1.130936e-200
# IFITM3        1.677613e-195   6.192558 0.975 0.049 2.300678e-191
# CFD           2.401156e-193   6.015172 0.938 0.038 3.292945e-189
# CD68          2.900384e-191   5.530330 0.926 0.035 3.977587e-187
# RP11-290F20.3 2.513244e-186   6.297999 0.840 0.017 3.446663e-182

# find markers for every cluster compared to all remaining cells,
# report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
#        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene
#  1 3.75e-112       1.21 0.912 0.592 5.14e-108 0       LDHB
#  2 9.57e- 88       2.40 0.447 0.108 1.31e- 83 0       CCR7
#  3 1.15e- 76       1.06 0.845 0.406 1.58e- 72 0       CD3D
#  4 1.12e- 54       1.04 0.731 0.4   1.54e- 50 0       CD3E
#  5 1.35e- 51       2.14 0.342 0.103 1.86e- 47 0       LEF1
#  6 1.94e- 47       1.20 0.629 0.359 2.66e- 43 0       NOSIP
#  7 2.81e- 44       1.53 0.443 0.185 3.85e- 40 0       PIK3IP1
#  8 6.27e- 43       1.99 0.33  0.112 8.60e- 39 0       PRKCQ-AS1
#  9 1.16e- 40       2.70 0.2   0.04  1.59e- 36 0       FHIT
# 10 1.34e- 34       1.96 0.268 0.087 1.84e- 30 0       MAL

# There are many different DE test options (e.g., ROC, DESeq2, wilcox_limma, LR, etc.).
# Default is 'wilcox'.
# The ROC test returns the ‘classification power’ for any individual
#      marker (ranging from 0 - random, to 1 - perfect).
cluster0.markers <- FindMarkers(
  pbmc,
  ident.1 = 0,
  logfc.threshold = 0.25,
  test.use = "roc",
  only.pos = TRUE
)

# Visualize marker expression =============
# Recommended functions: VlnPlot(), FeaturePlot(), RidgePlot(), CellScatter(), DotPlot()
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE) # raw counts
FeaturePlot(pbmc,
  features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
)
# Heatmap of  top 20 markers (or all markers if less than 20) for each cluster
top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
markerHeatmap <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()

pdf(file.path(plotsDir, "markerHeatmap.pdf"), width = 10, height = 10)
print(markerHeatmap)
dev.off()

# Assigning cell type identity to clusters -----------------------------------------
# Use canonical markers...
new.cluster.ids <- c(
  "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
  "FCGR3A+ Mono", "NK", "DC", "Platelet"
)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

umapCanonical <- DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE,
  label.size = 4.5
) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(
  plot = umapCanonical,
  filename = file.path(plotsDir, "umapCanonical.pdf"),
  height = 7, width = 12
)

# Save rds -----------------------------------------
saveRDS(pbmc, file = file.path(rDataDir, "pbmc3k_final.rds"))
