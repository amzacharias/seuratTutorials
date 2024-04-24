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





