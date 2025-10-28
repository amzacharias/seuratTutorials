#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Analyzing VISIUM Spatial Transcriptomics Data from 10x Genomics
# Author: Amanda Zacharias
# Date: 2025-10-28
# Email: 16amz1@queensu.ca
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.5.0
# https://satijalab.org/seurat/articles/spatial_vignette
# Data are sagittal sections from mouse brain (10x Genomics Visium v1)
#
# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# Source -----------------------------------------------

# Pathways -----------------------------------------------
## Input ===========

## Output ===========
proj_dir <- file.path("spatial_seqbased")
results_dir <- file.path(proj_dir, "results")
plots_dir <- file.path(results_dir, "plots")
rdata_dir <- file.path(results_dir, "rdata")

# Load data -----------------------------------------------
## spRNA-seq data ===========
# InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
brain2 <- LoadData("stxBrain", type = "posterior1")

# Visium data are stored as follows:
# 1. Spot by gene expression matrix (RNA assay slot)
#       columns = spots, not single-cells
# 2. Image of the tissue section (image slot)
# 3. Scaling factors relating high resolution image to lower resolution image

## scRNA-seq data ===========
# download.file(
#   url = "https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1",
#   destfile = file.path(proj_dir, "data", "allen_cortex.rds")
# )
allen_reference <- readRDS(file.path(proj_dir, "data", "allen_cortex.rds"))

# Preprocessing -----------------------------------------------
# Similar pre-processing to scRNA-seq.
# 1. Normalize data to account for sequencing depth differences across spots.
#   - This variance can be quite high for spots.
ncount_vln_plot <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
ncount_sp_plot <- SpatialPlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
ncount_plot <- wrap_plots(ncount_vln_plot, ncount_sp_plot)
ggsave(
  plot = ncount_plot, filename = "visium_ncount.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

# Depth appears to be influenced by tissue anatomy.
#   Spots with white matter have lower depth.
# Therefore, use a method that can hopefully minimize the technical and preserve bioloigical variation.
# Using SCTransform (Hafemeister and Satija, 2019) instead of typical LogNormalize().
# Tutorial demonstrates that LogNormalize fails for highly expressed genes.
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# Visualize gene expression of interest -----------------------------------------------
# Overlay molecular data on tissue histology.
# Hpca = marker of hippocampus
# Ttr marker of choroid plexus
## Emphasize the molecular data ===========
genes_of_interest <- c("Hpca", "Ttr")
for (gene in genes_of_interest) {
  p <- SpatialFeaturePlot(brain, features = gene, pt.size.factor = 1.6, alpha = c(1, 1))
  ggsave(
    plot = p, filename = sprintf("visium_%s_expr_emphasis.pdf", gene),
    path = plots_dir, width = 180, height = 100, units = "mm"
  )
}
## Emphasize the histology ===========
for (gene in genes_of_interest) {
  p <- SpatialFeaturePlot(brain, features = gene, pt.size.factor = 1, alpha = c(0.1, 1))
  ggsave(
    plot = p, filename = sprintf("visium_%s_hist_emphasis.pdf", gene),
    path = plots_dir, width = 180, height = 100, units = "mm"
  )
}

# Dimensional reduction and clustering -----------------------------------------------
# Same workflow as scRNA-seq.
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

umap_dim_plot <- DimPlot(brain, reduction = "umap", label = TRUE)
# facet.highlight makes clusters easier to see on spatial plot
# could also highlight specific clusters by providing a vector of identities to highlight (`CellsByIdentities()`)
umap_spatial_plot <- SpatialDimPlot(brain, label = TRUE, label.size = 3, facet.highlight = TRUE)
umap_plot <- wrap_plots(umap_dim_plot, umap_spatial_plot)
ggsave(
  plot = umap_plot, filename = "visium_umap.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

# It's nice that the clustering of spots seems to correspond to anatomical regions of the brain.

# Interactive visualization -----------------------------------------------
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
# for fun looked at Arntl, looks to be expressed everywhere, no specific region.
# were samples taken in the morning?
LinkedDimPlot(brain) # based on UMAP plot, highlight corresponding region on spatial plot

# Identify spatially variable features -----------------------------------------------
# Molecular features that correlate with location.
# 2 methods:
#     1. Differential expression based on pre-annotated anatomical regions.
#         - works well here b.c. cell clusters appear anatomically distinct.
#     2. Unsupervised identification of spatially variable features.
# Other methods exist including SpatialDE and Splotch.

## Method 1: DE based on clusters ===========
# E.g., what differentiates clusters 5 and 6?
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
de_markers_plot <- SpatialFeaturePlot(
  brain, features = rownames(de_markers)[1:3], ncol = 3, alpha = c(0.1, 1)
)
ggsave(
  plot = de_markers_plot, filename = "visium_de_markers.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

## Method 2: Unsupervised spatially variable features ===========
brain <- FindSpatiallyVariableFeatures(
  brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], selection.method = "moransi"
)
var_markers <- head(SpatiallyVariableFeatures(brain, method = "moransi"), 6)
var_markers_plot <- SpatialFeaturePlot(
  brain, features = var_markers, ncol = 3, alpha = c(0.1, 1)
)
ggsave(
  plot = var_markers_plot, filename = "visium_var_markers.pdf",
  path = plots_dir, width = 180, height = 200, units = "mm"
)

# Subset anatomical regions -----------------------------------------------
# Subset the ~ frontal cortex region. Allows integration with cortical scRNA-seq data.
# Requires using the column ("x") and row ("y") positions from
#    the image data, and adding them to metadata.

# Clusters of interest
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))

# Get spatial coordinates
centroids <- cortex[["anterior1"]]@boundaries$centroids
coords <- setNames(as.data.frame(centroids@coords), c("x", "y"))
rownames(coords) <- centroids@cells

# Add coords to the cortex metadata
cortex$x <- coords[colnames(cortex), "x"]
cortex$y <- coords[colnames(cortex), "y"]

# Spatial subsetting by image position.
# Get upper right quad
#    Invert means we'll get spots >= 2719 y and 7835 <= x.
cortex <- subset(cortex, y < 2719 | x > 7835, invert = TRUE)
# Remove lower right quad
#  Spots above the calculated descending line are removed.
m <- (9395 - 5747) / (4960 - 7715) # slope
b <- 5747 - m * 7715 # intercept, passes through x = 7715, y = 5747
cortex <- subset(cortex, y > m * x + b, invert = TRUE)

# Visualize
crop_spatial_plot <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 4)
full_spatial_plot <- SpatialDimPlot(
  cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3
)
crop_full_plot <- wrap_plots(crop_spatial_plot, full_spatial_plot)
ggsave(
  plot = crop_full_plot, filename = "visium_cortex_subset.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

# Integrate with scRNA-seq data -----------------------------------------------
## Normalize scRNA data ===========
# Speed up by learning noise on only 3000 cells. Apparently normalization is still good.
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) |>
  RunPCA(verbose = FALSE) |>
  RunUMAP(dims = 1:30, verbose = FALSE)
# Plot UMAP
# annotation is stored in 'subclass' metadata column for scRNA-seq
allen_umap_plot <- DimPlot(
  allen_reference, group.by = "subclass", label = TRUE
)
ggsave(
  plot = allen_umap_plot, filename = "visium_allen_umap.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

## Renormalize spatial data ===========
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) |>
  RunPCA(verbose = FALSE)

## Transfer cell-type annotations from allen brain atlas to spatial ===========
anchors <- FindTransferAnchors(
  reference = allen_reference,
  query = cortex,
  normalization.method = "SCT"
)
predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = allen_reference$subclass,
  prediction.assay = TRUE,
  weight.reduction = cortex[["pca"]],
  dims = 1:30
)
cortex[["predictions"]] <- predictions.assay

## Visualize predicted cell types ===========
DefaultAssay(cortex) <- "predictions"
# Explore laminar excitatory neurons in frontal cortex, sub-types
laminar_neurons_plot <- SpatialFeaturePlot(
  cortex,
  features = c("L2/3 IT", "L4"),
  pt.size.factor = 1.6,
  ncol = 2,
  crop = TRUE
)
ggsave(
  plot = laminar_neurons_plot, filename = "visium_laminar_neurons.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)
# Predict cell types based on cell type prediction scores, rather than gene expression
cortex <- FindSpatiallyVariableFeatures(
  cortex, assay = "predictions", selection.method = "moransi",
  features = rownames(cortex), r.metric = 5, slot = "data"
)
top.clusters <- head(SpatiallyVariableFeatures(cortex, method = "moransi"), 4)
pred_cells_plot <- SpatialPlot(
  object = cortex, features = top.clusters, ncol = 2
)
ggsave(
  plot = pred_cells_plot, filename = "visium_predicted_celltypes.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)
# Show spatial location of neuronal and non-neuronal cell types
spatial_neun_nonneuon_plot <- SpatialFeaturePlot(
  cortex,
  features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", "L6b", "Oligo"),
  pt.size.factor = 1,
  ncol = 2,
  crop = FALSE,
  alpha = c(0.1, 1)
)
ggsave(
  plot = spatial_neun_nonneuon_plot, filename = "visium_neun_nonneun.pdf",
  path = plots_dir, width = 180, height = 600, units = "mm"
)

# Working with multiple slices -----------------------------------------------
## Normalize a second slice from the posterior region of the brain ===========
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

## Merge with anterior slice ===========
brain.merge <- merge(brain, brain2)

## Joint dim redution and clustering ===========
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- brain.merge |>
  RunPCA(verbose = FALSE) |>
  FindNeighbors(reduction = "pca", dims = 1:30) |>
  FindClusters(verbose = FALSE) |>
  RunUMAP(reduction = "pca", dims = 1:30)

## Visualizations ===========
multi_dim_plot <- DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
multi_spatial_plot <- SpatialDimPlot(brain.merge)
multi_ft_plot <- SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))

multi_plot <- wrap_plots(
  multi_dim_plot, multi_spatial_plot, multi_ft_plot, ncol = 1
) + plot_layout(heights = c(1, 1, 2))
ggsave(
  plot = multi_plot, filename = "visium_multi_desc.pdf",
  path = plots_dir, width = 180, height = 300, units = "mm"
)

# Save final objects -----------------------------------------------
saveRDS(brain, file = file.path(rdata_dir, "visium_brain.rds"))
saveRDS(cortex, file = file.path(rdata_dir, "visium_cortex.rds"))
saveRDS(allen_reference, file = file.path(rdata_dir, "visium_allen_reference.rds"))
saveRDS(brain.merge, file = file.path(rdata_dir, "visium_multi_slice.rds"))
