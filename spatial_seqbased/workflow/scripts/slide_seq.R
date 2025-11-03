#!/usr/bin/env Rscript
# -*-coding: utf-8 -*-
#-----------------------------------------------
# Title: Analyzing Slide-Seq data
# Author: Amanda Zacharias
# Date: 2025-10-28
# Email: 16amz1@queensu.ca
# Notes -----------------------------------------------
# module load StdEnv/2023 r/4.5.0
# https://satijalab.org/seurat/articles/spatial_vignette
# Data are murine hippocampus section using Slide-Seq V2 technology.
#
# Options -----------------------------------------------

# Packages -----------------------------------------------
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)

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
# InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")

# Columns are beads, rows are genes

## scRNA-seq data ===========
# download.file(
#   url = "https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=1",
#   destfile = file.path(proj_dir, "data", "mouse_hippocampus_reference.rds")
# )
ref <- readRDS(file.path(proj_dir, "data", "mouse_hippocampus_reference.rds"))

# Preprocessing -----------------------------------------------
# Preview data
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
ncount_vln_plot <- VlnPlot(slide.seq, features = "log_nCount_Spatial", pt.size = 0) + NoLegend()
ncount_sp_plot <- SpatialPlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
ncount_plot <- wrap_plots(ncount_vln_plot, ncount_sp_plot)
ggsave(
  plot = ncount_plot, filename = "slideseq_ncount.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

# Normalize, transform counts, and find clusters
slide.seq <- slide.seq |>
  SCTransform(assay = "Spatial", ncells = 3000, verbose = FALSE) |>
  RunPCA(assay = "SCT", verbose = FALSE) |>
  RunUMAP(dims = 1:30) |>
  FindNeighbors(dims = 1:30) |>
  FindClusters(resolution = 0.3, verbose = FALSE)

# Visualize clustering
clus_dim_plot <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
clust_spatial_plot <- SpatialDimPlot(slide.seq, stroke = 0)
clust_ft_spatial_plot <- SpatialDimPlot(
  slide.seq, cells.highlight = CellsByIdentities(
    object = slide.seq, idents = c(1, 6, 13)
  ), facet.highlight = TRUE
)
clust_plot <- wrap_plots(clus_dim_plot, clust_spatial_plot, clust_ft_spatial_plot) +
  plot_layout(design = "\nAABB\nCCCC\n")
ggsave(
  plot = clust_plot, filename = "slideseq_clusters.pdf",
  path = plots_dir, width = 180, height = 200, units = "mm"
)
# Integrate with scRNA-seq data -----------------------------------------------
ref <- UpdateSeuratObject(ref)

## Predict cell types for beads ===========
anchors <- FindTransferAnchors(
  reference = ref,
  query = slide.seq,
  normalization.method = "SCT",
  npcs = 50 # Number of PCs to use
)
predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = ref$celltype,
  prediction.assay = TRUE, # return an Assay object
  weight.reduction = slide.seq[["pca"]],
  dims = 1:50
)
slide.seq[["predictions"]] <- predictions.assay

## Visualize predicted cell types ===========
DefaultAssay(slide.seq) <- "predictions"
# Show prediction scores
predicted_celltypes_plot <- SpatialFeaturePlot(
  slide.seq,
  features = c(
    "Dentate Principal cells", "CA3 Principal cells", "Entorhinal cortex",
    "Endothelial tip", "Ependymal", "Oligodendrocyte"
  ),
  alpha = c(0.1, 1)
)
ggsave(
  plot = predicted_celltypes_plot, filename = "slideseq_predicted_celltypes.pdf",
  path = plots_dir, width = 180, height = 200, units = "mm"
)
# Show beads with predicted cell types
slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
predicted_cellbeads_plot <- SpatialDimPlot(
  slide.seq,
  cells.highlight = CellsByIdentities(object = slide.seq, idents = c(
    "Dentate Principal cells", "CA3 Principal cells", "Entorhinal cortex",
    "Endothelial tip", "Ependymal", "Oligodendrocyte"
  )),
  facet.highlight = TRUE
)
ggsave(
  plot = predicted_cellbeads_plot, filename = "slideseq_predicted_cellbeads.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

# Identifying spatially variable features -----------------------------------------------
# Using second method, statistics that measure dependence of a feature on a spatial location.
# Moran's I statistic is used.
# Speeding up with binning the beads into a rectangular grid, and averaging a gene and location
#   within each bin.
DefaultAssay(slide.seq) <- "SCT"

# Calculate Moran's I for top 1000 variable features
slide.seq <- FindSpatiallyVariableFeatures(
  slide.seq, assay = "SCT",
  selection.method = "moransi",
  features = VariableFeatures(slide.seq)[1:1000],
  x.cuts = 100, y.cuts = 100
)

# Visualize
var_markers_plot <- SpatialFeaturePlot(
  slide.seq,
  features = head(SpatiallyVariableFeatures(slide.seq, selection.method = "moransi"), 6),
  ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95"
)
ggsave(
  plot = var_markers_plot, filename = "slideseq_spatially_variable_features.pdf",
  path = plots_dir, width = 180, height = 180, units = "mm"
)

# Spatial deconvolution using RCTD -----------------------------------------------
# RCTD = robust cell type decomposition (sounds similar to cell-composition deconvolution tools)
# Cable et al. Nat Biotechnol. 2022. https://pmc.ncbi.nlm.nih.gov/articles/PMC8606190/
# Allows for a single spot to be a mixture of cell types.
Idents(ref) <- "celltype"

# Get gene counts and cell type clusters from reference
ref_counts <- ref[["RNA"]]$counts
ref_cluster <- as.factor(ref$celltype)
names(ref_cluster) <- colnames(ref)
ref_nUMI <- ref$nCount_RNA
names(ref_nUMI) <- colnames(ref)
reference <- spacexr::Reference(ref_counts, ref_cluster, ref_nUMI)

# Setup query for spRNA
sp_counts <- slide.seq[["Spatial"]]@counts
sp_coords <- GetTissueCoordinates(slide.seq)
colnames(sp_coords) <- c("x", "y")
sp_coords[is.na(colnames(sp_coords))] <- NULL
query <- spacexr::SpatialRNA(sp_coords, sp_counts, colSums(sp_counts))

# Annotate slide.seq with RCTD results
RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet") # 2 cell types per pixel
slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)

# Plot results
rctd_1type_plot <- SpatialDimPlot(slide.seq, group.by = "first_type")
rctd_2type_plot <- SpatialDimPlot(slide.seq, group.by = "second_type")
rctd_plot <- wrap_plots(rctd_1type_plot, rctd_2type_plot)
ggsave(
  plot = rctd_plot, filename = "slideseq_rctd_celltypes.pdf",
  path = plots_dir, width = 180, height = 100, units = "mm"
)

# In the case of one spot having two cell types, how do we analyze that downstream?

# Save final objects -----------------------------------------------
saveRDS(slide.seq, file = file.path(rdata_dir, "slideseq_hippocampus.rds"))
saveRDS(ref, file = file.path(rdata_dir, "slideseq_ref.rds"))
