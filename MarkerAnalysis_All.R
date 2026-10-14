# ############################################################
# 1. PACKAGES AND CONFIGURATION
# Identify within-sample cluster markers and visualise spatial expression.
# ############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

sample_bank <- c(
  "V1", "V4", "V5", "V6", "V7", "V8",
  "C1A1", "C1D1", "C2A1", "C2D1", "C3D1",
  "C4A1", "C4D1", "C5A1", "C5D1"
)

# Select one sample from the bank to analyse in this run.
selected_sample <- "C2D1"

if (!selected_sample %in% sample_bank) {
  stop("selected_sample must be one of: ", paste(sample_bank, collapse = ", "))
}

input_base <- "C:/Spatial Transcriptomics/3_DimRed"
output_base <- "C:/Spatial Transcriptomics/4_MarkerAnalysis"

marker_assay <- "Spatial"       # Log-normalised expression created during normalisation
marker_min_pct <- 0.25         # At least 25% of spots in either group
marker_logfc <- 0.25
marker_adj_p_cutoff <- 0.05
n_top_markers <- 5

# These are spot-level expression panels, not definitive cell-type annotations.
gene_panels <- list(
  GIST = c("KIT", "ETV1", "ANO1", "IGF1R"),
  Fibroblast_CAF = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  Myeloid = c("CD68", "C1QA", "APOE", "CD163"),
  Endothelial = c("PECAM1", "VWF", "CDH5"),
  T_Cell = c("CD3D", "CD3E", "CD2")
)

# ############################################################
# 2. HELPER FUNCTIONS
# Save plots and skip genes that are absent from a sample.
# ############################################################

save_spatial_plot <- function(plot, filename, n_features = 1L, ncol = 2L) {
  n_rows <- ceiling(n_features / ncol)
  ggsave(
    filename = filename,
    plot = plot,
    width = 6 * ncol,
    height = 5 * n_rows,
    units = "in",
    limitsize = FALSE
  )
}

plot_present_genes <- function(visium, genes, filename, assay = "Spatial", ncol = 2L) {
  present <- unique(genes[genes %in% rownames(visium[[assay]])])
  missing <- setdiff(unique(genes), present)
  if (length(missing) > 0L) {
    message("  Missing genes (", basename(filename), "): ", paste(missing, collapse = ", "))
  }
  if (length(present) == 0L) return(invisible(NULL))
  DefaultAssay(visium) <- assay
  p <- SpatialFeaturePlot(visium, features = present, ncol = ncol, image.alpha = 0.2)
  save_spatial_plot(p, filename, n_features = length(present), ncol = min(ncol, length(present)))
  invisible(present)
}

# ############################################################
# 3. ANALYSE SELECTED SAMPLE
# Retain its saved PCA, UMAP and clustering parameters.
# ############################################################

for (sample_name in selected_sample) {
  message("\nProcessing sample: ", sample_name)

  input_file <- file.path(input_base, sample_name, paste0(sample_name, "_DimRed.rds"))
  output_dir <- file.path(output_base, sample_name)
  plot_dir <- file.path(output_dir, "Plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(input_file)) {
    warning("Missing input file; skipping ", sample_name, ": ", input_file)
    next
  }

  visium <- readRDS(input_file)
  if (!"final_clusters" %in% colnames(visium[[]])) {
    warning("Missing final_clusters metadata; skipping ", sample_name)
    next
  }
  if (!marker_assay %in% Assays(visium)) {
    warning("Missing ", marker_assay, " assay; skipping ", sample_name)
    next
  }

  Idents(visium) <- "final_clusters"
  DefaultAssay(visium) <- marker_assay

  # ##########################################################
  # 3A. CLUSTER MARKER IDENTIFICATION
  # Compare each spatial cluster against all other spots in this sample.
  # ##########################################################

  markers <- FindAllMarkers(
    object = visium,
    assay = marker_assay,
    slot = "data",
    only.pos = TRUE,
    min.pct = marker_min_pct,
    logfc.threshold = marker_logfc,
    test.use = "wilcox"
  )
  write.csv(markers, file.path(output_dir, "all_cluster_markers.csv"), row.names = FALSE)

  if (nrow(markers) > 0L) {
    significant_markers <- markers %>%
      filter(!is.na(p_val_adj), p_val_adj < marker_adj_p_cutoff)

    top5 <- significant_markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = n_top_markers, with_ties = FALSE) %>%
      arrange(cluster, desc(avg_log2FC)) %>%
      ungroup()

    top_by_pct <- significant_markers %>%
      filter(pct.1 == 1, pct.2 < 0.5)
  } else {
    significant_markers <- markers
    top5 <- markers
    top_by_pct <- markers
  }

  write.csv(significant_markers, file.path(output_dir, "significant_cluster_markers.csv"), row.names = FALSE)
  write.csv(top5, file.path(output_dir, "top5_markers_per_cluster.csv"), row.names = FALSE)
  write.csv(top_by_pct, file.path(output_dir, "markers_100pct_cluster_under50pct_other.csv"), row.names = FALSE)

  # ##########################################################
  # 3B. SPATIAL CLUSTER AND TOP-MARKER PLOTS
  # Plot all clusters and the top five genes for each cluster.
  # ##########################################################

  DefaultAssay(visium) <- marker_assay
  cluster_plot <- SpatialDimPlot(
    visium,
    group.by = "final_clusters",
    label = TRUE,
    label.size = 3,
    image.alpha = 0.2
  )
  ggsave(file.path(plot_dir, "all_clusters_spatial.png"), cluster_plot, width = 10, height = 8)

  cluster_ids <- levels(Idents(visium))
  for (cluster_id in cluster_ids) {
    cluster_genes <- top5 %>%
      filter(as.character(cluster) == as.character(cluster_id)) %>%
      pull(gene)

    if (length(cluster_genes) > 0L) {
      plot_present_genes(
        visium,
        genes = cluster_genes,
        filename = file.path(plot_dir, paste0("cluster_", cluster_id, "_top5_spatial.png")),
        assay = marker_assay,
        ncol = 2L
      )
    }

    highlighted_cells <- WhichCells(visium, idents = cluster_id)
    highlight_plot <- SpatialDimPlot(
      visium,
      cells.highlight = highlighted_cells,
      image.alpha = 0.2
    ) + ggtitle(paste("Cluster", cluster_id, "spatial location"))
    ggsave(
      file.path(plot_dir, paste0("cluster_", cluster_id, "_location.png")),
      highlight_plot,
      width = 9,
      height = 7
    )
  }

  if (nrow(top_by_pct) > 0L) {
    write.csv(
      top_by_pct %>% count(cluster, name = "marker_count"),
      file.path(output_dir, "high_specificity_marker_counts.csv"),
      row.names = FALSE
    )
    plot_present_genes(
      visium,
      genes = top_by_pct$gene,
      filename = file.path(plot_dir, "high_specificity_markers_spatial.png"),
      assay = marker_assay,
      ncol = 3L
    )
  }

  # ##########################################################
  # 3C. BIOLOGICALLY INFORMED GENE PANELS
  # Show expression of relevant GIST and microenvironment genes.
  # ##########################################################

  for (panel_name in names(gene_panels)) {
    plot_present_genes(
      visium,
      genes = gene_panels[[panel_name]],
      filename = file.path(plot_dir, paste0("panel_", panel_name, "_spatial.png")),
      assay = marker_assay,
      ncol = 3L
    )
  }

  # ##########################################################
  # 3D. SPATIALLY VARIABLE GENES (MORAN'S I)
  # Test spatial autocorrelation using the SCT assay.
  # ##########################################################

  if ("SCT" %in% Assays(visium)) {
    visium <- FindSpatiallyVariableFeatures(
      visium,
      assay = "SCT",
      selection.method = "moransi",
      nfeatures = 2000
    )

    spatial_genes <- SpatiallyVariableFeatures(visium, assay = "SCT", selection.method = "moransi")
    write.csv(
      data.frame(rank = seq_along(spatial_genes), gene = spatial_genes),
      file.path(output_dir, "moransi_ranked_genes.csv"),
      row.names = FALSE
    )

    if (length(spatial_genes) > 0L) {
      plot_present_genes(
        visium,
        genes = head(spatial_genes, 12L),
        filename = file.path(plot_dir, "moransi_top12_spatial.png"),
        assay = marker_assay,
        ncol = 3L
      )
    }
  } else {
    warning("No SCT assay for Moran's I in sample ", sample_name)
  }

  # ##########################################################
  # 3E. SAVE ANALYSIS OBJECT
  # Preserve clustering and append spatial variability results.
  # ##########################################################

  saveRDS(visium, file.path(output_dir, paste0(sample_name, "_MarkerAnalysis.rds")))
  message("Completed sample: ", sample_name)
}

message("\nMarker analysis completed.")
