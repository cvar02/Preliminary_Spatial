# ############################################################
# 1. PACKAGES
# Load packages used for dimensionality reduction and plotting.
# ############################################################
library(Seurat)
library(ggplot2)
library(patchwork)

# ############################################################
# 2. SAMPLE, PARAMETERS, AND PATHS
# Select a sample and set the dimensionality-reduction parameters and folders.
# ############################################################
# Run one sample at a time because suitable clustering parameters may differ between samples.

sample_bank <- c(
  "C2A1", "C2D1", "C4D1", "C1A1"
)
selected_sample <- "C2D1"

input_base <- "C:/Spatial Transcriptomics/2_Normalization"
output_base <- "C:/Spatial Transcriptomics/3_DimRed"

selected_pcs <- 40
selected_neighbors <- 17
selected_min_dist <- 0.05
selected_resolution <- 0.6

# ############################################################
# 3. LOAD DATA AND RUN PCA
# Load the normalized sample, run PCA, and inspect the variance plots.
# ############################################################
input_file <- file.path(
  input_base,
  selected_sample,
  paste0(selected_sample, "_normalized.rds")
)

output_dir <- file.path(
  output_base,
  selected_sample
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
visium <- readRDS(input_file)
DefaultAssay(visium) <- "SCT"

visium <- RunPCA(
  visium,
  assay = "SCT",
  npcs = 50,
  verbose = FALSE
)

elbow_plot <- ElbowPlot(
  visium,
  ndims = 50,
  reduction = "pca"
)

print(elbow_plot)

# Get the variance explained by each PC
ElbowPlot(
  visium,
  ndims = 50,
  plot_type = "cumulative_variance",
  reduction = 'pca'
) +
  geom_hline(yintercept = 90) +
  geom_vline(xintercept = 40)

# ############################################################
# 4. FIND NEIGHBORS, CLUSTERS, AND UMAP
# Build the neighbor graph, assign clusters, and calculate the final UMAP.
# ############################################################
visium <- FindNeighbors(
  visium,
  reduction = "pca",
  dims = 1:selected_pcs,
  verbose = FALSE
)

visium <- FindClusters(
  visium,
  resolution = selected_resolution,
  algorithm = 4,
  cluster.name = "final_clusters",
  verbose = FALSE
)

visium <- RunUMAP(
  visium,
  reduction = "pca",
  dims = 1:selected_pcs,
  n.neighbors = selected_neighbors,
  min.dist = selected_min_dist,
  reduction.name = "umap_final",
  verbose = FALSE
)

# ############################################################
# 5. PLOT, RECORD, AND SAVE RESULTS
# Visualize clustering in UMAP and tissue space, then save plots and outputs.
# ############################################################
umap_final <- DimPlot(
  visium,
  reduction = "umap_final",
  group.by = "final_clusters",
  label = TRUE
) +
  ggtitle(
    paste0(
      selected_sample,
      " | PCs=", selected_pcs,
      " | NN=", selected_neighbors,
      " | min.dist=", selected_min_dist,
      " | resolution=", selected_resolution
    )
  )

print(umap_final)

spatial_final <- SpatialDimPlot(
  visium,
  group.by = "final_clusters",
) +
  ggtitle(
    paste0(
      selected_sample,
      " spatial clusters"
    )
  )

print(spatial_final)

final_plots <- umap_final + spatial_final

ggsave(
  file.path(output_dir, "final_clustering.png"),
  final_plots,
  width = 14,
  height = 6,
  dpi = 300
)

visium@misc$final_clustering_parameters <- list(
  PCs = selected_pcs,
  n_neighbors = selected_neighbors,
  min_dist = selected_min_dist,
  resolution = selected_resolution
)

saveRDS(
  visium,
  file.path(
    output_dir,
    paste0(selected_sample, "_DimRed.rds")
  )
)

write.csv(
  data.frame(
    cluster = names(table(visium$final_clusters)),
    spots = as.numeric(table(visium$final_clusters))
  ),
  file.path(output_dir, "final_cluster_sizes.csv"),
  row.names = FALSE
)
