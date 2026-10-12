# ############################################################
# 1. PACKAGES
# Load the R packages used for spatial quality control and plotting.
# ############################################################
library(Seurat)
library(paletteer)
library(ggplot2)
library(sparseMatrixStats)
library(patchwork)

# ############################################################
# 2. SAMPLE AND QUALITY-CONTROL SETTINGS
# Define the samples, input/output folders, and filtering thresholds.
# ############################################################
sample_bank <- c("V1", "V4", "V5", "V6", "V7","V8",
                 "C1A1", "C1D1", "C2A1", "C2D1", "C3D1",
                 "C4A1", "C4D1", "C5A1", "C5D1")

samples_to_run <- c("C1D1")

input_base <- "C:/Spatial Transcriptomics/SpaceRanger"
output_base <- "C:/Spatial Transcriptomics/1_QC"

min_counts <- 1000
max_counts <- 40000

min_features <- 500
max_features <- 8000

max_mt <- 20

# ############################################################
# 3. LOAD SAMPLE AND CALCULATE QC METRICS
# Load each sample and calculate counts, detected features, and mitochondrial percentage.
# ############################################################
for (sample_name in samples_to_run) {

  sample_dir <- file.path(output_base, sample_name)
  figure_dir <- file.path(sample_dir, "Figures")
  data_dir <- file.path(sample_dir, "Data")
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  input_dir <- file.path(input_base, sample_name, "outs")
  visium <- Load10X_Spatial(
    data.dir = input_dir,
    filename = "filtered_feature_bc_matrix.h5"
  )

  visium[["nCount_Spatial"]] <- colSums(visium[["Spatial"]]$counts)
  visium[["nFeature_Spatial"]] <- colSums(visium[["Spatial"]]$counts > 0)
  visium[["percentMt_Spatial"]] <- PercentageFeatureSet(
    visium,
    assay = "Spatial",
    pattern = "^MT-"
  )
  visium$percentMt_Spatial[is.na(visium$percentMt_Spatial)] <- 0

  # ############################################################
  # 4. GENERATE PRE-FILTER QC FIGURES
  # Plot quality metrics and their distributions before applying the thresholds.
  # ############################################################
  QCvln <- VlnPlot(
    visium,
    features = c(
      "nCount_Spatial",
      "nFeature_Spatial",
      "percentMt_Spatial"
    )
  )

  ggsave(
    file.path(figure_dir, "QC_violin_plot.png"),
    plot = QCvln,
    width = 10,
    height = 5,
    dpi = 300
  )

  CorrPlot <- ggplot(
    visium@meta.data,
    aes(x = nCount_Spatial, y = nFeature_Spatial)
  ) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    labs(
      x = "Total counts",
      y = "Detected genes"
    )

  ggsave(
    file.path(figure_dir, "CorrPlot.png"),
    plot = CorrPlot,
    width = 10,
    height = 5,
    dpi = 300
  )

  SpatialpreQC <- SpatialPlot(
    visium,
    features = c(
      "nCount_Spatial",
      "nFeature_Spatial",
      "percentMt_Spatial"
    ),
    ncol = 2,
    image.alpha = 0.2
  )

  ggsave(
    file.path(figure_dir, "SpatialpreQC.png"),
    plot = SpatialpreQC,
    width = 10,
    height = 5,
    dpi = 300
  )

  png(
    file.path(figure_dir, "QC_cutoffs.png"),
    width = 10,
    height = 12,
    units = "in",
    res = 300
  )

  par(mfrow = c(3, 1))

  hist(
    visium$nCount_Spatial,
    breaks = 50,
    main = "Distribution of total counts",
    xlab = "Total counts"
  )

  abline(v = min_counts, col = "red", lwd = 2)
  abline(v = quantile(visium$nCount_Spatial, 0.01), col = "blue", lwd = 2)
  abline(v = max_counts, col = "red", lwd = 2)

  abline(
    v = quantile(
      visium$percentMt_Spatial,
      0.99,
      na.rm = TRUE
    ),
    col = "blue",
    lwd = 2
  )
  hist(
    visium$nFeature_Spatial,
    breaks = 50,
    main = "Distribution of detected genes",
    xlab = "Number of detected genes"
  )

  abline(v = min_features, col = "red", lwd = 2)
  abline(v = quantile(visium$nFeature_Spatial, 0.01), col = "blue", lwd = 2)
  abline(v = max_features, col = "red", lwd = 2)

  hist(
    visium$percentMt_Spatial,
    breaks = 50,
    main = "Distribution of mitochondrial percentage",
    xlab = "Percentage of mitochondrial genes"
  )
  abline(
    v = quantile(
      visium$percentMt_Spatial,
      0.99,
      na.rm = TRUE
    ),
    col = "blue",
    lwd = 2
  )
  dev.off()

  # ############################################################
  # 5. APPLY FILTERS AND REVIEW SPATIAL QC
  # Flag spots that meet any configured filtering threshold and plot the flags.
  # ############################################################
  visium$nCountFiltered <-
    visium$nCount_Spatial <= min_counts |
    visium$nCount_Spatial >= max_counts

  visium$percentMtFiltered <-
    visium$percentMt_Spatial >= max_mt

  visium$nFeatureFiltered <-
    visium$nFeature_Spatial <= min_features |
    visium$nFeature_Spatial >= max_features

  visium$QC_Filter <-

    visium$nCountFiltered |
    visium$percentMtFiltered |
    visium$nFeatureFiltered

  SpatialpostQC <- SpatialDimPlot(
    visium,
    group.by = c(
      "nCountFiltered",
      "percentMtFiltered",
      "nFeatureFiltered",
      "QC_Filter"
    ),
    cols = c("lightgrey", "red"),
    ncol = 2,
    image.alpha = 0.2
  )

  ggsave(
    file.path(figure_dir, "SpatialpostQC.png"),
    plot = SpatialpostQC,
    width = 10,
    height = 5,
    dpi = 300
  )

  # ############################################################
  # 6. SAVE FILTERED DATA AND QC SUMMARY
  # Save the QC object, filtered object, and per-sample filtering summary.
  # ############################################################
  filtered_visium <- subset(
    visium,
    subset = QC_Filter == FALSE
  )

  saveRDS(
    visium,
    file = file.path(data_dir, paste0(sample_name, "_filtered.rds"))
  )

  write.csv(
    data.frame(
      sample = sample_name,
      spots_before = ncol(visium),
      spots_removed = sum(visium$QC_Filter),
      spots_after = ncol(filtered_visium),
      percent_removed = 100 * sum(visium$QC_Filter) / ncol(visium)
    ),
    file = file.path(sample_dir, "QC_summary.csv"),
    row.names = FALSE
  )

  saveRDS(
    filtered_visium,
    file = file.path(data_dir, paste0(sample_name, "_QC_filtered.rds"))
  )
}
