# ############################################################
# 1. PACKAGES
# Load packages used for normalization and diagnostic plots.
# ############################################################
library(Seurat)
library(sparseMatrixStats)
library(ggplot2)

# ############################################################
# 2. SAMPLE AND PATH SETTINGS
# Define the input samples and the normalization output directory.
# ############################################################
#"C2A1", "C2D1", "C4D1", "C1A1"
sample_bank <- c(
  "C1A1", "C4D1"
)

input_base <- "C:/Spatial Transcriptomics/1_QC"
output_base <- "C:/Spatial Transcriptomics/2_Normalization"

dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
for (sample_name in sample_bank) {

  input_file <- file.path(
    input_base,
    sample_name,
    "Data",
    paste0(sample_name, "_QC_filtered.rds")
  )

  sample_dir <- file.path(output_base, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  # ############################################################
  # 3. LOAD QC DATA AND PLOT RAW COUNTS
  # Load the sample and visualize the mean-variance relationship in raw counts.
  # ############################################################
  visium <- readRDS(input_file)
  rawCts <- as.matrix(visium[["Spatial"]]$counts)
  df <- data.frame(
    mean = rowMeans(rawCts),
    var = matrixStats::rowVars(rawCts)
  )

  raw <- ggplot(df, aes(x = mean, y = var)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = "red"
    ) +
    theme_classic()

  ggsave(
    file.path(sample_dir, "rawMeanVariancePlot.png"),
    plot = raw,
    width = 10,
    height = 5,
    dpi = 300
  )

  # ############################################################
  # 4. RUN SCTransform AND PLOT SCALED DATA
  # Apply SCTransform while regressing mitochondrial percentage, then inspect scaled data.
  # ############################################################
  visium <- SCTransform(
    visium,
    assay = "Spatial",
    vars.to.regress = "percentMt_Spatial",
    verbose = FALSE
  )

  sctScaleData <- LayerData(
    visium,
    assay = "SCT",
    layer = "scale.data"
  )

  df <- data.frame(
    mean = rowMeans(sctScaleData),
    var = matrixStats::rowVars(sctScaleData)
  )

  SCTscale <- ggplot(df, aes(x = mean, y = var)) +
    geom_point(
      alpha = 0.4,
      size = 0.5
    ) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(
      intercept = 0,
      slope = 1,
      colour = "red"
    ) +
    theme_classic()

  ggsave(
    file.path(sample_dir, "SCT_MeanVariance.png"),
    plot = SCTscale,
    width = 10,
    height = 5,
    dpi = 300
  )

  # ############################################################
  # 5. NORMALIZE SPATIAL ASSAY AND SAVE
  # Normalize and scale the Spatial assay, then save the resulting object.
  # ############################################################
  DefaultAssay(visium) <- "Spatial"

  visium <- NormalizeData(
    visium,
    assay = "Spatial",
    verbose = FALSE
  )

  visium <- ScaleData(
    visium,
    assay = "Spatial",
    vars.to.regress = "percentMt_Spatial",
    verbose = FALSE
  )

  saveRDS(
    visium,
    file = file.path(
      sample_dir,
      paste0(sample_name, "_normalized.rds")
    ))}
