#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(future)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(dplyr)
  library(pheatmap)
  library(chromVARmotifs)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(argparse)
  library(ComplexHeatmap)
})

# ----------------------------
# Setup genome
# ----------------------------
addArchRGenome("mm10")

# ----------------------------
# Argument parsing
# ----------------------------
parser <- ArgumentParser(description = "ArchR Differential Gene Expression: Control vs OE")
parser$add_argument("-p", "--project", required = TRUE, help = "Path to ArchR project")
args <- parser$parse_args()
project_path <- args$project

groupBy <- "Sample"  # column to compare groups

# ----------------------------
# Load ArchR project
# ----------------------------
myProject <- loadArchRProject(path = project_path, force = FALSE, showLogo = TRUE)

# ----------------------------
# Check sample groups
# ----------------------------
sample_groups <- unique(getCellColData(myProject)[[groupBy]])
message("Found sample groups: ", paste(sample_groups, collapse = ", "))
if(!all(c("Control", "OE") %in% sample_groups)){
  stop("Project must contain both 'Control' and 'OE' samples in the 'Sample' metadata column")
}

# ----------------------------
# Differential Gene Expression
# ----------------------------
# Compute DE for all groups vs background (both Control and OE will be in the matrix)
features <- getMarkerFeatures(
  ArchRProj = myProject,
  useMatrix = "GeneExpressionMatrix",
  groupBy = groupBy,
  bias = c("Gex_nUMI","Gex_nGenes"),
  testMethod = "wilcoxon"
)

# ----------------------------
# Top DE genes for heatmap
# ----------------------------
log2fc <- assays(features)$Log2FC[,1]   # first group column
fdr    <- assays(features)$FDR[,1]

top_genes <- rowData(features)$name[order(fdr, -log2fc)][1:50]
subsetSE <- features[which(rowData(features)$name %in% top_genes),]

figure_name <- file.path(project_path, paste0("DEG_Heatmap.pdf"))
pdf(file = figure_name, width = 12)

# FIX: Create heatmap with all row labels
HeatmapObj <- plotMarkerHeatmap(
  seMarker = subsetSE,
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5",
  plotLog2FC = TRUE,
  labelRows = TRUE,  # Force all row labels
  nLabel = 50       # Label all 50 rows
)

draw(HeatmapObj, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

# Alternative: Manual heatmap with full control
figure_name2 <- file.path(project_path, paste0("DEG_Heatmap_FullLabels.pdf"))
pdf(file = figure_name2, width = 12)

# Extract the matrix for manual plotting
markerMatrix <- assays(subsetSE)$Mean
rownames(markerMatrix) <- rowData(subsetSE)$name

pheatmap(
  markerMatrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,  # Show ALL row names
  fontsize_row = 8,      # Adjust font size if needed
  main = "Differential Gene Expression - Control vs OE"
)

dev.off()

# ----------------------------
# Save DEGs to CSV
# ----------------------------
# Convert SummarizedExperiment assays to matrix and save
df <- data.frame(
  genes = rowData(features)$name,
  Log2FC = assays(features)$Log2FC,
  FDR = assays(features)$FDR,
  Mean = assays(features)$Mean,
  MeanDiff = assays(features)$MeanDiff,
  MeanBGD = assays(features)$MeanBGD,
  Pval = assays(features)$Pval
)

filename <- file.path(project_path, "DEGs_Control_OE.csv")
write.csv(df, filename, row.names = FALSE)

# ----------------------------
# Save ArchR project
# ----------------------------
saveArchRProject(ArchRProj = myProject, outputDirectory = project_path, load = FALSE)

cat("✅ DGE analysis complete. Heatmaps and CSV saved.\n")
