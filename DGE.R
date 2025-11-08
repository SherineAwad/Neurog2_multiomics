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
# Differential Gene Expression - FIXED BACKGROUND
# ----------------------------
features <- getMarkerFeatures(
  ArchRProj = myProject,
  useMatrix = "GeneExpressionMatrix",
  groupBy = groupBy,
  useGroups = "OE",
  bgdGroups = "Control",
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

# ONE HEATMAP WITH FULL LABELS
HeatmapObj <- plotMarkerHeatmap(
  seMarker = subsetSE,
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5",
  plotLog2FC = TRUE,
  labelRows = TRUE,
  nLabel = 50
)

draw(HeatmapObj, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

# ----------------------------
# Save DGEs to CSV
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

filename <- file.path(project_path, "DGEs_Control_OE.csv")
write.csv(df, filename, row.names = FALSE)

# ----------------------------
# Save ArchR project
# ----------------------------
saveArchRProject(ArchRProj = myProject, outputDirectory = project_path, load = FALSE)

cat("✅ DGE analysis complete. Heatmap and CSV saved.\n")
