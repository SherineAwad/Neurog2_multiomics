#!/usr/bin/env Rscript

library(ArchR)
library(argparse)
library(BSgenome.Mmusculus.UCSC.mm10)


# -------------------------------
# Arguments
# -------------------------------
parser <- ArgumentParser(description = "ArchR Peak Analysis for Neurog2")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
parser$add_argument("--threads", type = "integer", default = 8, help = "Number of threads")

args <- parser$parse_args()

# -------------------------------
# Setup
# -------------------------------
addArchRThreads(threads = args$threads)
set.seed(123)

cat("🔧 Parameters:\n")
cat("  Project:", args$project, "\n")
cat("  Threads:", args$threads, "\n\n")

# -------------------------------
# Load Project
# -------------------------------
cat("📁 Loading project...\n")
proj <- loadArchRProject(path = args$project, force = TRUE)

cat("✅ Project loaded successfully!\n")
cat("   Cells:", nCells(proj), "\n")
cat("   Samples:", paste(unique(proj$Sample), collapse = ", "), "\n\n")

# -------------------------------
# Step 1: Call Peaks using Cluster_Annotation
# -------------------------------
cat("1️⃣ ADDING GROUP COVERAGES using Cluster_Annotation...\n")
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Cluster_Annotation",
  minCells = 50,
  force = TRUE
)


cat("2️⃣ CALLING REPRODUCIBLE PEAKS...\n")
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Cluster_Annotation",
  path = file.path(getOutputDirectory(proj), "PeakCalls"),
  peakMethod = "Macs2",
  pathToMacs2 = "/nfs/turbo/umms-thahoang/sherine/miniconda/envs/archr_env/bin/macs2",
  cutOff = 0.05,
  additionalParams = "--nomodel --shift -100 --extsize 200",
  force = TRUE
)


cat("3️⃣ ADDING PEAK MATRIX...\n")
proj <- addPeakMatrix(proj, force = TRUE)

# Print peak information
peak_set <- getPeakSet(proj)
cat("✅ Peaks called:", length(peak_set), "\n\n")

# -------------------------------
# Step 2: Identify Marker Peaks between Clusters
# -------------------------------
cat("4️⃣ IDENTIFYING MARKER PEAKS between Cluster_Annotation...\n")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_Annotation",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 1000
)

# -------------------------------
# Step 3: Identify Differential Peaks between Control vs OE
# -------------------------------
cat("5️⃣ IDENTIFYING DIFFERENTIAL PEAKS between Control vs OE...\n")
markersPeaks_condition <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "OE",
  bgdGroups = "Control"
)

# -------------------------------
# Step 4: Identify Condition-Specific Peaks within Cell Types
# -------------------------------
cat("6️⃣ IDENTIFYING CONDITION-SPECIFIC PEAKS within cell types...\n")

# Use Cluster_Annotation for cell type-specific analysis
cell_types <- unique(proj$Cluster_Annotation)
markersPeaks_celltype_condition <- list()

for (cell_type in cell_types) {
  cat("   Processing:", cell_type, "\n")
  
  cells_subset <- getCellNames(proj)[which(proj$Cluster_Annotation == cell_type)]
  
  if (length(cells_subset) >= 100) { # Only analyze if sufficient cells
    proj_sub <- proj[cells_subset, ]
    
    # Check if both conditions are present
    condition_counts <- table(proj_sub$Sample)
    if (length(condition_counts) == 2 && min(condition_counts) >= 20) {
      markers <- getMarkerFeatures(
        ArchRProj = proj_sub,
        useMatrix = "PeakMatrix",
        groupBy = "Sample",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon",
        useGroups = "OE",
        bgdGroups = "Control"
      )
      markersPeaks_celltype_condition[[cell_type]] <- markers
    }
  }
}

# -------------------------------
# Step 5: Save Results
# -------------------------------
cat("7️⃣ SAVING RESULTS...\n")
output_dir <- file.path(getOutputDirectory(proj), "PeakAnalysis_Results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# A. Save cluster marker peaks
markerList_cluster <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
for (cluster_name in names(markerList_cluster)) {
  if (nrow(markerList_cluster[[cluster_name]]) > 0) {
    filename <- file.path(output_dir, paste0("MarkerPeaks_Cluster_", cluster_name, ".csv"))
    write.csv(as.data.frame(markerList_cluster[[cluster_name]]), filename, row.names = FALSE)
  }
}

# B. Save condition differential peaks
markerList_condition <- getMarkers(markersPeaks_condition, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
if (length(markerList_condition) > 0 && nrow(markerList_condition[[1]]) > 0) {
  filename <- file.path(output_dir, "DifferentialPeaks_OE_vs_Control_AllCells.csv")
  write.csv(as.data.frame(markerList_condition[[1]]), filename, row.names = FALSE)
  cat("   Saved overall OE vs Control peaks:", nrow(markerList_condition[[1]]), "\n")
}

# C. Save cell type-specific condition peaks
for (cell_type in names(markersPeaks_celltype_condition)) {
  marker_list <- getMarkers(markersPeaks_celltype_condition[[cell_type]], cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
  if (length(marker_list) > 0 && nrow(marker_list[[1]]) > 0) {
    filename <- file.path(output_dir, paste0("DifferentialPeaks_OE_vs_Control_", cell_type, ".csv"))
    write.csv(as.data.frame(marker_list[[1]]), filename, row.names = FALSE)
    cat("   Saved", cell_type, "specific peaks:", nrow(marker_list[[1]]), "\n")
  }
}

# -------------------------------
# Step 6: Create Visualizations
# -------------------------------
cat("8️⃣ CREATING VISUALIZATIONS...\n")

# A. Heatmap for cluster marker peaks
tryCatch({
  heatmap_peaks <- plotMarkerHeatmap(
    seMarker = markersPeaks,
    cutOff = "FDR <= 0.1 & Log2FC >= 1",
    transpose = TRUE
  )
  pdf(file.path(output_dir, "MarkerPeaks_Heatmap.pdf"), width = 12, height = 10)
  print(heatmap_peaks)
  dev.off()
  cat("   Saved cluster marker peaks heatmap\n")
}, error = function(e) cat("   Could not create cluster heatmap:", e$message, "\n"))

# B. Volcano plot for condition comparison
tryCatch({
  volcano_data <- metadata(markersPeaks_condition)$MarkerList$OE
  if (!is.null(volcano_data) && nrow(volcano_data) > 0) {
    pdf(file.path(output_dir, "VolcanoPlot_OE_vs_Control.pdf"), width = 10, height = 8)
    plot(volcano_data$Log2FC, -log10(volcano_data$FDR), 
         xlab = "Log2 Fold Change", ylab = "-log10(FDR)",
         main = "OE vs Control - Differential Peaks",
         pch = 19, col = ifelse(volcano_data$FDR < 0.1 & abs(volcano_data$Log2FC) > 0.5, "red", "grey"))
    abline(h = -log10(0.1), v = c(-0.5, 0.5), lty = 2)
    legend("topright", legend = c("FDR < 0.1 & |Log2FC| > 0.5"), col = "red", pch = 19)
    dev.off()
    cat("   Saved volcano plot\n")
  }
}, error = function(e) cat("   Could not create volcano plot:", e$message, "\n"))

# -------------------------------
# Step 7: Save Project
# -------------------------------
cat("9️⃣ SAVING UPDATED PROJECT...\n")
proj <- saveArchRProject(ArchRProj = proj)

# -------------------------------
# Final Summary
# -------------------------------
cat("\n✅ ANALYSIS COMPLETE!\n")
cat("📊 SUMMARY:\n")
cat("   - Peaks called:", length(peak_set), "\n")
cat("   - Output directory:", output_dir, "\n")
cat("   - Cluster comparisons: Cluster_Annotation (15 clusters)\n")
cat("   - Condition comparisons: OE vs Control\n")
cat("   - Cell type comparisons:", paste(cell_types, collapse = ", "), "\n")
cat("   - Available matrices:", paste(getAvailableMatrices(proj), collapse = ", "), "\n")

cat("\n🎯 NEXT STEPS:\n")
cat("   - Check CSV files in output directory for differential peaks\n")
cat("   - Use marker peaks for motif enrichment analysis\n")
cat("   - Integrate with your existing DGE results\n")
