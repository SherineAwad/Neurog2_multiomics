library(ArchR)
library(argparse)
library(BSgenome.Mmusculus.UCSC.mm10)

# -------------------------------
# Arguments
# -------------------------------
parser <- ArgumentParser(description = "ArchR Motif Enrichment Analysis")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
parser$add_argument("--threads", type = "integer", default = 8, help = "Number of threads")

args <- parser$parse_args()

# -------------------------------
# Analysis
# -------------------------------
addArchRThreads(threads = args$threads)

cat("1️⃣ LOADING PROJECT...\n")
proj <- loadArchRProject(path = args$project)

output_dir <- file.path(args$project, "MotifEnrichment")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("2️⃣ ADDING MOTIF ANNOTATIONS...\n")
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

cat("3️⃣ RECOMPUTING MARKER PEAKS...\n")
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Cluster_Annotation",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 1000
)

markersPeaks_condition <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "OE",
  bgdGroups = "Control"
)

cat("4️⃣ MOTIF ENRICHMENT FOR CLUSTER MARKER PEAKS...\n")
enrichMotifs_cluster <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

cat("5️⃣ MOTIF ENRICHMENT FOR CONDITION DIFFERENTIAL PEAKS...\n")
enrichMotifs_condition <- peakAnnoEnrichment(
  seMarker = markersPeaks_condition,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

cat("6️⃣ CREATING PLOTS...\n")
heatmapEM_cluster <- plotEnrichHeatmap(enrichMotifs_cluster, n = 5, transpose = TRUE)
heatmapEM_condition <- plotEnrichHeatmap(enrichMotifs_condition, n = 10, transpose = TRUE)

cat("7️⃣ SAVING RESULTS...\n")
saveRDS(enrichMotifs_cluster, file = file.path(output_dir, "motif_enrichments_cluster.rds"))
saveRDS(enrichMotifs_condition, file = file.path(output_dir, "motif_enrichments_condition.rds"))

plotPDF(
  heatmapEM_cluster,
  heatmapEM_condition,
  name = "Motif-Enrichment-Plots",
  width = 8,
  height = 6,
  ArchRProj = proj,
  addDOC = FALSE
)

cat("✅ MOTIF ENRICHMENT ANALYSIS COMPLETE!\n")
