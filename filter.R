library(ArchR)
library(patchwork)
library(argparse)

# -------------------------------
# Parse command line arguments
# -------------------------------
parser <- ArgumentParser(description = "ArchR project")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
parser$add_argument("--minTSS", type = "double", default = 10, help = "Minimum TSS enrichment")
parser$add_argument("--minFrags", type = "double", default = 5000, help = "Minimum number of fragments")
parser$add_argument("--minGexUMI", type = "double", default = 1000, help = "Minimum GEX UMI counts")
parser$add_argument("--maxGexUMI", type = "double", default = 30000, help = "Maximum GEX UMI counts")
parser$add_argument("--minGexGenes", type = "double", default = 500, help = "Minimum number of detected genes")
parser$add_argument("--maxGexGenes", type = "double", default = 7000, help = "Maximum number of detected genes")

args <- parser$parse_args()


# -------------------------------
# Load project
# -------------------------------
addArchRThreads(threads = 4)
addArchRGenome("mm10")

project_name = args$project
proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

# Print the ArchRProject object
print(proj_ALL)

# Get a summary
summary(proj_ALL)

# Show slot names
slotNames(proj_ALL)

# Access metadata
proj_ALL@sampleMetadata
proj_ALL@cellColData

# Before filtering
cat("Before filtering:\n")
summary(proj_ALL@cellColData[, c("TSSEnrichment", "nFrags", "Gex_nUMI", "Gex_nGenes")])

cat("Number of cells per sample before filtering:\n")
print(table(proj_ALL$Sample))


# -------------------------------
# Apply filtering dynamically
# -------------------------------
proj_ALL <- proj_ALL[
  proj_ALL$TSSEnrichment > args$minTSS &
  proj_ALL$nFrags > args$minFrags &
  !is.na(proj_ALL$Gex_nUMI)
]

proj_ALL <- proj_ALL[
  proj_ALL$Gex_nGenes > args$minGexGenes &
  proj_ALL$Gex_nGenes < args$maxGexGenes &
  proj_ALL$Gex_nUMI > args$minGexUMI &
  proj_ALL$Gex_nUMI < args$maxGexUMI
]

proj_ALL <- proj_ALL[complete.cases(proj_ALL@cellColData[,
                 c("TSSEnrichment","nFrags","Gex_nUMI","Gex_nGenes")])]

# -------------------------------
# After Filtering Summary
# -------------------------------
cat("Before filtering:\n")
summary(proj_ALL@cellColData[, c("TSSEnrichment", "nFrags", "Gex_nUMI", "Gex_nGenes")])

cat("Number of cells per sample after filtering:\n")
print(table(proj_ALL$Sample))

# -------------------------------
# QC Plots after filtering PNG
# -------------------------------
figure_name_png <- paste0(project_name, "_postFilterQC.png")

png(file = figure_name_png, width = 12, height = 8, units = "in", res = 300)

p1 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p2 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p3 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "Gex_nUMI", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p4 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "Gex_nGenes", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

(p1 | p2) / (p3 | p4)

dev.off()


saveArchRProject(ArchRProj = proj_ALL, outputDirectory = project_name, load = FALSE)

