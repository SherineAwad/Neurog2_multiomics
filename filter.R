library(ArchR)
library(patchwork)
library(argparse)

# -------------------------------
# Parse command line arguments
# -------------------------------
parser <- ArgumentParser(description = "Filter ArchR project cells and generate post-QC plots")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
parser$add_argument("--minTSS", type = "double", default = 4)
parser$add_argument("--minFrags", type = "double", default = 2000)
parser$add_argument("--minGexUMI", type = "double", default = 10)
parser$add_argument("--maxGexUMI", type = "double", default = 30000)
parser$add_argument("--minGexGenes", type = "double", default = 10)
parser$add_argument("--maxGexGenes", type = "double", default = 7000)
args <- parser$parse_args()

# -------------------------------
# Load project
# -------------------------------
addArchRThreads(threads = 4)
addArchRGenome("mm10")

proj_ALL <- loadArchRProject(path = args$project, force = FALSE, showLogo = TRUE)

# Print the ArchRProject object
print(proj_ALL)

# Get a summary
summary(proj_ALL)

# Show slot names
slotNames(proj_ALL)

# Access metadata
proj_ALL@sampleMetadata
proj_ALL@cellColData

valid_cells <- proj_ALL$cellNames[complete.cases(proj_ALL@cellColData[, c("TSSEnrichment","nFrags","Gex_nUMI","Gex_nGenes")])]
proj_ALL <- subsetArchRProject(proj_ALL, cells = valid_cells)


# -------------------------------
# Before Filtering Summary
# -------------------------------
cat("Number of cells per sample before filtering:\n")
print(table(proj_ALL$Sample))

figure_name_pdf <- paste0(args$project, "_beforeFilterQC.pdf")
pdf(file = figure_name_pdf, width = 12, height = 8)
p1 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "TSSEnrichment", plotAs = "violin", alpha = 0.4,
                 addBoxPlot = TRUE, na.rm = TRUE)
p2 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "log10(nFrags)", plotAs = "violin", alpha = 0.4,
                 addBoxPlot = TRUE, na.rm = TRUE)
p3 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "Gex_nUMI", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE,na.rm = TRUE)
p4 <- plotGroups(proj_ALL, groupBy = "Sample", colorBy = "cellColData",
                 name = "Gex_nGenes", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE,na.rm = TRUE)
(p1 | p2) / (p3 | p4)
dev.off()
if (FALSE)
{
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
cat("Number of cells per sample after filtering:\n")
print(table(proj_ALL$Sample))

# -------------------------------
# QC Plots after filtering (PDF + PNG)
# -------------------------------
figure_name_pdf <- paste0(args$project, "_postFilterQC.pdf")
pdf(file = figure_name_pdf, width = 12, height = 8)

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
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "Neurog2", load = FALSE)
}

