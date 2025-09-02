library(ArchR)
library(patchwork)
library(argparse)

# -------------------------------
# Parse command line arguments
# -------------------------------
parser <- ArgumentParser(description = "Filter ArchR project cells and generate post-QC plots")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
addArchRThreads(threads = 4)
addArchRGenome("mm10")

args <- parser$parse_args()

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

# Logical vectors
has_rna <- !is.na(proj_ALL$Gex_nUMI) & !is.na(proj_ALL$Gex_nGenes)
has_atac <- !is.na(proj_ALL$nFrags)  # basically all cells in ArchRProject have ATAC

# Make categories
cell_type <- ifelse(has_rna & has_atac, "ATAC+RNA",
                    ifelse(has_atac & !has_rna, "ATAC_only", "Other"))

# Create a table
table(cell_type)


# -------------------------------
# Before Filtering Summary
# -------------------------------
cat("Number of cells per sample before filtering:\n")
print(table(proj_ALL$Sample))

# Identify cells with RNA
rna_cells <- proj_ALL$cellNames[!is.na(proj_ALL$Gex_nUMI) & !is.na(proj_ALL$Gex_nGenes)]

# Now plot only these cells for RNA metrics
p3 <- plotGroups(proj_ALL[rna_cells, ], groupBy="Sample", colorBy="cellColData",
                 name="Gex_nUMI", plotAs="violin", alpha=0.4, addBoxPlot=TRUE)
p4 <- plotGroups(proj_ALL[rna_cells, ], groupBy="Sample", colorBy="cellColData",
                 name="Gex_nGenes", plotAs="violin", alpha=0.4, addBoxPlot=TRUE)

# ATAC metrics: can still use all cells
p1 <- plotGroups(proj_ALL, groupBy="Sample", colorBy="cellColData",
                 name="TSSEnrichment", plotAs="violin", alpha=0.4, addBoxPlot=TRUE, na.rm=TRUE)
p2 <- plotGroups(proj_ALL, groupBy="Sample", colorBy="cellColData",
                 name="log10(nFrags)", plotAs="violin", alpha=0.4, addBoxPlot=TRUE, na.rm=TRUE)

# Save PDF
figure_name_pdf <- paste0(args$project, "_beforeFilterQC.pdf")
pdf(file=figure_name_pdf, width=12, height=8)
(p1 | p2) / (p3 | p4)
dev.off()





