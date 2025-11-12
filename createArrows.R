library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(pheatmap)
library(chromVARmotifs)
library(argparse)
library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(rtracklayer)

setwd(getwd())
addArchRThreads(threads = 6)

# Load custom BSgenome properly
cat("Loading custom BSgenome...\n")
library(BSgenome.Mmusculus.custom.neurog2)

# Get the actual genome object - it's called Mmusculus per your DESCRIPTION
cat("Available objects in package:\n")
print(ls("package:BSgenome.Mmusculus.custom.neurog2"))

# Use the actual genome object name
custom_genome <- Mmusculus
cat("Genome object class:", class(custom_genome), "\n")
cat("Genome organism:", organism(custom_genome), "\n")

# Create genome annotation using the ACTUAL genome object
cat("Creating custom genome annotation...\n")
genomeAnnotation <- createGenomeAnnotation(
  genome = custom_genome,  # Use the object, not the package name
  blacklist = NULL,
  filter = FALSE
)

# Create gene annotation from your neurog2.gtf file
cat("Creating gene annotation from neurog2.gtf...\n")
gtf_file <- "neurog2.gtf"
if (file.exists(gtf_file)) {
  gtf <- import(gtf_file)
  genes <- gtf[gtf$type == "gene"]
  exons <- gtf[gtf$type == "exon"]
  TSS <- resize(genes, width = 1, fix = "start")
  
  geneAnnotation <- createGeneAnnotation(
    genes = genes,
    exons = exons,
    TSS = TSS
  )
  cat("✅ Gene annotation created successfully from", gtf_file, "\n")
} else {
  stop("GTF file not found: ", gtf_file)
}

parser <- ArgumentParser(description = "ArchR project")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
args <- parser$parse_args()
project_name = args$project

atacFiles <- c("Control" = "TH1_atac_fragments.tsv.gz", "OE" = "TH2_atac_fragments.tsv.gz")
rnaFiles <- c("Control" = "TH1_filtered_feature_bc_matrix.h5", "OE" = "TH2_filtered_feature_bc_matrix.h5")
all.equal(names(atacFiles), names(rnaFiles))

cat("Creating ArrowFiles with custom genome annotations...\n")
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 4,
  minFrags = 500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = FALSE
)

cat("Creating ArchR project...\n")
proj_ALL <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = project_name,
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

seRNA_A <- import10xFeatureMatrix(input = "TH1_filtered_feature_bc_matrix.h5", names = "Control")
seRNA_B <- import10xFeatureMatrix(input = "TH2_filtered_feature_bc_matrix.h5", names = "OE")

cat("\n--- Before any processing seRNA_A ---\n")
cat("Total genes in seRNA_A:", nrow(seRNA_A), "\n")

cat("\n=== SKIPPING GENE FILTERING ===\n")

# Keep ALL genes - no filtering
seRNA_A_filtered <- seRNA_A
seRNA_B_filtered <- seRNA_B

# Only filter cells (this is necessary) - now using proj_ALL
k_cells_control <- which(rownames(proj_ALL@cellColData[proj_ALL$Sample == "Control",]) %in% colnames(seRNA_A_filtered) == TRUE)
k_cells_oe <- which(rownames(proj_ALL@cellColData[proj_ALL$Sample == "OE",]) %in% colnames(seRNA_B_filtered) == TRUE)

# Subset the main project
cells_to_keep <- c(rownames(proj_ALL@cellColData[proj_ALL$Sample == "Control",])[k_cells_control],
                   rownames(proj_ALL@cellColData[proj_ALL$Sample == "OE",])[k_cells_oe])
proj_ALL <- proj_ALL[cells_to_keep, ]

# Combine RNA
seRNAcombined <- cbind(assay(seRNA_A_filtered), assay(seRNA_B_filtered))
seRNA_all <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA_A_filtered))

# Add to ArchR project
proj_ALL <- addGeneExpressionMatrix(input = proj_ALL, seRNA = seRNA_all)

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = project_name, load = FALSE)

cat("🎉 ArchR analysis completed successfully with custom genome!\n")
