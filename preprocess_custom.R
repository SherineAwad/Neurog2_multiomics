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

# Parse arguments
parser <- ArgumentParser(description = "ArchR project")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
args <- parser$parse_args()
project_name = args$project

# Define files
atacFiles <- c("Control" = "TH1_atac_fragments.tsv.gz", "OE" = "TH2_atac_fragments.tsv.gz")
rnaFiles <- c("Control" = "TH1_filtered_feature_bc_matrix.h5", "OE" = "TH2_filtered_feature_bc_matrix.h5")

# Load custom BSgenome
library(BSgenome.Mmusculus.custom.neurog2)
custom_genome <- Mmusculus

# Create gene annotation
gtf_file <- "neurog2.gtf"
gtf <- import(gtf_file)
genes <- gtf[gtf$type == "gene"]
exons <- gtf[gtf$type == "exon"]
TSS <- resize(genes, width = 1, fix = "start")

geneAnnotation <- createGeneAnnotation(
  genes = genes,
  exons = exons,
  TSS = TSS
)

# Create genome annotation
genomeAnnotation <- createGenomeAnnotation(genome = custom_genome, filter=FALSE, filterChr=NULL)
# Create ArrowFiles
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 4,
  minFrags = 500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

# Create project
proj_ALL <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = project_name,
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

# Import RNA data
seRNA_A <- import10xFeatureMatrix(input = "TH1_filtered_feature_bc_matrix.h5", names = "Control")
seRNA_B <- import10xFeatureMatrix(input = "TH2_filtered_feature_bc_matrix.h5", names = "OE")

# Keep all genes
seRNA_A_filtered <- seRNA_A
seRNA_B_filtered <- seRNA_B

# Filter cells
k_cells_control <- which(rownames(proj_ALL@cellColData[proj_ALL$Sample == "Control",]) %in% colnames(seRNA_A_filtered))
k_cells_oe <- which(rownames(proj_ALL@cellColData[proj_ALL$Sample == "OE",]) %in% colnames(seRNA_B_filtered))

cells_to_keep <- c(rownames(proj_ALL@cellColData[proj_ALL$Sample == "Control",])[k_cells_control],
                   rownames(proj_ALL@cellColData[proj_ALL$Sample == "OE",])[k_cells_oe])
proj_ALL <- proj_ALL[cells_to_keep, ]

# Combine RNA data
seRNAcombined <- cbind(assay(seRNA_A_filtered), assay(seRNA_B_filtered))
seRNA_all <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA_A_filtered))

# Add to ArchR project
proj_ALL <- addGeneExpressionMatrix(input = proj_ALL, seRNA = seRNA_all)

# Save project
saveArchRProject(ArchRProj = proj_ALL, outputDirectory = project_name, load = FALSE)

cat("Done\n")
