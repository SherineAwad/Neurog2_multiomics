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

# Parse arguments first
parser <- ArgumentParser(description = "ArchR project")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
args <- parser$parse_args()
project_name = args$project

# Define files FIRST - before they're used
atacFiles <- c("Control" = "TH1_atac_fragments.tsv.gz", "OE" = "TH2_atac_fragments.tsv.gz")
rnaFiles <- c("Control" = "TH1_filtered_feature_bc_matrix.h5", "OE" = "TH2_filtered_feature_bc_matrix.h5")
all.equal(names(atacFiles), names(rnaFiles))

# Create gene annotation from your neurog2.gtf file FIRST
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

# Load mm10 genome to create proper genome annotation
cat("Loading mm10 genome for annotation structure...\n")
library(BSgenome.Mmusculus.UCSC.mm10)
mm10_genome <- BSgenome.Mmusculus.UCSC.mm10

# Create proper mm10 genome annotation
mm10_genomeAnnotation <- SimpleList(
  genome = "mm10",
  chromSizes = GRanges(
    seqnames = names(seqlengths(mm10_genome)),
    ranges = IRanges(start = 1, end = seqlengths(mm10_genome)),
    seqinfo = seqinfo(mm10_genome)
  ),
  blacklist = GRanges()  # Empty blacklist
)

cat("Using mm10 genome for initial Arrow file creation...\n")
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = mm10_genomeAnnotation,  # Use the mm10 genome annotation we created
  minTSS = 4,
  minFrags = 500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = FALSE
)

# Create project with mm10
cat("Creating ArchR project with mm10...\n")
proj_ALL <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = project_name,
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = mm10_genomeAnnotation
)

# NOW load your custom BSgenome and swap it in
cat("Loading custom BSgenome and swapping into project...\n")
library(BSgenome.Mmusculus.custom.neurog2)
custom_genome <- Mmusculus

# Create custom genome annotation
custom_genomeAnnotation <- SimpleList(
  genome = "BSgenome.Mmusculus.custom.neurog2",
  chromSizes = GRanges(
    seqnames = names(seqlengths(custom_genome)),
    ranges = IRanges(start = 1, end = seqlengths(custom_genome)),
    seqinfo = seqinfo(custom_genome)
  ),
  blacklist = GRanges()
)

# Replace the genome annotation in the project
proj_ALL@genomeAnnotation <- custom_genomeAnnotation

# Verify the swap worked
cat("Project genome after swap:", proj_ALL@genomeAnnotation$genome, "\n")
cat("Custom genome organism:", organism(custom_genome), "\n")

# Continue with the rest of your analysis...
seRNA_A <- import10xFeatureMatrix(input = "TH1_filtered_feature_bc_matrix.h5", names = "Control")
seRNA_B <- import10xFeatureMatrix(input = "TH2_filtered_feature_bc_matrix.h5", names = "OE")

cat("\n--- RNA data summary ---\n")
cat("Total genes in seRNA_A:", nrow(seRNA_A), "\n")
cat("Total genes in seRNA_B:", nrow(seRNA_B), "\n")

cat("\n=== SKIPPING GENE FILTERING ===\n")

# Keep ALL genes - no filtering
seRNA_A_filtered <- seRNA_A
seRNA_B_filtered <- seRNA_B

# Filter cells to match those in the ArchR project
k_cells_control <- which(rownames(proj_ALL@cellColData[proj_ALL$Sample == "Control",]) %in% colnames(seRNA_A_filtered))
k_cells_oe <- which(rownames(proj_ALL@cellColData[proj_ALL$Sample == "OE",]) %in% colnames(seRNA_B_filtered))

# Subset the main project
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

cat("🎉 ArchR analysis completed successfully with custom genome!\n")
