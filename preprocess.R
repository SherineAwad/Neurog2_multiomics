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
library(GenomicRanges)
library(IRanges)

setwd(getwd())
addArchRThreads(threads = 6)
addArchRGenome("mm10")
parser <- ArgumentParser(description = "ArchR project")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
args <- parser$parse_args()
project_name = args$project

atacFiles <- c("Control" = "TH1_atac_fragments.tsv.gz", "OE" = "TH2_atac_fragments.tsv.gz")
rnaFiles <- c("Control" = "TH1_filtered_feature_bc_matrix.h5", "OE" = "TH2_filtered_feature_bc_matrix.h5")
all.equal(names(atacFiles), names(rnaFiles))

ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4,
  minFrags = 500,  addTileMat = TRUE,
  addGeneScoreMat = TRUE, force=FALSE
)

ArrowFiles <- c("Control.arrow","OE.arrow")
project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = project_name, copyArrows = FALSE)

proj_A <- ArchRProject(ArrowFiles[1], outputDirectory = "Control", copyArrows = FALSE)
proj_B <- ArchRProject(ArrowFiles[2], outputDirectory = "OE", copyArrows = FALSE)

seRNA_A <- import10xFeatureMatrix(input="TH1_filtered_feature_bc_matrix.h5", names="Control")
seRNA_B <- import10xFeatureMatrix(input="TH2_filtered_feature_bc_matrix.h5", names="OE")

# Genes to keep
genes_to_keep <- c("Neurog2_S9A", "mScarlet")

# ---------------- Print all genes before filtering ----------------
writeLines(rownames(seRNA_A), paste0(project_name, "_seRNA_A_before_filtering.txt"))
writeLines(rownames(seRNA_B), paste0(project_name, "_seRNA_B_before_filtering.txt"))

cat("\n--- Before any processing seRNA_A ---\n")
cat("Total genes in seRNA_A:", nrow(seRNA_A), "\n")

# ---------------- SKIP GENE FILTERING ----------------
cat("\n=== SKIPPING GENE FILTERING ===\n")

# Keep ALL genes - no filtering
seRNA_A_filtered <- seRNA_A
seRNA_B_filtered <- seRNA_B

# Only filter cells (this is necessary)
k_cells_A <- which(rownames(proj_A@cellColData) %in% colnames(seRNA_A_filtered) == TRUE)
proj_A <- proj_A[k_cells_A]

k_cells_B <- which(rownames(proj_B@cellColData) %in% colnames(seRNA_B_filtered) == TRUE)
proj_B <- proj_B[k_cells_B]

# ---------------- Combine RNA ----------------
seRNAcombined <- cbind(assay(seRNA_A_filtered), assay(seRNA_B_filtered))
seRNA_all <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA_A_filtered))

# ---------------- Print all genes after combining ----------------
writeLines(rownames(seRNA_all), paste0(project_name, "_seRNA_all_after_combining.txt"))

# ---------------- FIX GENOMIC COORDINATES ----------------
cat("\n=== FIXING GENOMIC COORDINATES ===\n")

# Get valid chromosomes from the EXISTING ArchR project
valid_chroms <- seqnames(project_ALL@genomeAnnotation$chromSizes)

# Get the current rowRanges
current_ranges <- rowRanges(seRNA_all)

# Count how many genes have invalid chromosomes
all_chroms <- as.character(seqnames(current_ranges))
invalid_chroms <- !all_chroms %in% valid_chroms
cat("Genes with invalid chromosomes:", sum(invalid_chroms), "/", length(all_chroms), "\n")

# FAST APPROACH: Use vectorized operations
new_seqnames <- as.character(seqnames(current_ranges))
new_starts <- start(current_ranges)
new_ends <- end(current_ranges)
new_strands <- as.character(strand(current_ranges))

# Fix invalid chromosomes
for(i in which(invalid_chroms)) {
  gene_name <- names(current_ranges)[i]
  
  if(gene_name == "Neurog2_S9A") {
    new_seqnames[i] <- "chr2"
    new_starts[i] <- 157595000
    new_ends[i] <- 157596000
    new_strands[i] <- "+"
  } else if(gene_name == "mScarlet") {
    new_seqnames[i] <- "chr6"
    new_starts[i] <- 113066000
    new_ends[i] <- 113067000
    new_strands[i] <- "+"
  } else {
    new_seqnames[i] <- "chr19"
    new_starts[i] <- 10000000 + (i * 1000)
    new_ends[i] <- new_starts[i] + 1000
    new_strands[i] <- "+"
  }
}

# Create new GRanges object all at once (fast)
new_ranges <- GRanges(
  seqnames = new_seqnames,
  ranges = IRanges(start = new_starts, end = new_ends),
  strand = new_strands
)
names(new_ranges) <- names(current_ranges)

# Update the rowRanges
rowRanges(seRNA_all) <- new_ranges

cat("✅ Successfully fixed genomic coordinates for", sum(invalid_chroms), "genes\n")

# ---------------- Add to ArchR project ----------------
proj_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = project_name, copyArrows = TRUE)
proj_ALL <- addGeneExpressionMatrix(input = proj_ALL, seRNA = seRNA_all)

# Final verification - FIXED FUNCTION NAME
cat("\n=== FINAL VERIFICATION ===\n")
gene_scores <- getMatrixFromProject(proj_ALL, "GeneScoreMatrix")
custom_in_final <- genes_to_keep[genes_to_keep %in% rownames(gene_scores)]
cat("Custom genes in final ArchR project:", if(length(custom_in_final)==0) "None" else paste(custom_in_final, collapse = ", "), "\n")

if(length(custom_in_final) > 0) {
  cat("✅ SUCCESS: Custom genes preserved in final project!\n")
} else {
  cat("❌ FAILURE: Custom genes not found in final project\n")
}

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = project_name, load = FALSE)
