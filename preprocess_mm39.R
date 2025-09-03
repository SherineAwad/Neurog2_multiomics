library(ArchR)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(pheatmap)
library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.mm39)      # Genome sequences (chrom sizes)
library(ensembldb)


setwd("/nfs/turbo/umms-thahoang/sherine/Neurog2_multiomics")
addArchRThreads(threads = 6)
project_name ="Neurog2_mm39"
atacFiles <- c("ctrl" = "TH1_atac_fragments.tsv.gz", "ko" = "TH2_atac_fragments.tsv.gz")
rnaFiles <- c("ctrl" = "TH1_filtered_feature_bc_matrix.h5", "ko" = "TH2_filtered_feature_bc_matrix.h5")
all.equal(names(atacFiles), names(rnaFiles))


# Create a gene annotation object from EnsDb
edb <- makeEnsembldbFromGtf(
  gtf = "Mus_musculus.GRCm39.103.gtf",
  genome = "Mus_musculus.GRCm39.dna.primary_assembly.fa"
)

# Create a gene annotation object from EnsDb
geneAnnotation <- createGeneAnnotation(
  edb = edb,
  genome = BSgenome.Mmusculus.UCSC.mm39
)

addGeneAnnotation(geneAnnotation)

ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4,
  minFrags = 1000,
  geneAnnotation = geneAnnotation,   # <--- use your EnsDb annotation
  addGeneScoreMat = TRUE,
  force = FALSE
)


ArrowFiles <- c("NControl.arrow","NKO.arrow")
project_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "Neurog2_mm39", copyArrows = FALSE)

proj_A <- ArchRProject(ArrowFiles[1], outputDirectory = "Control", copyArrows = FALSE)
proj_B <- ArchRProject(ArrowFiles[2], outputDirectory = "KO", copyArrows = FALSE)

seRNA_A <- import10xFeatureMatrix(input="TH1_filtered_feature_bc_matrix.h5", names="Control")
seRNA_B <- import10xFeatureMatrix(input="TH2_filtered_feature_bc_matrix.h5", names="KO")


k <- which(seqnames(rowRanges(seRNA_A)) %in% seqnames(proj_A@genomeAnnotation$chromSizes) == T)
seRNA_A = seRNA_A[k,]
k = which(rownames(proj_A@cellColData) %in% colnames(seRNA_A) == T)
proj_A = proj_A[k]

k <- which(seqnames(rowRanges(seRNA_B)) %in% seqnames(proj_B@genomeAnnotation$chromSizes) == T)
seRNA_B = seRNA_B[k,]
k = which(rownames(proj_B@cellColData) %in% colnames(seRNA_B) == T)
proj_B = proj_B[k]
seRNAcombined <- cbind(assay(seRNA_A), assay(seRNA_B))

seRNA_all <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA_A))
proj_ALL <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "Neurog2_mm39", copyArrows = TRUE)
proj_ALL <- addGeneExpressionMatrix(input = project_ALL, seRNA = seRNA_all)


# 1️⃣ Clean ATAC barcodes from proj_ALL
atac_clean <- gsub(".*#", "", rownames(proj_ALL@cellColData))
cat("Number of ATAC cells:", length(atac_clean), "\n")

# 2️⃣ Clean RNA barcodes from the combined SE
rna_clean <- gsub("-1$", "", colnames(seRNAcombined))
cat("Number of RNA barcodes:", length(rna_clean), "\n")

# 3️⃣ Optional: count how many are in both
common_cells <- intersect(atac_clean, rna_clean)
cat("Number of barcodes present in both ATAC and RNA:", length(common_cells), "\n")

saveArchRProject(ArchRProj = proj_ALL, outputDirectory = "Neurog2_mm39", load = FALSE)
