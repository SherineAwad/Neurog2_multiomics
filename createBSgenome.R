#!/usr/bin/env Rscript

# =============================================================================
# Script: createBSgenome.R
# Description: Create BSgenome package ONLY (no installation)
# =============================================================================

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  if (!require("BSgenome", quietly = TRUE)) {
    BiocManager::install("BSgenome")
  }

  library(BSgenome)
  library(Biostrings)
  library(rtracklayer)
})

# =============================================================================
# Configuration
# =============================================================================

fasta_file <- "neurog2.fa"
package_name <- "BSgenome.Mmusculus.custom.neurog2"

# =============================================================================
# Step 1: Clean previous attempts
# =============================================================================

cat("Step 1: Cleaning previous attempts...\n")

if(dir.exists(package_name)) unlink(package_name, recursive = TRUE)
tar_files <- list.files(pattern = paste0(package_name, ".*\\.tar\\.gz"))
if(length(tar_files) > 0) file.remove(tar_files)

# =============================================================================
# Step 2: Create package structure
# =============================================================================

cat("Step 2: Creating package structure...\n")

pkg_dir <- package_name
dir.create(pkg_dir)
dir.create(file.path(pkg_dir, "R"))
dir.create(file.path(pkg_dir, "inst"))
dir.create(file.path(pkg_dir, "inst", "extdata"))
dir.create(file.path(pkg_dir, "man"))

# =============================================================================
# Step 3: Process FASTA file and get EXACT names
# =============================================================================

cat("Step 3: Processing FASTA file...\n")

genome_sequences <- readDNAStringSet(fasta_file)
chrom_names <- names(genome_sequences)

# Use EXACT names from FASTA - no modifications
cat("Using EXACT sequence names from FASTA:\n")
print(chrom_names[1:5])  # Show first 5 to verify
cat("Total sequences:", length(chrom_names), "\n")

# Create 2bit file
two_bit_file <- file.path(pkg_dir, "inst", "extdata", "single_sequences.2bit")
export(genome_sequences, two_bit_file, format = "2bit")
cat("Created 2bit file\n")

# =============================================================================
# Step 4: Create DESCRIPTION file
# =============================================================================

cat("Step 4: Creating DESCRIPTION file...\n")

desc_content <- sprintf('Package: %s
Title: Full genome sequences for Mus musculus (neurog2)
Description: Full genome sequences for Mus musculus (neurog2)
Version: 1.0.0
Author: Custom
Maintainer: Custom <custom@example.com>
License: Artistic-2.0
Depends: R (>= 4.0.0), BSgenome
Imports: BSgenome
Suggests: BiocStyle
organism: Mus musculus
common_name: Mouse
genome: neurog2
provider: Custom
provider_version: neurog2
release_date: %s
release_name: neurog2
BSgenomeObjname: Mmusculus', package_name, format(Sys.Date(), "%Y/%m/%d"))

writeLines(desc_content, file.path(pkg_dir, "DESCRIPTION"))

# =============================================================================
# Step 5: Create NAMESPACE file
# =============================================================================

cat("Step 5: Creating NAMESPACE file...\n")

namespace_content <- 'export("Mmusculus")
export("BSgenome.Mmusculus.custom.neurog2")'

writeLines(namespace_content, file.path(pkg_dir, "NAMESPACE"))

# =============================================================================
# Step 6: Create R file with EXACT sequence names
# =============================================================================

cat("Step 6: Creating R file...\n")

# Create the seqnames vector manually to ensure exact match
seqnames_vector <- paste0('c("', paste(chrom_names, collapse = '", "'), '")')

r_code <- sprintf('
.Mmusculus <- BSgenome(
    organism = "Mus musculus",
    common_name = "Mouse",
    genome = "neurog2",
    provider = "Custom", 
    provider_version = "neurog2",
    release_date = "%s",
    release_name = "neurog2",
    source_url = "",
    seqnames = %s,
    circ_seqs = character(0),
    mseqnames = NULL,
    seqs_pkgname = "%s",
    seqs_dirpath = system.file("extdata", package = "%s")
)

# Set the pkgname slot
.Mmusculus@pkgname <- "BSgenome.Mmusculus.custom.neurog2"

# Set COMPLETE metadata
metadata(.Mmusculus) <- list(
    genome = "neurog2",
    organism = "Mus musculus",
    provider = "Custom",
    provider_version = "neurog2", 
    release_date = "%s",
    release_name = "neurog2",
    source_url = "",
    common_name = "Mouse",
    BSgenomeObjname = "Mmusculus"
)

Mmusculus <- .Mmusculus
BSgenome.Mmusculus.custom.neurog2 <- .Mmusculus
',
format(Sys.Date(), "%Y/%m/%d"),
seqnames_vector,
package_name,
package_name,
format(Sys.Date(), "%Y/%m/%d"))

writeLines(r_code, file.path(pkg_dir, "R", paste0(package_name, ".R")))

# =============================================================================
# Step 7: Build the package only (no install)
# =============================================================================

cat("Step 7: Building package...\n")

build_cmd <- paste("R CMD build", pkg_dir)
system(build_cmd)

# Check if tar file was created
tar_file <- list.files(pattern = paste0(package_name, ".*\\.tar\\.gz"))[1]
if(file.exists(tar_file)) {
  cat("✅ SUCCESS: Package built as", tar_file, "\n")
  cat("To install, run: R CMD INSTALL --no-staged-install", tar_file, "\n")
} else {
  cat("❌ FAILED: Package tar file not created\n")
}

cat("\n=== BSGENOME PACKAGE CREATION COMPLETE ===\n")
