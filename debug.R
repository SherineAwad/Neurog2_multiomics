suppressPackageStartupMessages({
  library(ArchR)
  library(argparse)
})

# -------------------------------
# Arguments
# -------------------------------
parser <- ArgumentParser(description = "List embeddings and clusters in ArchR project")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
args <- parser$parse_args()

proj_path <- args$project

# -------------------------------
# Load ArchR project
# -------------------------------
cat("Loading ArchR project from:", proj_path, "\n")
proj <- loadArchRProject(path = proj_path, force = FALSE)

# -------------------------------
# List embeddings
# -------------------------------
embeddings <- names(proj@reducedDims)
cat("\n=== AVAILABLE EMBEDDINGS ===\n")
if (length(embeddings) == 0) {
  cat("❌ No embeddings found.\n")
} else {
  for (emb in embeddings) {
    cat(" -", emb, "\n")
  }
}

# -------------------------------
# List all cluster columns
# -------------------------------
cat("\n=== AVAILABLE CLUSTER COLUMNS ===\n")
all_columns <- names(getCellColData(proj))
cluster_columns <- grep("cluster|Cluster", all_columns, value = TRUE, ignore.case = TRUE)

if (length(cluster_columns) == 0) {
  cat("❌ No cluster columns found.\n")
} else {
  for (clust_col in cluster_columns) {
    cat("📊 Column:", clust_col, "\n")
    clusters <- table(getCellColData(proj)[[clust_col]])
    cat("   Clusters:", paste(names(clusters), collapse = ", "), "\n")
    cat("   Counts:", paste(clusters, collapse = ", "), "\n")
    cat("   Total cells:", sum(clusters), "\n\n")
  }
}

# -------------------------------
# List all sample/condition columns
# -------------------------------
cat("\n=== AVAILABLE SAMPLE/CONDITION COLUMNS ===\n")
sample_columns <- grep("sample|Sample|condition|Condition|group|Group", all_columns, value = TRUE, ignore.case = TRUE)

if (length(sample_columns) == 0) {
  cat("❌ No sample/condition columns found.\n")
} else {
  for (sample_col in sample_columns) {
    cat("🔬 Column:", sample_col, "\n")
    samples <- table(getCellColData(proj)[[sample_col]])
    cat("   Values:", paste(names(samples), collapse = ", "), "\n")
    cat("   Counts:", paste(samples, collapse = ", "), "\n\n")
  }
}

# -------------------------------
# Available matrices
# -------------------------------
cat("\n=== AVAILABLE MATRICES ===\n")
matrices <- getAvailableMatrices(proj)
if (length(matrices) == 0) {
  cat("❌ No matrices found.\n")
} else {
  for (mat in matrices) {
    cat(" -", mat, "\n")
  }
}

cat("\n=== DEBUG COMPLETE ===\n")
