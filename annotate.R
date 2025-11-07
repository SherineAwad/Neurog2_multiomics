#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(argparse)
  library(readr)
  library(ggplot2)
})

# -------------------------------
# Arguments
# -------------------------------
parser <- ArgumentParser(description = "Annotate ArchR clusters, remove specific clusters, and plot UMAP")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
parser$add_argument("--annotations", required = TRUE, help = "CSV with Cluster,Annotation columns")
args <- parser$parse_args()

proj_path <- args$project
annotation_file <- args$annotations

# -------------------------------
# Load ArchR project
# -------------------------------
cat("Loading ArchR project from:", proj_path, "\n")
proj <- loadArchRProject(path = proj_path, force = FALSE)

# -------------------------------
# Remove unwanted clusters (in-memory)
# -------------------------------
clusters_to_remove <- c("C1", "C2", "C3", "C4", "C17", "C18", "C21", "C22", "C23", "C24")
cat("Removing clusters:", paste(clusters_to_remove, collapse = ", "), "\n")

if(!"Clusters_Combined" %in% colnames(proj@cellColData)){
  stop("Column 'Clusters_Combined' not found in project metadata.")
}

# Filter the metadata and cell list directly
keep_cells <- proj$cellNames[!proj$Clusters_Combined %in% clusters_to_remove]
proj <- proj[keep_cells, ]  # <-- in-memory filtering

cat(length(keep_cells), "cells kept after filtering\n")

# -------------------------------
# Load cluster annotations (fixed)
# -------------------------------
cat("Loading annotations from:", annotation_file, "\n")
annotations <- read_csv(annotation_file, col_types = cols(
  Cluster = col_character(),      # treat cluster IDs as character
  Annotation = col_character()
))

# Convert project clusters to character
proj_clusters <- as.character(proj$Clusters_Combined)

# Create mapping
cluster_annotation_map <- setNames(annotations$Annotation, annotations$Cluster)

# Assign annotations to cells
proj$Cluster_Annotation <- cluster_annotation_map[proj_clusters]

# Any missing clusters -> Unknown
proj$Cluster_Annotation[is.na(proj$Cluster_Annotation)] <- "Unknown"

# -------------------------------
# Plot UMAP
# -------------------------------
cat("Generating UMAP plot...\n")
umap <- getEmbedding(proj, embedding = "UMAP_Combined")
umap_df <- data.frame(
  UMAP1 = umap[,1],
  UMAP2 = umap[,2],
  Annotation = proj$Cluster_Annotation
)

pdf(file = file.path(proj_path, "UMAP_Cluster_Annotation.pdf"), width = 6, height = 5)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Annotation)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "UMAP of ArchR Clusters", color = "Annotation")
dev.off()

# -------------------------------
# Save back to same project
# -------------------------------
cat("Saving project back to:", proj_path, "\n")
saveArchRProject(proj, outputDirectory = proj_path, load = FALSE)

cat("✅ Project saved successfully — annotations included, plot generated.\n")

