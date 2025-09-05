# -------------------------------
# Load required packages
# -------------------------------
library(ArchR)
library(argparse)
library(ggplot2)

# -------------------------------
# Parse command line arguments
# -------------------------------
parser <- ArgumentParser(description = "ArchR project")
parser$add_argument("--project", required = TRUE, help = "Path to ArchR project")
parser$add_argument("--markers", required = TRUE, help = "List of marker genes")
args <- parser$parse_args()

project_name <- args$project
markers_file <- args$markers

# -------------------------------
# Set threads and genome
# -------------------------------
addArchRThreads(threads = 4)
addArchRGenome("mm10")

# -------------------------------
# Load ArchR Project
# -------------------------------
proj_ALL <- loadArchRProject(path = project_name, force = FALSE, showLogo = TRUE)

# -------------------------------
# Read marker genes
# -------------------------------
markerGenes <- readLines(markers_file)
markerGenes <- trimws(markerGenes)

# -------------------------------
# Define plotting function with visible UMAP structure and axis labels
# -------------------------------
plot_gene_umap_zscore <- function(proj, gene_symbol) {
  
  message("Plotting gene: ", gene_symbol)
  flush.console()
  
  # Get expression matrix
  expr_matrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
  log_norm_expr <- assay(expr_matrix)
  gene_names <- rowData(expr_matrix)$name
  if(is.null(gene_names)) {
    warning(paste("GeneExpressionMatrix has no rownames. Skipping gene", gene_symbol))
    return(NULL)
  }
  rownames(log_norm_expr) <- gene_names
  
  # Find gene index case-insensitive
  gene_idx <- which(tolower(rownames(log_norm_expr)) == tolower(gene_symbol))
  if(length(gene_idx) == 0) {
    warning(paste("Gene", gene_symbol, "not found in expression matrix. Skipping."))
    return(NULL)
  }
  
  expr_values <- log_norm_expr[gene_idx, ]
  
  # Get UMAP embedding
  embedding <- getEmbedding(proj, embedding = "UMAP_Combined", returnDF = TRUE)
  if(ncol(embedding) == 2) {
    colnames(embedding) <- c("UMAP_1", "UMAP_2")
  } else {
    stop("Embedding does not have exactly 2 columns! Please check embedding name or data.")
  }
  
  # Match cells
  common_cells <- intersect(rownames(embedding), colnames(log_norm_expr))
  if(length(common_cells) == 0) stop("No common cells between embedding and expression matrix!")
  
  embedding_sub <- embedding[common_cells, , drop = FALSE]
  expr_sub <- expr_values[common_cells]
  
  # Z-score normalization
  expr_zscore <- (expr_sub - mean(expr_sub)) / sd(expr_sub)
  
  # Data frame for plotting
  df <- data.frame(embedding_sub, ExpressionZ = expr_zscore)
  
  # Split expressed vs non-expressed for coloring
  df$ExpressionCategory <- ifelse(df$ExpressionZ > 0, "expressed", "not_expressed")
  
  # Create ggplot
  p <- ggplot() +
    # All points in light grey
    geom_point(data = df[df$ExpressionCategory == "not_expressed", ],
               aes(x = UMAP_1, y = UMAP_2),
               color = "lightgrey", size = 0.5) +
    # Expressed points in gradient from light red to deep red
    geom_point(data = df[df$ExpressionCategory == "expressed", ],
               aes(x = UMAP_1, y = UMAP_2, color = ExpressionZ),
               size = 0.5) +
    scale_color_gradient(low = "#FFA07A", high = "#8B0000") +  # light salmon to deep red
    theme_classic() +
    ggtitle(paste0(project_name, "_UMAP_Zscore_", gene_symbol)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

# -------------------------------
# Loop through marker genes and save PNG
# -------------------------------
for (gene in markerGenes) {
  p <- plot_gene_umap_zscore(proj_ALL, gene)
  
  # Skip if plot returned NULL (gene not found)
  if (is.null(p)) next
  
  # Save PNG
  figure_name <- paste0(project_name, "_UMAP_Zscore_", gene, ".png")
  ggsave(
    filename = figure_name,
    plot = p,
    width = 5,
    height = 5,
    units = "in",
    dpi = 300
  )
}

