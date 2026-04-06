#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(immunedeconv)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Gene x sample TPM matrix"),
  make_option(c("-o", "--output"), type = "character", help = "Output TSV"),
  make_option(c("--sample_col"), type = "character", default = NULL,
              help = "Optional: not used now; placeholder"),
  make_option(c("--plot"), type = "character", default = NULL,
              help = "Optional heatmap PDF output")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Both --input and --output are required")
}

# Read matrix
mat <- fread(opt$input, data.table = FALSE)

# Assume first column contains HGNC gene symbols
rownames(mat) <- mat[[1]]
mat[[1]] <- NULL

# Force numeric matrix
mat <- as.matrix(mat)
mode(mat) <- "numeric"

# Clean rows
keep <- !is.na(rownames(mat)) & rownames(mat) != ""
mat <- mat[keep, , drop = FALSE]
mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]

# Run MCP-counter
res <- deconvolute(mat, method = "mcp_counter")

# Save table
fwrite(as.data.frame(res), file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

# Optional quick heatmap
if (!is.null(opt$plot)) {
  suppressPackageStartupMessages(library(pheatmap))

  df <- as.data.frame(res)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  hm <- as.matrix(df)
  mode(hm) <- "numeric"

  pdf(opt$plot, width = 10, height = 6)
  pheatmap(hm, scale = "row")
  dev.off()
}
