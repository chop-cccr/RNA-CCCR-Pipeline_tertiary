nextflow.enable.dsl=2

/*
 * This module provides a named subworkflow:
 *   RSEM_POSTREPORT_MIN( quant_ch )
 * where quant_ch: tuples (sid, genes.results, isoforms.results)
 * Emits:
 *   report_pdf -> path to rsem_report.pdf
 *   summaries  -> tuples (sid, <sid>.rsem.summary.tsv)
 */

process RSEM_SUMMARY_MIN {
  tag { sid }
  label 'r_report'
  publishDir "${params.outdir ?: 'results'}/report/summaries", mode: 'copy', overwrite: true

  input:
    tuple val(sid), path(genes), path(isoforms)

  output:
    tuple val(sid), path("${sid}.rsem.summary.tsv")

  shell:
  '''
  set -euo pipefail

  # inline R summarizer
  cat > summarize.R <<'RSCRIPT'

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(argparse)
})
p <- argparse::ArgumentParser()
p$add_argument("--sample", required=TRUE)
p$add_argument("--genes",  required=TRUE)
p$add_argument("--out",    required=TRUE)
a <- p$parse_args()

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NA) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}
sample_id <- arg_val("--sample")
genes_fp  <- arg_val("--genes")
out_fp    <- arg_val("--out")

if (is.na(sample_id) || is.na(genes_fp) || is.na(out_fp)) {
  stop("Usage: Rscript rsem_summarize_simple.R --sample <ID> --genes <genes.results> --out <output.tsv>")
}

# -------- read the RSEM genes.results (TSV) --------
genes <- read.table(
  genes_fp,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE
)

# normalize header names; fill blanks
nm <- tolower(gsub("[^a-z0-9]+", "_", names(genes)))
blank <- which(nchar(nm) == 0)
if (length(blank)) nm[blank] <- paste0("col", blank)
names(genes) <- nm

# helper to pick a column by name (case-normalized) with index fallback
pick_col <- function(df, name_candidates, fallback_idx = NA_integer_) {
  nn <- names(df)
  hit <- which(nn %in% tolower(name_candidates))
  if (length(hit) >= 1) {
    nn[hit[1]]
  } else if (!is.na(fallback_idx) && fallback_idx <= ncol(df)) {
    nn[fallback_idx]
  } else {
    NA_character_
  }
}

# RSEM canonical order: ... length [3], effective_length [4], expected_count [5], TPM [6], FPKM [7]
exp_col <- pick_col(genes, c("expected_count","expected.count"), fallback_idx = 5)
tpm_col <- pick_col(genes, c("tpm"),                            fallback_idx = 6)

if (is.na(exp_col) || is.na(tpm_col)) {
  stop("Could not locate required columns (expected_count, TPM). Have: ", paste(names(genes), collapse = ", "))
}

# numeric vectors
to_num <- function(x) suppressWarnings(as.numeric(x))
tpm_vec <- to_num(genes[[tpm_col]])
exp_vec <- to_num(genes[[exp_col]])

# -------- metrics --------
genes_detected        <- sum(tpm_vec > 0,  na.rm = TRUE)
genes_tpm_ge_0_1      <- sum(tpm_vec >= 0.1, na.rm = TRUE)
genes_tpm_ge_1        <- sum(tpm_vec >= 1,    na.rm = TRUE)
total_expected_counts <- sum(exp_vec,         na.rm = TRUE)
median_tpm            <- median(tpm_vec,      na.rm = TRUE)
q75_tpm               <- as.numeric(stats::quantile(tpm_vec, 0.75, na.rm = TRUE))

# -------- write one-row TSV --------
out_df <- data.frame(
  sample                = sample_id,
  genes_detected        = genes_detected,
  genes_tpm_ge_0_1      = genes_tpm_ge_0_1,
  genes_tpm_ge_1        = genes_tpm_ge_1,
  total_expected_counts = total_expected_counts,
  median_tpm            = median_tpm,
  q75_tpm               = q75_tpm,
  stringsAsFactors      = FALSE
)

write.table(out_df, file = out_fp, sep = "\t", quote = FALSE, row.names = FALSE)





RSCRIPT

  Rscript summarize.R \
    --sample !{sid} \
    --genes  !{genes} \
    --out    !{sid}.rsem.summary.tsv
  '''
}

process RSEM_REPORT_MIN {
  label 'r_report'
  publishDir "${params.outdir ?: 'results'}/report", mode: 'copy', overwrite: true

  input:
    path(summary_files)

  output:
    path "rsem_report.pdf"

  shell:
  '''
  set -euo pipefail
  mkdir -p _tmp

  # Merge per-sample TSVs; keep header from the first file only
  awk 'FNR==1 && NR!=1 {next} {print}' !{summary_files} > _tmp/rsem_summaries.tsv

  # inline Rmd
  cat > report.Rmd <<'RMD'
  ---
  title: "RNA-seq Quantification Report"
  output:
    pdf_document:
      toc: true
      number_sections: true
  params:
    summaries_path: "_tmp/rsem_summaries.tsv"
  ---

  ```{r setup, include=FALSE}
  suppressPackageStartupMessages({
    library(readr); library(dplyr); library(ggplot2); library(scales); library(tidyr); library(knitr)
  })
  theme_set(theme_minimal(base_size = 11))
  summ <- readr::read_tsv(params$summaries_path, show_col_types = FALSE)

summ |>
  arrange(sample) |>
  mutate(across(where(is.numeric), ~round(., 3))) |>
  knitr::kable()

ggplot(summ, aes(x = genes_detected)) +
  geom_histogram(bins = 30) +
  labs(x="Genes with TPM>0", y="Samples")

ggplot(summ, aes(x = median_tpm)) +
  geom_histogram(bins = 30) +
  labs(x="Median TPM per sample", y="Samples")


RMD

Rscript -e "rmarkdown::render('report.Rmd', output_file='rsem_report.pdf')"
#test -s rsem_report.pdf
'''

}

//
// Named subworkflow (NOT an anonymous workflow {})
//
workflow RSEM_POSTREPORT_MIN {
take:
quant_ch // tuples: (sid, genes.results, isoforms.results)

main:
summaries_ch = RSEM_SUMMARY_MIN( quant_ch )
report_pdf_ch = RSEM_REPORT_MIN( summaries_ch.map{ it[1] }.collect() )

emit:
report_pdf = report_pdf_ch
summaries = summaries_ch
}
