nextflow.enable.dsl=2

/*
 * RSCRIPT_PLOTS
 * Expects tuples: (sid, star_log, reads_per_gene?, rsem_genes, rsem_isoforms?)
 * Only 'star_log' and 'rsem_genes' are strictly required here; others may be empty.
 * Pass your R script via --rscript; outputs land in ${params.outdir}/plots.
 */

process rscript_plots {
  tag   "${sid}"
  label 'rplots'

  publishDir "${params.outdir}/plots", mode: 'copy', overwrite: true

  conda """
    channels:
      - conda-forge
      - bioconda
    dependencies:
      - r-base=4.3
      - r-argparse
      - r-optparse
      - r-tidyverse
      - r-data.table
      - r-readr
      - r-ggplot2
      - r-cowplot
      - r-reshape2
  """

  cpus  2
  memory '4 GB'
  time   '2h'

  input:
  tuple val(sid), path(star_log), path(reads_per_gene), path(rsem_genes), path(rsem_isoforms)

  output:
  path "${sid}_*.pdf",  optional: true, emit: pdf
  path "${sid}_*.png",  optional: true, emit: png
  path "${sid}_report.html", optional: true, emit: html

  when:
  params.rscript

  script:
  def rscript = file(params.rscript)
  if (!rscript.exists()) { exit 1, "R script not found: ${params.rscript}" }

  def rpg_arg  = reads_per_gene ? "--reads_per_gene ${reads_per_gene}" : ""
  def riso_arg = rsem_isoforms  ? "--rsem_isoforms ${rsem_isoforms}"   : ""

  """
  Rscript ${rscript} \
    --sid ${sid} \
    --star_log ${star_log} \
    --rsem_genes ${rsem_genes} \
    ${rpg_arg} \
    ${riso_arg} \
    --outdir .
  """
}

