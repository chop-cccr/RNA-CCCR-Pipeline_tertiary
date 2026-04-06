nextflow.enable.dsl=2

/*
 * Module: RMATS_SINGLE_BAM
 * Purpose: Run rMATS (turbo) on ONE coordinate-sorted BAM + GTF
 * Mode: single-sample, no stats (uses --statoff)
 *
 * Input channel item:
 *   tuple val(sid), path(bam), path(gtf)
 *
 * Outputs:
 *   - Directory with rMATS results (SE, A5SS, A3SS, MXE, RI, etc.)
 *   - Emitted as: outdir (path)
 *
 * Important:
 *   - BAM should be coordinate-sorted; .bai alongside is recommended
 *   - Set strandness via params.rmats_libType:
 *       fr-unstranded | fr-firststrand | fr-secondstrand
 */

params.rmats_outdir              = params.rmats_outdir              ?: "rmats_out"
params.rmats_threads             = params.rmats_threads             ?: 8
params.rmats_readType            = params.rmats_readType            ?: "paired"     // "paired" or "single"
params.rmats_readLen             = params.rmats_readLen             ?: 100
params.rmats_variableReadLength  = params.rmats_variableReadLength  ?: false        // if true, ignore readLen
params.rmats_libType             = params.rmats_libType             ?: "fr-unstranded"
params.rmats_additionalArgs      = params.rmats_additionalArgs      ?: ""           // e.g. "--novelSS --mil 1 --mel 1"
params.rmats_container           = params.rmats_container           ?: null         // e.g. "quay.io/biocontainers/rmats:4.3.0--py39..."
params.rmats_conda               = params.rmats_conda               ?: "bioconda::rmats=4.3.0"

process RMATS_single_sample {

  tag { sid }
  cpus params.rmats_threads
  memory { 8.GB * Math.max(1, (task.cpus.intdiv(4))) }
  time '24h'

  // dynamic publishDir per sample
  publishDir { "${params.rmats_outdir}/${sid}" }, mode: 'copy'

  // Prefer container if provided, else conda
  if (params.rmats_container) {
    container params.rmats_container
  } else {
    conda params.rmats_conda
  }

  input:
  tuple val(sid), path(bam), path(gtf)

  output:
  // emit the per-sample result directory
  path "rmats_${sid}", emit: outdir

  script:
  def varRL = params.rmats_variableReadLength ? "--variable-read-length" : "--readLength ${params.rmats_readLen}"

  """
  set -euo pipefail

  OUTDIR="rmats_${sid}"
  TMPDIR="\${OUTDIR}/tmp"
  mkdir -p "\$OUTDIR" "\$TMPDIR"

  # rMATS expects a list file; one BAM -> one line
  echo "${bam}" > b1.txt

  rmats.py \
    --b1 b1.txt \
    --gtf "${gtf}" \
    --od  "\$OUTDIR" \
    --tmp "\$TMPDIR" \
    -t "${params.rmats_readType}" \
    --nthread ${task.cpus} \
    --libType "${params.rmats_libType}" \
    ${varRL} \
    --statoff \
    --task both \
    ${params.rmats_additionalArgs}
  """
}
