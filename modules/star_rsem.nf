// STAR emits both genomic (QC) and transcriptome BAMs

process STAR_ALIGN {
  tag { "STAR ${sid}" }
  publishDir "${params.outdir ?: 'results'}/star/${sid}",
             mode: 'copy',
             pattern: "*.{bam,Log.final.out,SJ.out.tab,ReadsPerGene.out.tab}"

  input:
  tuple val(sid), path(r1), path(r2), val(rg)
  val star_index

  output:
  tuple val(sid), path("${sid}.Aligned.sortedByCoord.out.bam"),   emit: bam
  tuple val(sid), path("${sid}.Aligned.toTranscriptome.out.bam"), emit: tx_bam
  path "${sid}.Log.final.out",                                    emit: logs
  path "${sid}.SJ.out.tab",                                       emit: splice
  path "${sid}.ReadsPerGene.out.tab", 				  emit: genecounts

  cpus   (params.star_cpus ?: 8)
  memory (params.star_mem  ?: '32 GB')

  script:
  def readFilesCmd = r1.name.endsWith('.gz') ? '--readFilesCommand zcat' : ''
  def readIn       = r2 ? "${r1} ${r2}" : "${r1}"

  """
  set -euo pipefail
  STAR \\
    --runThreadN ${task.cpus} \\
    --genomeDir ${star_index} \\
    --readFilesIn ${readIn} \\
    ${readFilesCmd} \\
    --outFileNamePrefix ${sid}. \\
    --outSAMtype BAM SortedByCoordinate \\
    --quantMode TranscriptomeSAM GeneCounts
  """
}




// RSEM from STAR’s transcriptome BAM


process RSEM_CALCULATE {
  tag { "RSEM ${sid}" }
  publishDir "${params.outdir ?: 'results'}/rsem/${sid}",
             mode: 'copy',
             pattern: "${sid}.*.results"

  input:
  tuple val(sid), path(tx_bam)     // ★ transcriptome BAM from STAR
  val rsem_ref                     // prefix, e.g. /refs/rsem/mm10/mm10

  output:
  tuple val(sid), path("${sid}.genes.results"), path("${sid}.isoforms.results")

  cpus   (params.rsem_cpus ?: 8)
  memory (params.rsem_mem  ?: '32 GB')
  // conda "rsem=1.3.3 samtools=1.19"
  // or container 'quay.io/biocontainers/rsem:1.3.3--pl5321h7ff8a90_7'

  /*
   * Compute flags in Groovy (NOT bash) to avoid 'bad substitution'
   */
  script:
  // Strandedness mapping
  def stranded = (params.strandedness ?: 'unstranded').toString()
  def STRAND_FLAG = stranded == 'forward' ? '--forward-prob 1.0'
                   : stranded == 'reverse' ? '--forward-prob 0.0'
                   : ''

  """
  set -euo pipefail

  # Auto-detect paired-end from BAM flags (bit 0x1)
  if [ \$(samtools view -c -f 1 ${tx_bam}) -gt 0 ]; then
    PE_FLAG="--paired-end"
  else
    PE_FLAG=""
  fi

  rsem-calculate-expression \\
    --bam \\
    --no-bam-output \\
    --paired-end \\
    --num-threads ${task.cpus} \\
    ${tx_bam} \\
    ${rsem_ref} \\
    ${sid}
  """
}

