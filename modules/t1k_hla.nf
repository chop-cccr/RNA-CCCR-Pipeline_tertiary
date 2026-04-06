process T1K_TYPING {
  tag "$sample_id"
  cpus 4
  publishDir "${params.outdir ?: 'results'}/T1K_HLA/${sample_id}",
             mode: 'copy', overwrite : true


  conda 'bioconda::t1k=1.0'

  input:
    tuple val(sample_id), path(r1), path(r2)

  output:
    tuple val(sample_id), path("${sample_id}.t1k_genotype.tsv"), emit: genotype
    path("${sample_id}.t1k_out"), type: 'dir', emit: outdir

  script:
  """
  set -euo pipefail

  run-t1k \
    -1 ${r1} -2 ${r2} \
    --preset hla \
    -f ${params.t1k_hla_fasta} \
    -t ${task.cpus} \
    --od ${sample_id}.t1k_out \
    -o ${sample_id}


 cp ${sample_id}.t1k_out/${sample_id}_genotype.tsv ${sample_id}.t1k_genotype.tsv
 """
}

