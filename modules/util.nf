process samtools_index{
tag { "samtools index ${bam}" }
publishDir "${params.outdir ?: 'results'}/star/${sid}", 
           mode: 'copy', 
           pattern: "${sid}.*.bai"


input:
tuple val(sid), path(bam)

output:
tuple val(sid), path("${bam}.bai"),	emit: bai

conda (params.use_conda ? 'envs/rna.yml' : null)
container params.container


script:
"""
samtools index -@ ${task.cpus} ${bam}
"""
}
