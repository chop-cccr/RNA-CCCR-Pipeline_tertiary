process PREP_STAR_INDEX {
output:
path 'star_index', emit: star_index


when:
params.build_index


conda (params.use_conda ? 'envs/rna.yml' : null)
container params.container


script:
def fasta = reference.fasta
def gtf = reference.gtf
def sjdb = params.sjdbOverhang ?: 100
def threads = task.cpus
def out = 'star_index'
"""
mkdir -p ${out}
STAR \
--runThreadN ${threads} \
--runMode genomeGenerate \
--genomeDir ${out} \
--genomeFastaFiles ${fasta} \
--sjdbGTFfile ${gtf} \
--sjdbOverhang ${sjdb}
"""
}


process PREP_RSEM_INDEX {
tag "RSEM reference (${reference.species})"
publishDir params.outdir, mode: 'copy', pattern: 'rsem_ref/**'


input:
val reference


output:
path 'rsem_ref', emit: rsem_ref


when:
params.build_index


conda (params.use_conda ? 'envs/rna.yml' : null)
container params.container


script:
def fasta = reference.fasta
def gtf = reference.gtf
def prefix = 'rsem_ref/' + (reference.species == 'human' ? 'GRCh38' : reference.species == 'mouse' ? 'GRCm39' : 'custom')
"""
mkdir -p rsem_ref
rsem-prepare-reference \
--gtf ${gtf} \
--star \
${fasta} \
${prefix}
"""
}
