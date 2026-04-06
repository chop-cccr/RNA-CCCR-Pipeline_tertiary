nextflow.enable.dsl=2

process FASTQC {
    tag "${sid}"
    label 'fastqc'

     publishDir "${params.outdir}/fastqc/${sid}", mode: 'copy', overwrite: true


    input:
    tuple val(sid), path(r1), path(r2), val(rg)
    //tuple val(sample), path(fq)

    output:
    path "${sid}_fastqc.html", emit: html
    path "${sid}_fastqc.zip",  emit: zip

    script:
    def readIn = "${r1}"


    //conda     'bioconda::fastqc=0.12.1'
    //container 'biocontainers/fastqc:v0.12.1_cv8'

    script:
    """
    fastqc  --outdir . ${readIn}

    # FastQC names outputs based on the FASTQ basename; standardize to sample-based
    html_file=\$(ls *_fastqc.html | head -n 1)
    zip_file=\$(ls *_fastqc.zip  | head -n 1)

    mv "\$html_file" "${sid}_fastqc.html"
    mv "\$zip_file"  "${sid}_fastqc.zip"
    """
}

