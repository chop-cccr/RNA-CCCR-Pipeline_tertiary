process MCP_COUNTER {
    tag "${matrix.baseName}"

    publishDir "${params.outdir}/immune_deconv", mode: 'copy'

    conda "${projectDir}/envs/immunedeconv.yml"

    input:
    path matrix

    output:
    path "mcp_counter.tsv", emit: mcp_tsv
    path "mcp_counter_heatmap.pdf", emit: mcp_heatmap

    script:
    """
    Rscript ${projectDir}/scripts/run_mcp_counter.R \
      --input ${matrix} \
      --output mcp_counter.tsv \
      --plot mcp_counter_heatmap.pdf
    """
}
