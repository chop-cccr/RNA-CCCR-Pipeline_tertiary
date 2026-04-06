nextflow run main_final.nf -profile slurm \
  -c ./nextflow.ALIGN.config \
  --samplesheet ./assets/samplesheet.example.csv \
  --read_group TEST \
  --outdir ../TEST_HPC_runs \
  --genomeDir /mnt/isilon/cccr_bfx/CCCR_Pipelines/RNA-Seq_Pipeline_HPC/reference/Human/hg38_gencode_v39_may2024/ \
  --rsem_ref /mnt/isilon/cccr_bfx/CCCR_Pipelines/RNA-Seq_Pipeline_HPC/reference/Mouse/rsem_index/rsem_mm10

