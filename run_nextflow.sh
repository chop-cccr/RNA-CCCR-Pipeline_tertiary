#!/bin/sh
#If reference folder was downloaded in the same directory 
#as nextflow. If not, change the --genomeDir and --rsem_ref
    nextflow run main_final.nf \
    -c  ./nextflow.ALIGN.config \
    --samplesheet ./assets/samplesheet.example.csv \
    --read_group TEST  \
    --outdir ./TEST_HPC_runs \
    --species Human \
    --t1k_hla_fasta /mnt/isilon/dbhi_bfx/mishrap/CCCR/RNA-FEB/t1k_ref/hlaidx/hlaidx_rna_seq.fa \
    --genomeDir /mnt/isilon/dbhi_bfx/mishrap/CCCR/RNA-FEB/RNA-CCCR-Pipeline_WITH_T1K/reference/Human/hg38_gencode_v39_may2024/  \
    --rsem_ref /mnt/isilon/dbhi_bfx/mishrap/CCCR/RNA-FEB/RNA-CCCR-Pipeline_WITH_T1K/reference/Human/rsem_ensemble_index/rsem_ensemble  \
    --threads 10 -profile slurm,singularity --resume 

