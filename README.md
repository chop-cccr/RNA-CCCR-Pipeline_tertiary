

# STAR + RSEM (human/mouse) – Nextflow DSL2


A simple, reliable pipeline for aligning RNA‑seq reads with [STAR] and quantifying with [RSEM]. Works on  HPC with Docker/Singularity/Conda. Also a built in docker image can be used for deployment.

# Introduction:
STAR (Spliced Transcripts Alignment to a Reference) [1] is an ultrafast splice-aware aligner for RNA-seq that maps reads to genomes, detects canonical/non-canonical splice junctions, and supports two-pass mapping to improve novel junction discovery; it’s commonly used upstream of quantifiers like RSEM. 

RSEM (RNA-Seq by Expectation-Maximization) [2] estimates gene- and isoform-level expression from RNA-seq data, modeling fragment/length effects and outputting TPM/FPKM and expected counts. Works with reference transcriptomes (or de novo assemblies) and supports single/paired-end and stranded protocols

Using this pipeline, user can either create new indexes of genomes or use pre-existing ones. 


## Features
- ✅ Human or mouse presets; or fully custom references
- ✅ Optional **on‑the‑fly index** building from FASTA+GTF
- ✅ Slurm‑ready config; resuming supported (`-resume`)
- ✅ Paired‑end or single‑end reads
- ✅ Containers (Docker/Singularity) or Conda env

## Steps to make sure it runs smoothly
1. Make sure you have loaded singularty.
2. Update conda to activate FASTQC environment.
3. Make sure you have also loaded Java version higher or similar to mentioned in the nextflow instructions manual
4. Make sure you have installed nextflow. You can follow instructions on the nextflow webpage:https://www.nextflow.io/docs/latest/install.html
5. Add the path where nextflow is to your bash (Follow instructions in Step 2).
    For bash :
    export PATH= "NEXTFLOW_INSTALLATION_PATH":$PATH
6. User can input unzipped FASTQ file or gzipped files. The pipeline will detect it and make commands accordingly.
7. FASTQC can be installed using conda. This is the easiest way to make sure FASTQC works.
8. **Reference download instructions:**
	- Run the `fetch_reference.sh` script to copy reference files.
	- You can specify a destination path, or it will copy to the current directory.
	- **Command:** `./fetch_reference.sh`


## SETUP before running nextflow pipeline
Load the following modules. These are necessary modules and you may run into errors if any of them are not loaded properly.
```bash

# To install and run FASTQC
*Install FastQC using envs/fastqc.yml using these commands*

conda  env create -f fastqc.yml
conda activate FastQC_nf

#FASTQC can be installed via conda:
conda install bioconda::fastqc

#Load Java for singularity
module load Java-17.0.6

#Export nextflow PATH
export PATH= "NEXTFLOW_INSTALLATION_PATH":$PATH

#To run STAR and RSEM
module load RSEM
module load STAR
module load samtools
module load singularity

** Also  make sure you load Java after activating FastQC environment

#Load all modules using script on HPC
source load_stuff.sh

will load all above. Please run this after you load FastQC

To run DOCKER: 
nextflow run main_final.nf -profile slurm,singularity \
  --container "docker://ghcr.io/chop-cccr/RNA-CCCR-Pipeline:v0.1.0"

```

Download reference files

```bash
To download reference files for mouse and humans, use script fetch_reference.sh
./fetch_reference.sh 

This will download the reference folder to the currect directory

./fetch_reference <DESTINATION_PATH>

This will download the reference folder to the user defined path

```

## Input files
In the assets folder, you will find an example samplesheet that the pipeline can accept.
``` bash

sample,fastq1,fastq2,read_group
USER_SAMPLE_ID,/mnt/isilon/cccr_bfx/Pipelines/Independent/data/SRR2557083_1.fastq.gz,/mnt/isilon/cccr_bfx/Pipelines/Independent/data/SRR2557083_2.fastq.gz,ID:ANY_ID
```
USER_SAMPLE_ID and ANY_ID is user input and cannot be skipped. They can be same or different.

## Quick start
```bash
## For HPC (Load singularity)

    nextflow run main_final.nf 
    -c  nextflow.ALIGN.config 
    --samplesheet ./assets/samplesheet.example.csv 
    --read_group TEST  
    --outdir ./HELLO_mouse
    --genomeDir ./reference/Mouse/Reference/mm10_gencode/mm10STAR/  
    --rsem_ref ./reference/Mouse/rsem_index/rsem_mm10  
    --threads 10 -profile singularity -resume


** In case of --genomeDir :  add the path to the correct species STAR reference 
              --rsem_ref  :  add the path to the correct species RSEM reference 

** Also make sure your config file paths in nextflow.config points to the directory conf (in the repository)

(Under construction)
With Docker
nextflow run main.nf -profile docker 
  --samplesheet my_samples.csv 
  --species mouse


Inputs
--samplesheet CSV with columns:

sample_id (unique)
fastq_1 (path to R1)
fastq_2 (path to R2)
read_group (optional; default added if missing)


At this point this pipeline only runs PE (Paired end)
--species: human or mouse. For other species, use custom references.

Reference options

Presets (edit in main.nf): hard‑coded cluster paths for GRCh38/GRCm39.

Custom prebuilt indices:

--star_index /path/to/STAR/index \
--rsem_ref /path/to/rsem/prefix   # NOTE: RSEM prefix string, not a folder

Build on the fly from FASTA+GTF:

--build_index --fasta /ref/genome.fa --gtf /ref/genes.gtf [--sjdbOverhang 100]
```
``` bash
Outputs
RESULT_directory/
├── star/
│   ├── <sample>.Aligned.sortedByCoord.out.bam
|   ├── <sample>.Aligned.toTranscriptome.out.bam
│   ├── <sample>.Aligned.sortedByCoord.out.bam.bai
│   └── <sample>.Log.final.out
│   └── <sample>.SJ.out.tab
│   └── <sample>.ReadsPerGene.out.tab
├── rsem/
│   ├── <sample>.genes.results
│   └── <sample>.isoforms.results
├── report/
    ├── summaries/
│   ├── <sample>.rsem.summary.tsv 
├── star_index/        # if build_index=true
└── rsem_ref/          # if build_index=true

```

Notes & tips

Detects .gz reads and uses zcat automatically.



---


## How to adapt to your GitHub repo
1. Drop these files into your repo (or replace your existing `main.nf` & configs).
2. Edit the preset reference paths in `main.nf` under `reference_resolver()` to match your HPC.
3. Build the Docker image (optional) and push to GHCR:
   ```bash
   docker build -t ghcr.io/<youruser>/star-rsem:latest -f containers/Dockerfile .
   docker push ghcr.io/<youruser>/star-rsem:latest

Run with -profile slurm,singularity on HPC or -profile docker locally.






## Reference:
1. Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29(1):15–21 (2013). doi:10.1093/bioinformatics/bts635.
2. Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 12, 323 (2011). doi:10.1186/1471-2105-12-323
