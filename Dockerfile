# Example pattern – adapt to your existing Dockerfile
FROM mambaorg/micromamba:1.5.8

USER root
RUN apt-get update && apt-get install -y \
      openjdk-17-jre-headless \
      pigz \
      gzip \
      bzip2 \
      curl \
      wget \
      perl \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Create env with required tools
RUN micromamba create -y -n rnaseq-env -c conda-forge -c bioconda \
    fastqc \
    star \
    rsem \
    samtools \
    && micromamba clean -a -y





ENV MAMBA_DEFAULT_ENV=rnaseq-env
ENV PATH=/opt/conda/envs/rnaseq-env/bin:$PATH

WORKDIR /workspace

