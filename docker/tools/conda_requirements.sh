#!/bin/bash

# Channels
CHANNELS="-c bioconda -c conda-forge"

# Quality control
conda create --yes --prefix /opt/conda-envs/busco $CHANNELS busco=2.0
conda create --yes --prefix /opt/conda-envs/cutadapt $CHANNELS cutadapt=4.4
conda create --yes --prefix /opt/conda-envs/deeptools $CHANNELS deeptools=3.5.4
conda create --yes --prefix /opt/conda-envs/fmlrc2 $CHANNELS fmlrc2=0.1.7
conda create --yes --prefix /opt/conda-envs/multiqc $CHANNELS multiqc=1.14
conda create --yes --prefix /opt/conda-envs/nanopack $CHANNELS nanopack
conda create --yes --prefix /opt/conda-envs/picard $CHANNELS picard=2.19
conda create --yes --prefix /opt/conda-envs/rseqc $CHANNELS rseqc=2.6.4
conda create --yes --prefix /opt/conda-envs/transrate $CHANNELS transrate=1.0.3

# BAM & BED tools
conda create --yes --prefix /opt/conda-envs/joint $CHANNELS samtools=1.9 htslib=1.9 bowtie2=2.2.9 minimap2=2.16
conda create --yes --prefix /opt/conda-envs/bedtools $CHANNELS bedtools=2.29.0
conda create --yes --prefix /opt/conda-envs/bedops $CHANNELS bedops=2.4.41
conda create --yes --prefix /opt/conda-envs/ucsc_tools $CHANNELS ucsc_tools
conda create --yes --prefix /opt/conda-envs/fastqpair $CHANNELS fastq-pair=0.3

# Aligners
conda create --yes --prefix /opt/conda-envs/blast $CHANNELS blast=2.4.0
conda create --yes --prefix /opt/conda-envs/cdhit $CHANNELS cd-hit=4.6.8
conda create --yes --prefix /opt/conda-envs/gmap $CHANNELS gmap=2023.07.20
conda create --yes --prefix /opt/conda-envs/pblat $CHANNELS pblat=2.5.1

# Assemblers
conda create --yes --prefix /opt/conda-envs/trinity $CHANNELS trinity=2.8.4
conda create --yes --prefix /opt/conda-envs/spades $CHANNELS spades=3.14.1

# lncRNAs
conda create --yes --prefix /opt/conda-envs/feelnc $CHANNELS feelnc=0.2
conda create --yes --prefix /opt/conda-envs/ezlncpred $CHANNELS python=3.6 samtools bcftools pysam numpy