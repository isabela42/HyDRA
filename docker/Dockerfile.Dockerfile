FROM ubuntu:22.04

LABEL maintainer="mb.isabela42@gmail.com"

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    tar \
    curl \
    gcc \
    python3 \
    python3-pip \
    python2 \
    perl \
    default-jre \
    git \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    libssl-dev \
    pkg-config \
    libfreetype6-dev \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy test dataset into container
COPY . .

# Install tools that require downloading with wget 
RUN bash tools/wget_requirements.sh
ENV PATH="/tools/bbmap:$PATH"
ENV PATH="/opt/conda/bin:$PATH"
RUN echo '#!/bin/bash\njava -jar /usr/local/bin/trimmomatic.jar "$@"' > /usr/local/bin/trimmomatic && \
    chmod +x /usr/local/bin/trimmomatic

# Install a Bioconda packages
RUN conda install -y mamba -n base -c conda-forge
RUN bash tools/conda_requirements.sh
RUN /opt/conda/bin/conda run --prefix /opt/conda-envs/ezlncpred pip3 install ezlncpred
ENV PATH="/opt/conda-envs/busco/bin:\
/opt/conda-envs/cutadapt/bin:\
/opt/conda-envs/deeptools/bin:\
/opt/conda-envs/fmlrc2/bin:\
/opt/conda-envs/multiqc/bin:\
/opt/conda-envs/nanopack/bin:\
/opt/conda-envs/picard/bin:\
/opt/conda-envs/rseqc/bin:\
/opt/conda-envs/transrate/bin:\
/opt/conda-envs/joint/bin:\
/opt/conda-envs/bedtools/bin:\
/opt/conda-envs/bedops/bin:\
/opt/conda-envs/ucsc_tools/bin:\
/opt/conda-envs/fastqpair/bin:\
/opt/conda-envs/blast/bin:\
/opt/conda-envs/cdhit/bin:\
/opt/conda-envs/gmap/bin:\
/opt/conda-envs/pblat/bin:\
/opt/conda-envs/trinity/bin:\
/opt/conda-envs/spades/bin:\
/opt/conda-envs/feelnc/bin:\
/opt/conda-envs/ezlncpred/bin:\
$PATH"
ENV AUGUSTUS_CONFIG_PATH=/opt/conda/pkgs/augustus-*/config/

# Install tools that require git clone 
RUN bash tools/git_requirements.sh

# Set entrypoint for container (default)
ENTRYPOINT ["/bin/bash", "-l", "-c"]