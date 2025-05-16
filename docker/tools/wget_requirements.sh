#!/bin/bash

cd tools/

# Miniconda
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -u -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    conda clean -afy

# pip2 is deprecated in Ubuntu 22.04
wget -O get-pip.py https://bootstrap.pypa.io/pip/2.7/get-pip.py && \
    python2 get-pip.py && \
    rm get-pip.py

# FastQC v0.12.1
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod +x FastQC/fastqc && \
    ln -s "$(pwd)/FastQC/fastqc" /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip

# BBMap suite v39.01
wget https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz && \
    tar -xzf BBMap_39.01.tar.gz && \
    rm BBMap_39.01.tar.gz

# Chopper-0.5.0
wget https://github.com/wdecoster/chopper/releases/download/v0.5.0/chopper-linux.zip && \
    unzip chopper-linux.zip && \
    chmod +x chopper && \
    mv chopper /usr/local/bin/

# Trimmomatic v0.36
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
    unzip Trimmomatic-0.36.zip && \
    mv Trimmomatic-0.36 /opt/trimmomatic && \
    ln -s /opt/trimmomatic/trimmomatic-0.36.jar /usr/local/bin/trimmomatic.jar && \
    rm Trimmomatic-0.36.zip

