#!/bin/bash

cd tools/

# Seqtk v1.3-r106
git clone https://github.com/lh3/seqtk.git --branch v1.3 && \
    cd seqtk && \
    make && \
    chmod +x seqtk && \
    ln -s "$(pwd)/seqtk" /usr/local/bin/seqtk && \
    cd ../

# Rcorrector v1.0.5 - problematic
git clone https://github.com/mourisl/Rcorrector --branch v1.0.5 && \
    cd Rcorrector && \
    make && \
    cd ../

# FilterUncorrectabledPEfastq.py - HyDRA was developed with Python2 version, this runs on Python3
git clone https://github.com/harvardinformatics/TranscriptomeAssemblyTools.git

# Porechop v0.2.4
git clone https://github.com/rrwick/Porechop.git && \
    cd Porechop && \
    python3 setup.py install && \
    chmod +x porechop && \
    cd ../

# RopeBWT2 (r187)
git clone https://github.com/lh3/ropebwt2.git && \
    cd ropebwt2 && \
    make && \
    ln -s "$(pwd)/ropebwt2" /usr/local/bin/ropebwt2 && \
    cd ../


