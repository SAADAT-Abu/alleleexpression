FROM nfcore/base:2.1

LABEL \
    author="Abu Saadat" \
    description="Docker image containing phaser tools for ASENext pipeline" \
    maintainer="Abu Saadat"

# Install system dependencies including all bioinformatics tools phaser might need
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        python3-dev \
        python3-pip \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        samtools \
        bcftools \
        tabix \
        bedtools \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install phaser requirements
RUN pip3 install --no-cache-dir \
    numpy \
    scipy \
    pysam \
    intervaltree \
    pandas

# Clone phaser repository
RUN git clone https://github.com/secastel/phaser.git /opt/phaser

# Create wrapper scripts in /usr/local/bin
RUN echo '#!/bin/bash\npython3 /opt/phaser/phaser/phaser.py "$@"' > /usr/local/bin/phaser.py && \
    echo '#!/bin/bash\npython3 /opt/phaser/phaser_gene_ae/phaser_gene_ae.py "$@"' > /usr/local/bin/phaser_gene_ae.py && \
    chmod +x /usr/local/bin/phaser.py && \
    chmod +x /usr/local/bin/phaser_gene_ae.py

# Test installation
RUN phaser.py --help && phaser_gene_ae.py --help

WORKDIR /data/
