FROM nfcore/base:2.1

LABEL \
    author="Abu Saadat" \
    description="Docker image containing phaser tools for ASENext pipeline" \
    maintainer="Abu Saadat"

# Install dependencies
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
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Clone phaser repository
RUN git clone https://github.com/secastel/phaser.git /opt/phaser

# Install phaser requirements
RUN pip3 install --no-cache-dir \
    numpy \
    scipy \
    pysam \
    intervaltree

# Add phaser to PATH
ENV PATH="/opt/phaser/phaser:/opt/phaser/phaser_gene_ae:${PATH}"

# Test installation
RUN phaser.py --help && phaser_gene_ae.py --help

WORKDIR /data/
