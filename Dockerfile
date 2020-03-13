# This container will install IVA from master
#
FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /tmp

ENV INSTALL_DIRECTORY=/opt


# Install required dependancies
RUN apt-get update -y -qq \
    && apt-get install -y -qq \
      openjdk-8-jdk \
      python3-pip \
      git \
      wget \
      unzip \
      zlib1g-dev \
      libncurses5-dev \
      libbz2-dev \
      liblzma-dev
  
ENV JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64"

ARG KMC_VERSION=3.0.0
ARG MUMMER_VERSION=3.23
ARG SAMTOOLS_VERSION=1.3
ARG SMALT_VERSION=0.7.6
ARG TRIMMOMATIC_VERSION=0.38
ARG KRAKEN_VERSION=1.0
ARG BLAST_VERSION=2.5.0

# kmc
RUN mkdir /opt/kmc-${KMC_VERSION} \
    && cd /opt/kmc-${KMC_VERSION} \
    && wget --progress=dot:giga https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC3.linux.tar.gz \
    && tar xf KMC3.linux.tar.gz \
    && rm KMC3.linux.tar.gz \
    && chmod -R 755 /opt/kmc-${KMC_VERSION}
ENV PATH=/opt/kmc-${KMC_VERSION}:$PATH

# MUMmer
RUN cd /opt \
    && wget --progress=dot:giga "http://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz" \
    && tar xf MUMmer${MUMMER_VERSION}.tar.gz \
    && rm MUMmer${MUMMER_VERSION}.tar.gz \
    && cd MUMmer${MUMMER_VERSION} \
    && make
ENV PATH=/opt/MUMmer${MUMMER_VERSION}:$PATH

# samtools
RUN cd /opt \
    && wget --progress=dot:giga "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    && tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && make
ENV PATH=/opt/samtools-${SAMTOOLS_VERSION}:$PATH
 
# smalt
RUN cd /opt \
    && wget --progress=dot:giga http://downloads.sourceforge.net/project/smalt/smalt-${SMALT_VERSION}-bin.tar.gz \
    && tar xf smalt-${SMALT_VERSION}-bin.tar.gz \
    && rm smalt-${SMALT_VERSION}-bin.tar.gz \
    && cd smalt-${SMALT_VERSION}-bin \
    && ln -fs smalt_x86_64 smalt
ENV PATH=/opt/smalt-${SMALT_VERSION}-bin:$PATH
    
# Trimmomatic
RUN cd / \
    && wget --progress=dot:giga http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
    && rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip

# Kraken
RUN mkdir /tmp/KRAKEN \
    && cd /tmp/KRAKEN \
    && wget --progress=dot:giga http://ccb.jhu.edu/software/kraken/dl/kraken-${KRAKEN_VERSION}.tgz \
    && tar -xvzf kraken-${KRAKEN_VERSION}.tgz \
    && cd kraken-${KRAKEN_VERSION} \
    && ./install_kraken.sh /opt/kraken-${KRAKEN_VERSION} \
    && rm -rf /tmp/KRAKEN
ENV PATH=/opt/kraken-${KRAKEN_VERSION}/:$PATH

# ncbi blast
RUN cd /opt \
    && wget --progress=dot:giga ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && tar -xf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
    && rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
ENV PATH=/opt/ncbi-blast-${BLAST_VERSION}+/bin:$PATH

# Install IVA
ARG BUILD_DIR=/tmp/IVA
COPY . $BUILD_DIR
RUN cd ${BUILD_DIR} \
    && pip3 install cython \
    && python3 setup.py test \
    && python3 setup.py install \
    && rm -rf ${BUILD_DIR}
