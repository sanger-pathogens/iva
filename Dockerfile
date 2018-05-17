# This container will install IVA from master
#
FROM debian:testing

# Install required dependancies
RUN apt-get -qq update && apt-get install -y openjdk-8-jdk python3-pip git wget unzip zlib1g-dev libncurses5-dev
ENV JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64"
RUN git clone https://github.com/sanger-pathogens/iva.git
RUN cd iva && ./install_dependencies.sh
ENV PATH /iva/build/kmc-2.3.0:/iva/build/samtools-1.3:/iva/build/smalt-0.7.6-bin:/iva/build/samtools-1.3:/iva/build/MUMmer3.23:/iva/build/SPAdes-3.7.1-Linux/bin:$PATH
RUN export PATH

# Install optional dependencies
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && unzip Trimmomatic-0.38.zip

RUN wget http://ccb.jhu.edu/software/kraken/dl/kraken-1.0.tgz && tar -xvzf kraken-1.0.tgz
RUN cd kraken-1.0 && ./install_kraken.sh ../kraken_install
ENV PATH /kraken_install/:$PATH
RUN export PATH

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.5.0/ncbi-blast-2.5.0+-x64-linux.tar.gz && tar -xzvf ncbi-blast-2.5.0+-x64-linux.tar.gz
ENV PATH /ncbi-blast-2.5.0+:$PATH
RUN export PATH

# Install IVA
RUN cd iva && python3 setup.py install