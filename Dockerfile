# This container will install IVA from master
#
FROM debian:testing

# Install the dependancies
RUN apt-get update -qq && apt-get install -y python3-pip git wget unzip zlib1g-dev libncurses5-dev

RUN git clone https://github.com/sanger-pathogens/iva.git
RUN cd iva && ./install_dependencies.sh
ENV PATH /iva/build/kmc-2.3.0:/iva/build/samtools-1.3:/iva/build/smalt-0.7.6-bin:/iva/build/samtools-1.3:/iva/build/MUMmer3.23:/iva/build/SPAdes-3.7.1-Linux/bin:$PATH
RUN export PATH
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && unzip Trimmomatic-0.38.zip
RUN cd iva && python3 setup.py install