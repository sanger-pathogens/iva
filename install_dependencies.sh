#!/bin/bash
set -e
set -x

start_dir=$(pwd)

KMC_VERSION=3.0.0
MUMMER_VERSION=3.23
SAMTOOLS_VERSION=${SAMTOOLS_VERSION:-"1.3"}
BOWTIE2_VERSION=2.2.8

#KMC_DOWNLOAD_URL=http://sun.aei.polsl.pl/REFRESH/kmc/downloads/${KMC_VERSION}/linux/kmc
#KMCDUMP_DOWNLOAD_URL=http://sun.aei.polsl.pl/REFRESH/kmc/downloads/${KMC_VERSION}/linux/kmc_dump
KMC3_DOWNLOAD_URL="https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC3.linux.tar.gz"
MUMMER_DOWNLOAD_URL="http://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz"
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
BOWTIE2_DOWNLOAD_URL=https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip

# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}



# --------------- kmc -----------------
kmc_dir="$build_dir/kmc-${KMC_VERSION}"
rm -fr $kmc_dir
mkdir $kmc_dir
cd $kmc_dir
#download $KMC3_DOWNLOAD_URL "kmc"
#download $KMCDUMP_DOWNLOAD_URL "kmc_dump"
#chmod +x kmc kmc_dump

download "${KMC3_DOWNLOAD_URL}" "KMC3.linux.tar.gz"
tar xzf KMC3.linux.tar.gz
chmod +x kmc
chmod +x kmc_tools
chmod +x kmc_dump


# --------------- mummer ------------------
cd $build_dir
download $MUMMER_DOWNLOAD_URL "MUMmer${MUMMER_VERSION}.tar.gz"
mummer_dir="$build_dir/MUMmer${MUMMER_VERSION}"
tar -zxf MUMmer${MUMMER_VERSION}.tar.gz
cd $mummer_dir
make


# --------------- samtools -----------------
cd $build_dir
download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd $samtools_dir
make


# --------------- bowtie2 -----------------
cd $build_dir
download $BOWTIE2_DOWNLOAD_URL "bowtie2-${BOWTIE2_VERSION}.zip"
unzip -qq bowtie2-${BOWTIE2_VERSION}.zip
bowtie2_dir="$build_dir/bowtie2-${BOWTIE2_VERSION}"
rm bowtie2-${BOWTIE2_VERSION}.zip


cd $start_dir

update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${kmc_dir}
update_path ${mummer_dir}
update_path ${samtools_dir}
update_path ${bowtie2_dir}


