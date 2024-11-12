#!/bin/bash

set -x
set -e

start_dir=$(pwd)

MUMMER_VERSION="3.23"

MUMMER_DOWNLOAD_URL="http://sourceforge.net/projects/mummer/files/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz/download"

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

download $MUMMER_DOWNLOAD_URL "mummer-${MUMMER_VERSION}.tgz"

# Build all the things
cd $build_dir

## Mummer
mummer_dir=$(pwd)/MUMmer${MUMMER_VERSION}
if [ ! -d $mummer_dir ]; then
  tar xzf mummer-${MUMMER_VERSION}.tgz
fi
cd $mummer_dir
if [ -e "${mummer_dir}/mummer" ]; then
  echo "Already built Mummer; skipping build"
else
  make check
  make
  cd ./src/tigr 
  make
fi

cd $build_dir

# Setup environment variables
update_path () {
  new_dir=$1
  export PATH=${PATH:-$new_dir}
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${mummer_dir}
update_path ${mummer_dir}/src/tigr

cd $start_dir

set +x
set +e
