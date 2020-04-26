#!/bin/bash
set -eo pipefail

echo Install Singularity dependencies.. >> artifacts/test_artifact.log
sudo apt-get update
DEBIAN_FRONTEND=noninteractive sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup

# check that go is already installed in the github runner machine
echo $(which go)
go version

# install Singularity
export VERSION=3.5.3
echo Install Singularity version $VERSION .. >> artifacts/test_artifact.log
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
tar -xzf singularity-${VERSION}.tar.gz
cd singularity
./mconfig
make -C builddir
sudo make -C builddir install
cd ..

echo $(which singularity)
singularity --version
