Bootstrap: docker
From: continuumio/miniconda3:latest
Stage: condabuild

%files
environments/illumina/environment.yml /environment.yml
environments/extras.yml /extras.yml

%labels
authors="Matt Bull" 
description="Docker image containing all requirements for the ARTIC project's ncov2019 pipeline"

%post
/opt/conda/bin/conda install mamba -c conda-forge && \
/opt/conda/bin/mamba env create -f /environment.yml #&& \
/opt/conda/bin/mamba env update -f /extras.yml -n artic-ncov2019-illumina


Bootstrap: docker
From: debian:buster-slim
Stage: final

%environment
export PATH=/opt/conda/envs/artic-ncov2019-illumina/bin:$PATH
export LC_ALL=C.UTF-8
export LANG=C.UTF-8


%files from condabuild
  /opt/conda/envs/artic-ncov2019-illumina /opt/conda/envs/artic-ncov2019-illumina

%post
apt-get update 
apt-get install -y g++ git make procps rsync locales 
apt-get clean -y
