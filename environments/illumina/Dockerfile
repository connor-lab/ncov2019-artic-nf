#FROM continuumio/miniconda3:latest
LABEL authors="Matt Bull" \
      description="Docker image containing all requirements for an Illumina ncov2019 pipeline"

COPY environments/extras.yml /extras.yml
COPY environments/illumina/environment.yml /environment.yml
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda
RUN apt-get update && apt-get install -y curl g++ git make procps && apt-get clean -y
RUN /opt/conda/bin/conda env create -f /environment.yml 
RUN /opt/conda/bin/conda env update -f /extras.yml -n artic-ncov2019-illumina && /opt/conda/bin/conda clean -a
ENV PATH=/opt/conda/envs/artic-ncov2019-illumina/bin:$PATH
