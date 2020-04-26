FROM continuumio/miniconda3:latest
LABEL authors="Matt Bull" \
      description="Docker image containing all requirements for an Illumina ncov2019 pipeline"

COPY environments/extras.txt /extras.txt
COPY environments/illumina/environment.yml /environment.yml
RUN apt-get update && apt-get install -y curl g++ git make procps && apt-get clean -y
RUN /opt/conda/bin/conda env create -f /environment.yml 
RUN /opt/conda/envs/artic-ncov2019-illumina/bin/pip install -r /extras.txt && /opt/conda/bin/conda clean -a
ENV PATH=/opt/conda/envs/artic-ncov2019-illumina/bin:$PATH
