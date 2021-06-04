FROM continuumio/miniconda3:latest AS condabuild
LABEL authors="Matt Bull" \
      description="Docker image containing all requirements for the ARTIC project's ncov2019 pipeline"

COPY environments/extras.yml /extras.yml
COPY environments/nanopore/environment.yml /environment.yml
RUN /opt/conda/bin/conda install mamba -c conda-forge && \
/opt/conda/bin/mamba env create -f /environment.yml

FROM debian:buster-slim
RUN apt-get update && \
apt-get install -y git procps && \
apt-get clean -y
COPY --from=condabuild /opt/conda/envs/artic /opt/conda/envs/artic
ENV PATH=/opt/conda/envs/artic/bin:$PATH
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENTRYPOINT ["/opt/conda/envs/artic/bin/artic"]
