FROM continuumio/miniconda3:latest
LABEL authors="Matt Bull" \
      description="Docker image containing all requirements for the ARTIC project's ncov2019 pipeline"

COPY extras.yml /extras.yml
RUN apt-get update && apt-get install -y curl g++ git make procps && apt-get clean -y
RUN curl -o /environment.yml -fsSL 'https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/environment.yml'
RUN /opt/conda/bin/conda env create -f /environment.yml
RUN /opt/conda/bin/conda env update -f /extras.yml -n artic-ncov2019
ENV PATH /opt/conda/envs/artic-ncov2019/bin:$PATH
# RUN conda install nextflow
RUN git clone --recurse-submodules https://github.com/artic-network/artic-ncov2019 scheme
ENTRYPOINT ["/opt/conda/envs/artic-ncov2019/bin/artic"]
