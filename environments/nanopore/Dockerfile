FROM continuumio/miniconda3:latest
LABEL authors="Matt Bull" \
      description="Docker image containing all requirements for the ARTIC project's ncov2019 pipeline"

COPY environments/extras.txt /extras.txt
RUN apt-get update && apt-get install -y curl g++ git make procps && apt-get clean -y
RUN curl -o /environment.yml -fsSL 'https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/environment.yml'
RUN /opt/conda/bin/conda env create -f /environment.yml 
RUN /opt/conda/envs/artic-ncov2019/bin/pip install -r /extras.txt && /opt/conda/bin/conda clean -a
ENV PATH /opt/conda/envs/artic-ncov2019/bin:$PATH
ENTRYPOINT ["/opt/conda/envs/artic-ncov2019/bin/artic"]
