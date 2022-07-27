FROM ubuntu:xenial

# install needed software for conda install; cleanup apt garbage
RUN apt-get update && apt-get install -y \
 wget \
 git \
 build-essential && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/* 

# get miniconda and the pangolin repo
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
 bash ./Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b  && \
 rm Miniconda3-latest-Linux-x86_64.sh && \
 wget "https://github.com/cov-lineages/pangolin/archive/refs/tags/v3.1.14.tar.gz" && \
 tar -xf v3.1.14.tar.gz && \
 rm v3.1.14.tar.gz && \
 mv -v pangolin-* pangolin

# set the environment
ENV PATH="/miniconda/bin:$PATH" 
ENV LC_ALL C.UTF-8 
ENV LANG C.UTF-8

# modify environment.yml to pin specific versions during install
# create the conda environment using modified environment.yml and set as default
RUN sed -i "s|usher.*|usher=${USHER_VER}|" /pangolin/environment.yml && \
 sed -i "s|pangoLEARN.git|pangoLEARN.git@${PANGOLEARN_VER}|" /pangolin/environment.yml && \
 sed -i "s|scorpio.git|scorpio.git@${SCORPIO_VER}|" /pangolin/environment.yml && \
 sed -i "s|constellations.git|constellations.git@${CONSTELLATIONS_VER}|" /pangolin/environment.yml && \
 sed -i "s|pango-designation.git|pango-designation.git@${PANGO_DESIGNATION_VER}|" /pangolin/environment.yml && \
 conda env create -f /pangolin/environment.yml && \
 echo "source activate pangolin" > /etc/bash.bashrc
ENV PATH /miniconda/envs/pangolin/bin:$PATH

# took this out from below step: python setup.py install (seems like the install instructions changed?)
# replaced with pip install .
RUN cd pangolin && \
 pip install . && \
 conda clean -a -y && \
 mkdir /data && \
 pangolin -v && pangolin -pv && pangolin -dv

WORKDIR /data