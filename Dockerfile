#!/bin/bash

# Dockerfile to create celligner image
#
# Run build_docker.sh

FROM python:3.9

# install conseq to get the conseq-helper command which conseq needs to run commands inside
# containers and transfer data in and out of the container.
ADD https://github.com/broadinstitute/conseq/archive/refs/tags/v1.29.0.tar.gz /install/conseq.tar.gz
RUN cd /install && \
    tar xzf /install/conseq.tar.gz && \
    python -m venv /opt/conseq && \
    /opt/conseq/bin/python -m pip install poetry && \
    bash -c 'source /opt/conseq/bin/activate  && cd /install/conseq-1.29.0 && poetry install'

#add R and CMAKE
RUN apt-get update && apt-get install -y r-base cmake
# install
RUN R -e 'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos="http://cran.us.r-project.org")};BiocManager::install("limma");' 

#install requirements
COPY requirements.txt .
RUN pip install --upgrade pip &&\
	pip install -r requirements.txt

RUN python -m pip install git+https://github.com/DeKegel/mnnpy.git

RUN mkdir -p "/var/env/script"

COPY run_celligner.py /var/env/script/
COPY celligner /var/env/script/celligner/


WORKDIR "/var/env/script"
