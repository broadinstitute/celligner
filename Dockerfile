#!/bin/bash

# Dockerfile to create celligner image
#
# Run build_docker.sh

FROM python:3.8

#add R and CMAKE
RUN apt-get update && apt-get install -y r-base cmake
# install
RUN R -e 'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos="http://cran.us.r-project.org")};BiocManager::install("limma");' 

#install requirements
COPY requirements.txt .
RUN pip install --upgrade pip &&\
	pip install -r requirements.txt