#!/bin/bash

#install submodules
echo "INSTALLING SUBMODULES..."

#upgrade pip
pip install --upgrade pip

#setup other dependencies
pip install taigapy
cd mnnpy; pip install .; cd ..

#run QC
echo "RUNNING CELLIGNER..."
python "$@"