#!/bin/bash

#install submodules
echo "INSTALLING SUBMODULES..."

#upgrade pip
pip install --upgrade pip


#setup other dependencies
# pip install taigapy
# echo "installing mnnpy"
# cd mnnpy; pip install .; cd ..


#run QC
echo "RUNNING CELLIGNER..."
python "$@" --input inputs.json