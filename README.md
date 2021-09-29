# Celligner

![](docs/typical_celligner.webp)

celligner is a computational project to align multiple cancer datasets across sequencing modalities, tissue conditions (media, perturbations..) and format (CL/tumor/organoids/spheroids)

see our latest paper on aligning CCLE cell lines with TCGA tumors:
[2020 paper](https://www.nature.com/articles/s41467-020-20294-x)


## Install

> TO see the old R package installation instruction, see the `R/`folder.

`pip install celligner`

a dockerized version is available at `jkobject:pycelligner`

## for developers

see `CONTRIBUTING.md`

## run Celligner

see `Celligner_demo.ipynb` for an example of usage.
  
## Multidataset alignment

see `Celligner_demo.ipynb` for an example of usage.

one can use addToFit(), addToPredict() depending on whether they want to align their dataset to another or align another dataset to theirs.

if you have a very small dataset and want to align to CCLE or CGA, use the parameter `doAdd=True` in order to not rerun the entire pipeline and use cached information.


Please use _github issues_ for any problem related to the tool.

Maintainers:

Jérémie Kalfon @jkobject
James McFarland
Javad Noorbak @jnoorbak

Created by: 

Allie Warren @awarren

