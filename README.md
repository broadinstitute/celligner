# Celligner

![](docs/typical_celligner.webp)

celligner is a computational project to align multiple cancer datasets across sequencing modalities, tissue conditions (media, perturbations..) and format (CL/tumor/organoids/spheroids)

see our latest paper on aligning CCLE cell lines with TCGA tumors:
[2020 paper](https://www.nature.com/articles/s41467-020-20294-x)

## Remark

Celligner is initially an R project that you can find in the `R/` folder.

A python version was made that is the exact same version as the R verion. However one should not expect the same plot.

### Why?

The plot some users have been used to is a unique run of UMAP on the celligner realignment data. This is done by fixing the seed of the UMAP algorithm. You can still do that for the python version but it is disabled by default and not recommended. We recommend users to play with the UMAP parameter and make multiple plots. 
This helps to prevent reading too much into UMAP's output. Things that don't stay the same are not necessarily true attributes of the data.

Additionally we also advice users to completement assumptions by applying methods like differential expression analysis across clusters to find any meaningful information.

### Is there any algorithmic difference?

Celligner is composed of 4 key steps:

1. A Louvain clustering: this version is the ScanPy implementation of this method while Celligner is using Seurat's. There might be some slight implementation differences.
2. A limma diff expression analysis to find key variance genes across clusters for each dataset: this version is 100% similar to the R version of celligner.
3. A cPCA to remove tumor impurity signal. This method is exactly the same except that the python version does exact PCA computation while the R version does an approximate version.
4. An MNN allignment: this version is 100% similar to the R version of celligner in its output.

### Is there any other differences?

Overall improvements yes:

1. A “pre-fitted” model is available to download here: gs://celligner/model.pkl (on request for now)
2. Using your own dataset and adding new dataset is super simple now with fit(), transform() syntax
3. You don’t need to rerun the entire model when adding new (adding 600 new samples take only 5mns to run)
4. The model takes much less memory to run and can run on any machine now (you don’t need 64Gb of RAM anymore), and it also takes less than an hour to fully run (on a good machine).
5. There is now an interactive plot using bokeh to better visualise your samples of interest.
6. You can now easily choose parameters and even choose between 2 different versions of MNN.


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

