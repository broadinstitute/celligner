# Celligner

![](docs/typical_celligner.webp)

__Celligner__ is a computational project to align multiple cancer datasets across sequencing modalities, tissue conditions (media, perturbations..) and format (CL/tumor/organoids/spheroids)

See our latest paper on aligning CCLE cell lines with TCGA tumors:
[2020 paper](https://www.nature.com/articles/s41467-020-20294-x)

## Remark

__Celligner__ is initially an R project that you can find in the `R/` folder.

A Python version was made that performs the same computations as the R version. However one should not expect the exact same plot for a couple reasons:

#### UMAP

The plot shown in the paper / the DepMap portal is the result of a unique run of UMAP on the __Celligner__ aligned data.. This is done by fixing the seed of the UMAP algorithm. You can still do that for the python version but it is disabled by default and not recommended. We recommend users to play with the UMAP parameter and make multiple plots. This helps to prevent reading too much into UMAP's output. Things that don't stay the same are not necessarily true attributes of the data.

Learn more here: [distill](https://distill.pub/2016/misread-tsne/), [Lior's twittorial](https://twitter.com/lpachter/status/1431325969411821572).

#### Algorithmic differences

There are likely some slight implementation differences in the Louvain clustering and contrastive PCA steps.

## Overview

 A reference expression dataset (CCLE cell lines in the original implementation) should be fit using the `fit()` function, and a target expression dataset (TCGA+ tumor samples) can then be aligned to this reference using the `transform()` function. 

Go here for the production version: [https://depmap.org/portal/celligner/](https://depmap.org/portal/celligner/)

See the `run_celligner.py` script for example usage.

## Install

> To see the old R package installation instruction, see the `R/`folder.

Before running pip, make sure that you have R installed.

`pip install celligner`

Even with R, some platform might not have all the required packages already installed (thanks R for being so easy to work with!)

In that case, please refer to our docker image:s

A dockerized version is available at `jkobject:pycelligner`

to install the latest unstaged version of Celligner in dev mode, do:

```bash
git clone https://github.com/broadinstitute/celligner.git
cd celligner
pip install -e .
```

## For developers

see `CONTRIBUTING.md`

## Using Celligner

The python version has `fit()` and `transform()` functions in the style of scikit-learn models.

A reference expression dataset (e.g. CCLE cell lines TPM expression) should first be fit:

```python
from celligner import Celligner

my_celligner = Celligner()
my_celligner.fit(CCLE_expression, CCLE_annotation)
```

A target expression dataset (e.g. TCGA+ tumor samples) can then be aligned to this reference using the transform function:

```python
my_celligner.transform(TCGA_expression, TCGA_annotation)
```

The combined transformed expression matrix can be accessed via `my_celligner.combined_output` and metrics for this matrix can be computed with `my_celligner.computeMetricsForOutput()`. There are also functions to save/load a fitted model as a .pkl file and many other data attributes.

### Computational complexity

Depending on the dataset, Celligner can be quite memory hungry.
for TCGA, expect at least _50-60Gb_ of memory being used. You might need a powerfull computer, lots of _swap_ and to increase R's default _maximum allowed memory_.

You can also use the `low_memory=True` option to reduce the memory used by celligner in the memory intensive `PCA` & `cPCA` methods.

## Multidataset alignment

For multidataset alignment, there are two broad options:
- The new dataset can be concatenated with an existing reference or target dataset before re-running Celligner.
- The previously transformed reference and target datasets can be used as the new reference against which a new (third) expression dataset is aligned (shown below).

```python
my_celligner.makeNewReference()
my_celligner.transform(your_tpm, your_annotations)
```

# R Celligner

For the original R version of celligner, please check the R/README.md file here: [https://github.com/broadinstitute.org/celligner/tree/master/R/README.md](https://github.com/broadinstitute.org/celligner/tree/master/R/README.md)

Please use _github issues_ for any problem related to the tool.

---

__Initial Project:__

Allie Warren @awarren

__Maintainer:__

Jérémie Kalfon @jkobject
