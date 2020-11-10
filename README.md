# celligner

celligner is a R package to run the celligner pipeline and incorporate new DepMap data

# Install
``` r
library(devtools)
devtools::install_github("acwarren/celligner")
```

The package can then be loaded by calling
``` r
library(celligner)
```

# Functions

The entire pipeline can be run by calling run_Celligner(). 
The parameters are:
  - cell_line_data_name: if cell_line_taiga = TRUE, then the data.name of the taiga file containing the cell line expression data, 
  if cell_line_taiga=FALSE, then the file path to the local folder containing the cell line expression data. To run the pipeline on
  new DepMap data this is the only parameter that should need to be updated (to refer to the new virtual dataset for the relevant release).
  - cell_line_data_file: if cell_line_taiga = TRUE, then the data.file of the taiga file containing the cell line expression data,
  if cell_line_taiga=FALSE, then the name of the file of cell line expression data. By default uses the virtual dataset data file 'CCLE_expression_full'.
  - cell_line_version: (optional) parameter to specify the version to pull from taiga for the cell line expression data, default set to NULL
  - cell_line_taiga: if TRUE (default) then pulls the cell line expression data from taiga, if FALSE then finds cell line expression data in local folder
  - cell_line_ann_name: if cell_line_taiga = TRUE, then the data.name of the taiga file containing the cell line annotations,
  if cell_line_taiga=FALSE, then the file path to the local folder containing the cell line annotations. By default pulls the arxspan data from taiga.
  - cell_line_ann_file: if cell_line_taiga = TRUE, then the data.file of the taiga file containing the cell line annotations,
  if cell_line_taiga=FALSE, then the name of the file of cell line annotations. If pulling from taiga (default), assumes that the file is the arxspan
  file (could also use virtual release sample_info file), if not then assumes it is a local file containing the columns sampleID, lineage, subtype, and type=='CL'.
  - cell_line_ann_version: (optional) parameter to specify the version to pull from taiga for the cell line annotations, default set to NULL
  - cell_line_ann_taiga: if TRUE (default) then pulls the cell line annotations from taiga, if FALSE then finds cell line annotations in local folder
  - tumor_data_name: if tumor_taiga = TRUE (default), then the data.name of the taiga file containing the tumor expression data,
  if tumor_taiga=FALSE, then the file path to the local folder containing the tumor expression data. By default, pulls the TCGA+ (TCGA, TARGET, & Treehouse data 
  downloaded from xena browser)
  - tumor_data_file: if tumor_taiga = TRUE (default), then the data.file of the taiga file containing the tumor expression data,
  if tumor_taiga=FALSE, then the name of the file the tumor expression data
  - tumor_version: (optional) parameter to specify the version to pull from taiga for the tumor expression data, default set to NULL
  - tumor_taiga: if TRUE (default) then pulls the tumor expression data from taiga, if FALSE then finds tumor expression data in local folder
  - tumor_ann_name: if tumor_taiga = TRUE (default), then the data.name of the taiga file containing the tumor annotations,
  if tumor_taiga=FALSE, then the file path to the local folder containing the tumor annotations
  - tumor_ann_file: if tumor_ann_taiga = TRUE (default), then the data.file of the taiga file containing the tumor annotations,
  if tumor_ann_taiga=FALSE, then the name of the file the tumor annotations. If pulling from taiga, assumes that the file is the already 
  created Celligner info file used in the Celligner manuscript, if not then assumes it is a local file containing the columns 
  sampleID, lineage, subtype, and type=='tumor'.
  - tumor_ann_version: (optional) parameter to specify the version to pull from taiga for the tumor annotations, default set to NULL
  - tumor_ann_taiga: if TRUE (default) then pulls the tumor annotations from taiga, if FALSE then finds tumor annotations in local folder
  - additional_annotations_name: if additional_annotations_taiga = TRUE (default), then the data.name of the taiga file containing the additional annotations,
  if additional_annotations_taiga=FALSE, then the file path to the local folder containing the additional annotations. Used to add more fine-grain subtype annotations
  for the cell lines. If NULL, assumes there are no additional annotations. 
  - additional_annotations_file: if additional_annotations_taiga = TRUE (default), then the data.file of the taiga file containing the additional annotations,
  if additional_annotations_taiga=FALSE, then the name of the file the additional annotations. If null, assumes there are
  no additional annotations. By default pulls the Celligner_info file, used in the Celligner manuscript.
  - additional_annotations_version: (optional) parameter to specify the version to pull from taiga for the tumor annotations, default set to NULL
  additional_annotations_taiga: if TRUE (default) then pulls the tumor annotations from taiga, if FALSE then finds tumor annotations in local folder
  - hgnc_data_name: if hgnc_taiga = TRUE (default), then the data.name of the taiga file containing the HGNC gene annotations,
  if hgnc_taiga=FALSE, then the file path to the local folder containing the HGNC gene annotations
  - hgnc_data_file: if hgnc_taiga = TRUE (default), then the data.file of the taiga file containing the HGNC gene annotations,
  if hgnc_taiga=FALSE, then the name of the file the HGNC gne annotations
  - hgnc_version: (optional) parameter to specify the version to pull from taiga for the HGNC gene annotations, default set to NULL
  - hgnc_taiga: if TRUE (default) then pulls the HGNC gene annotations from taiga, if FALSE then finds HGNC gene annotations in local folder
  - save_output: by default is NULL and won't save output, to save output pass in a filepath of where to save the output

It returns:
  - a seurat object containing the Celligner-aligned data, the UMAP dimensionality reduction of the Celligner-aligned data, clustering, sample metadata
  - plots
  - if save_output is given a filepath then output files will be saved to the folder specified
