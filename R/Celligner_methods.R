library(magrittr)
library(tidyverse)


#' @param cell_line_data_name: if cell_line_taiga = TRUE, then the data.name of the taiga file containing the cell line expression data,
#' if cell_line_taiga=FALSE, then the file path to the local folder containing the cell line expression data
#' @param cell_line_data_file: if cell_line_taiga = TRUE, then the data.file of the taiga file containing the cell line expression data,
#' if cell_line_taiga=FALSE, then the name of the file of cell line expression data
#' @param cell_line_version: parameter to specify the version to pull from taiga for the cell line expression data, default set to NULL
#' @param cell_line_taiga: if TRUE then pulls the cell line expression data from taiga, if FALSE then finds cell line expression data in local folder
#' @param cell_line_ann_name: if cell_line_taiga = TRUE, then the data.name of the taiga file containing the cell line annotations,
#' if cell_line_taiga=FALSE, then the file path to the local folder containing the cell line annotations
#' @param cell_line_ann_file: if cell_line_taiga = TRUE, then the data.file of the taiga file containing the cell line annotations,
#' if cell_line_taiga=FALSE, then the name of the file of cell line annotations. If pulling from taiga, assumes that the file is the arxspan
#' file (could also use virtual release sample_info file), if not then assumes it is a local file containing the columns sampleID, lineage, subtype, and type=='CL'.
#' @param cell_line_ann_version: parameter to specify the version to pull from taiga for the cell line annotations, default set to NULL
#' @param cell_line_ann_taiga: if TRUE then pulls the cell line annotations from taiga, if FALSE then finds cell line annotations in local folder
#' @param tumor_data_name: if tumor_taiga = TRUE, then the data.name of the taiga file containing the tumor expression data,
#' if tumor_taiga=FALSE, then the file path to the local folder containing the tumor expression data.
#' If pulling from taiga, assumes that the file is the already create Celligner info file used in the Celligner manuscript,
#'  if not then assumes it is a local file containing the columns sampleID, lineage, subtype, and type=='tumor'.
#' @param tumor_data_file: if tumor_taiga = TRUE, then the data.file of the taiga file containing the tumor expression data,
#' if tumor_taiga=FALSE, then the name of the file the tumor expression data
#' @param tumor_version: parameter to specify the version to pull from taiga for the tumor expression data, default set to NULL
#' @param tumor_taiga: if TRUE then pulls the tumor expression data from taiga, if FALSE then finds tumor expression data in local folder
#' @param tumor_ann_name: if tumor_taiga = TRUE, then the data.name of the taiga file containing the tumor annotations,
#' if tumor_taiga=FALSE, then the file path to the local folder containing the tumor annotations
#' @param tumor_ann_file: if tumor_ann_taiga = TRUE, then the data.file of the taiga file containing the tumor annotations,
#' if tumor_ann_taiga=FALSE, then the name of the file the tumor annotations
#' @param tumor_version: parameter to specify the version to pull from taiga for the tumor annotations, default set to NULL
#' @param tumor_ann_taiga: if TRUE then pulls the tumor annotations from taiga, if FALSE then finds tumor annotations in local folder
#' @param additional_annotations_name: if additional_annotations_taiga = TRUE, then the data.name of the taiga file containing the additional annotations,
#' if additional_annotations_taiga=FALSE, then the file path to the local folder containing the additional annotations. Used to add more fine-grain subtype annotations
#' for the cell lines. If null, assumes there are no additional annotations.
#' @param additional_annotations_file: if additional_annotations_taiga = TRUE, then the data.file of the taiga file containing the additional annotations,
#' if additional_annotations_taiga=FALSE, then the name of the file the additional annotations. If null, assumes there are
#' no additional annotations.
#' @param additional_annotations_version: parameter to specify the version to pull from taiga for the tumor annotations, default set to NULL
#' @param additional_annotations_taiga: if TRUE then pulls the tumor annotations from taiga, if FALSE then finds tumor annotations in local folder
#' @param hgnc_data_name: if hgnc_taiga = TRUE, then the data.name of the taiga file containing the HGNC gene annotations,
#' if hgnc_taiga=FALSE, then the file path to the local folder containing the HGNC gene annotations
#' @param hgnc_data_file: if hgnc_taiga = TRUE, then the data.file of the taiga file containing the HGNC gene annotations,
#' if hgnc_taiga=FALSE, then the name of the file the HGNC gne annotations
#' @param hgnc_version: parameter to specify the version to pull from taiga for the HGNC gene annotations, default set to NULL
#' @param hgnc_taiga: if TRUE then pulls the HGNC gene annotations from taiga, if FALSE then finds HGNC gene annotations in local folder
#'
#' @importFrom magrittr "%>%"
#'
#' @description load expression and annotation files for cell lines and tumors
#' @export load_data
load_data <- function(cell_line_data_name, cell_line_data_file, cell_line_version, cell_line_taiga,
                      cell_line_ann_name, cell_line_ann_file,cell_line_ann_version, cell_line_ann_taiga,
                      tumor_data_name, tumor_data_file, tumor_version, tumor_taiga,
                      tumor_ann_name, tumor_ann_file, tumor_ann_version, tumor_ann_taiga,
                      additional_annotations_name, additional_annotations_file, additional_annotations_version, additional_annotations_taiga,
                      hgnc_data_name, hgnc_data_file, hgnc_version, hgnc_taiga) {

  if(hgnc_taiga) {
    hgnc.complete.set <- taigr::load.from.taiga(data.name = hgnc_data_name, data.version = hgnc_version, data.file = hgnc_data_file)
    if(is.null(hgnc.complete.set)) {
      stop("HGNC gene file input does not exist on taiga")
    }
  } else {
    if(file.exists(file.path(hgnc_data_name, hgnc_data_file))) {
      hgnc.complete.set <- data.table::fread(file.path(hgnc_data_name, hgnc_data_file)) %>%
        as.data.frame()
    } else {
      stop('HGNC gene file input does not exist')
    }
  }

  if(!all(c('symbol', 'ensembl_gene_id', 'locus_group') %in% colnames(hgnc.complete.set))) {
    stop('HGNC gene file does not contain expected columns (symbol, ensembl_gene_id, & locus_group)')
  }

  if(tumor_taiga) {
    TCGA_mat <- taigr::load.from.taiga(data.name = tumor_data_name, data.version = tumor_version, data.file = tumor_data_file)
    if(is.null(TCGA_mat)) {
      stop("tumor expression data file input does not exist on taiga")
    }
  } else {
    if(file.exists(file.path(tumor_data_name, tumor_data_file))) {
      TCGA_mat <-  readr::read_tsv(file.path(tumor_data_name, tumor_data_file)) %>%
        as.data.frame() %>%
        tibble::column_to_rownames('Gene') %>%
        as.matrix() %>%
        t()
    } else {
      stop('tumor expression data file input does not exist')
    }
  }


  if(cell_line_taiga) {
    CCLE_mat <- taigr::load.from.taiga(data.name = cell_line_data_name, data.version = cell_line_version, data.file = cell_line_data_file)
    if(is.null(CCLE_mat)) {
      stop("cell line expression data file input does not exist on taiga")
    }
    } else {
      if(file.exists(file.path(cell_line_data_name, cell_line_data_file))) {
        CCLE_mat <-  readr::read_csv(file.path(cell_line_data_name, cell_line_data_file)) %>%
          as.data.frame() %>%
          tibble::column_to_rownames('X1') %>%
          as.matrix()
      } else {
        stop('cell line data file input does not exist')
      }
  }

  # subset gene names to just ensembl IDs
  # add test for this
  colnames(CCLE_mat) <- stringr::str_match(colnames(CCLE_mat), '\\((.+)\\)')[,2]

  # convert tumor gene names to ensembl IDs, if needed
  if(length(grep('ENS', colnames(TCGA_mat))) != ncol(TCGA_mat)) {
    print('converting TCGA column names from HGNC ids to ensembl ids')
    common_genes <- dplyr::intersect(colnames(TCGA_mat), hgnc.complete.set$symbol)
    if(length(common_genes) < 10000) {
      sprint('only %s genes in overlapping between genes in columns of the tumor data and hgnc dataset')
      warning('low overlap of genes in tumor data and gene symbol, either tumor data
              or gene file may not be in correct format')
    }
    TCGA_mat <- TCGA_mat[,common_genes]
    hgnc.complete.set <- dplyr::filter(hgnc.complete.set, symbol %in% common_genes)
    hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$symbol),]
    rownames(hgnc.complete.set) <- hgnc.complete.set$symbol
    hgnc.complete.set <- hgnc.complete.set[common_genes,]
    colnames(TCGA_mat) <- hgnc.complete.set$ensembl_gene_id
  }


  if(cell_line_ann_taiga) {
    CCLE_ann <- taigr::load.from.taiga(data.name = cell_line_ann_name, data.version = cell_line_ann_version, data.file = cell_line_ann_file)
    column_names <- c('arxspan_id', 'lineage', 'lineage_subtype')
    if('DepMap_ID' %in% colnames(CCLE_ann)) {
      column_names[1] <- 'DepMap_ID'
    }
    if(is.null(CCLE_ann)) {
      warning('cell line annotation file does not exist on taiga, creating default annotations')
      CCLE_ann <- data.frame(sampleID =  rownames(CCLE_mat),
                             lineage = NA,
                             subtype = NA,
                             type = 'CL')
    }
    if(!all(column_names %in% colnames(CCLE_ann))) {
      warning('cell line annotation file does not contain expected columns (arxspan_id or DepMap_ID, lineage, & lineage_subtype), creating default annotation file')
      CCLE_ann <- data.frame(sampleID =  rownames(CCLE_mat),
                             lineage = NA,
                             subtype = NA,
                             type = 'CL')
    } else {
      CCLE_ann <- CCLE_ann[,column_names]
      colnames(CCLE_ann) <- c('sampleID', 'lineage', 'subtype')
      CCLE_ann$type <- 'CL'
    }
  } else {
    if(file.exists(file.path(cell_line_ann_name, cell_line_ann_file))) {
      CCLE_ann <- data.table::fread(file.path(cell_line_ann_name, cell_line_ann_file)) %>%
        as.data.frame()
    } else {
      warning('cell line annotation file does not exist, creating default annotations')
      CCLE_ann <- data.frame(sampleID =  rownames(CCLE_mat),
                             lineage = NA,
                             subtype = NA,
                             type = 'CL')
    }
  }

  if(!all(c('sampleID', 'lineage', 'subtype', 'type') %in% colnames(CCLE_ann))) {
    warning('cell line annotation file does not contain expected columns (sampleID, lineage, subtype & type), creating default annotations')
    CCLE_ann <- data.frame(sampleID =  rownames(CCLE_mat),
                           lineage = NA,
                           subtype = NA,
                           type = 'CL')
  }

  if(tumor_ann_taiga) {
    TCGA_ann <- taigr::load.from.taiga(data.name = tumor_ann_name, data.version = tumor_ann_version, data.file = tumor_ann_file)
    tumor_column_names <- c('sampleID', 'lineage', 'subtype')
    if(is.null(TCGA_ann)) {
      warning('tumor annotation file does not exist on taiga, creating default annotations')
      TCGA_ann <- data.frame(sampleID =  rownames(TCGA_mat),
                             lineage = NA,
                             subtype = NA,
                             type = 'tumor')
    }
    if(!all(tumor_column_names %in% colnames(TCGA_ann))) {
      warning('tumor annotation file does not contain expected columns (sampleID, lineage, & subtype), creating default tumor annotations')
      TCGA_ann <- data.frame(sampleID =  rownames(TCGA_mat),
                             lineage = NA,
                             subtype = NA,
                             type = 'tumor')
    } else {
      TCGA_ann <- TCGA_ann[,tumor_column_names]
      TCGA_ann$type <- 'tumor'
    }
  } else {
    if(file.exists(file.path(tumor_ann_name, tumor_ann_file))) {
      TCGA_ann <- data.table::fread(file.path(tumor_ann_name, tumor_ann_file)) %>%
        as.data.frame()
    } else {
      warning('tumor annotation file does not exist, creating default annotations')
      TCGA_ann <- data.frame(sampleID =  rownames(TCGA_mat),
                             lineage = NA,
                             subtype = NA,
                             type = 'tumor')
    }
    if(!all(c('sampleID', 'lineage', 'subtype', 'type') %in% colnames(TCGA_ann))) {
      warning('tumor annotation file does not contain expected columns (sampleID, lineage, subtype & type), creating default annotations')
      TCGA_ann <- data.frame(sampleID =  rownames(TCGA_mat),
                             lineage = NA,
                             subtype = NA,
                             type = 'tumor')
    }
  }

  if(!(is.null(additional_annotations_name) | is.null(additional_annotations_file))) {
    if(additional_annotations_taiga) {
      add_ann <- taigr::load.from.taiga(data.name = additional_annotations_name, data.version = additional_annotations_version, data.file = additional_annotations_file)
      tumor_column_names <- c('sampleID', 'lineage', 'subtype', 'type')
      if(is.null(add_ann)) {
        warning('additional annotation file does not exist on taiga, no additional annotations used')
      }
      if(!all(c('sampleID', 'subtype') %in% colnames(add_ann))) {
        warning('additional annotation file does not contain expected columns (sampleID & subtype), no additional annotations used')
      } else {
        shared_samples <- intersect(CCLE_ann$sampleID, add_ann$sampleID)
        CCLE_ann[match(shared_samples, CCLE_ann$sampleID),'subtype'] <- add_ann[match(shared_samples, add_ann$sampleID),'subtype']
      }
    } else {
    if(file.exists(file.path(additional_annotations_name, additional_annotations_file))) {
      add_ann <- data.table::fread(file.path(additional_annotations_name, additional_annotations_file)) %>%
        as.data.frame()
      if(!all(c('sampleID', 'subtype') %in% colnames(add_ann))) {
        warning('additional annotation file does not contain expected columns (sampleID & subtype), no additional annotations used')
      } else {
        shared_samples <- intersect(CCLE_ann$sampleID, add_ann$sampleID)
        CCLE_ann[match(shared_samples, CCLE_ann$sampleID),'subtype'] <- add_ann[match(shared_samples, add_ann$sampleID),'subtype']
      }
      } else {
        warning('additional annotation file does not exist, no additional annotations used')

      }
    }
  }


  # subset to samples in both the annotation and gene expression matrices, and match ordering between them
  common_cls <- intersect(rownames(CCLE_mat), CCLE_ann$sampleID)
  if(length(setdiff(rownames(CCLE_mat), CCLE_ann$sampleID))>0) {
    sprintf('Missing annotations for these cell lines: %s', paste(rownames(CCLE_mat), CCLE_ann$sampleID, collapse=", "))
  }

  CCLE_mat <- CCLE_mat[common_cls,]
  CCLE_ann <- CCLE_ann[match(common_cls,CCLE_ann$sampleID),]

  common_tumors <- intersect(rownames(TCGA_mat), TCGA_ann$sampleID)
  if(length(setdiff(rownames(TCGA_mat), common_tumors))>0) {
    sprintf('Missing annotations for these tumors: %s', paste(rownames(TCGA_mat), common_tumors, collapse=", "))
  }
  TCGA_mat <- TCGA_mat[common_tumors,]
  TCGA_ann <- TCGA_ann[match(common_tumors,TCGA_ann$sampleID),]

  # subset genes to functional genes
  func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id
  genes_used <- intersect(colnames(TCGA_mat), colnames(CCLE_mat))
  genes_used <- intersect(genes_used, func_genes)

  TCGA_mat <- TCGA_mat[,genes_used]
  CCLE_mat <- CCLE_mat[,genes_used]

  return(list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, CCLE_mat = CCLE_mat, CCLE_ann = CCLE_ann))
}

#' @param dat: data object containing tumor and cell line expression data and annotations produced by running load_data
#' @param hgnc_data_name: if hgnc_taiga = TRUE, then the data.name of the taiga file containing the HGNC gene annotations,
#' if hgnc_taiga=FALSE, then the file path to the local folder containing the HGNC gene annotations
#' @param hgnc_data_file: if hgnc_taiga = TRUE, then the data.file of the taiga file containing the HGNC gene annotations,
#' if hgnc_taiga=FALSE, then the name of the file the HGNC gne annotations
#' @param hgnc_version: parameter to specify the version to pull from taiga for the HGNC gene annotations, default set to NULL
#' @param hgnc_taiga: if TRUE then pulls the HGNC gene annotations from taiga, if FALSE then finds HGNC gene annotations in local folder
#'
#' @description calculate the average gene expression and variance
#' @export cluster_data
calc_gene_stats <- function(dat, hgnc_data_name, hgnc_data_file, hgnc_version, hgnc_taiga) {
  common_genes <- intersect(colnames(dat$TCGA_mat), colnames(dat$CCLE_mat))

  if(hgnc_taiga) {
    hgnc.complete.set <- taigr::load.from.taiga(data.name = hgnc_data_name, data.version = hgnc_version, data.file = hgnc_data_file)
    if(is.null(hgnc.complete.set)) {
      stop("HGNC gene file input does not exist on taiga")
    }
  } else {
    if(file.exists(file.path(hgnc_data_name, hgnc_data_file))) {
      hgnc.complete.set <- data.table::fread(file.path(hgnc_data_name, hgnc_data_file)) %>%
        as.data.frame()
    } else {
      stop('HGNC gene file input does not exist')
    }
  }

  if(!all(c('symbol', 'ensembl_gene_id', 'locus_group') %in% colnames(hgnc.complete.set))) {
    stop('HGNC gene file does not contain expected columns (symbol, ensembl_gene_id, & locus_group)')
  }

  hgnc.complete.set <- hgnc.complete.set %>%
    dplyr::select(Gene = ensembl_gene_id, Symbol = symbol) %>%
    filter(Gene %in% common_genes)
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$Gene),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$Gene
  hgnc.complete.set <- hgnc.complete.set[common_genes,]

  gene_stats <- data.frame(
    Tumor_SD = apply(dat$TCGA_mat, 2, sd, na.rm=T),
    CCLE_SD = apply(dat$CCLE_mat, 2, sd, na.rm=T),
    Tumor_mean = colMeans(dat$TCGA_mat, na.rm=T),
    CCLE_mean = colMeans(dat$CCLE_mat, na.rm=T),
    Gene = common_genes,
    stringsAsFactors = F) %>%
  dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) #add avg and max SD per gene

  gene_stats <- left_join(hgnc.complete.set, gene_stats, by = "Gene")

  return(gene_stats)

  }


#' @param exp_mat: matrix of samples by genes, where genes are ensembl gene IDs. Data should be log2(X+1) TPM data.
#' @param ann: matrix of sample anntoations. Expects column 'sampleID' which matches the rownames of exp_mat.
#' @param type: optional parameter, string specifying the data type of the current data (ex. 'tumor'), which is added to the annotation matrix.
#' @description create Seurat object of expression data and annotations and run dimensionality reduction.
#' Dimensionality reductions will be run with the parameters (n_PC_dims, umap_n_neighbors, umap_min_dist, distance_metric) specified in global.
#' @export cluster_data
create_Seurat_object <- function(exp_mat, ann, type = NULL) {
  seu_obj <- Seurat::CreateSeuratObject(t(exp_mat),
                                         min.cells = 0,
                                         min.features = 0,
                                         meta.data = ann %>%
                                           magrittr::set_rownames(ann$sampleID))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)

  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                               features = rownames(Seurat::GetAssayData(seu_obj)),
                               npcs = global$n_PC_dims, verbose = F)

  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                reduction = 'pca',
                                n.neighbors = global$umap_n_neighbors,
                                min.dist =  global$umap_min_dist,
                                metric = global$distance_metric, verbose=F)

  return(seu_obj)
}

#' @param seu_obj: seurat object containing expression data and sample annotations.
#' Expects PCA for the seurat object has already been calculated.
#' @description cluster data in seurat object, using default Seurat clustering method. Clsuters data
#' within PCA space using the number of dimensions provided in global$n_PC_dims (default is 70)
#'
#' @export cluster_data
cluster_data <- function(seu_obj) {
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = 'pca',
                                    dims = 1:global$n_PC_dims,
                                    k.param = 20,
                                    force.recalc = TRUE,
                                    verbose = FALSE)

  seu_obj %<>% Seurat::FindClusters(reduction = 'pca',
                                     resolution = global$mod_clust_res)

  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters

  return(seu_obj)

  }

#' @param seu_obj: seurat object containing expression data and sample annotations. Expects data in the Seurat object
#' slot scale.data and a column 'seurat_clusters' within the meta.data of the Seurat object.
#' @description find genes that are differentially expressed between clusters within the expression data
#'
#' @export find_differentially_expressed_genes
find_differentially_expressed_genes <- function(seu_obj) {
  if(nrow(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data'))==0) {
    stop("Seurat object doesn't have expression data at scale.data, run 'create_Seurat_object' first")
  }
  if(!'seurat_clusters' %in% colnames(seu_obj@meta.data)) {
    stop("Seurat object doesn't contain the column 'seurat_clusters', run 'cluster_data' first")
  }
  n_clusts <- nlevels(seu_obj@meta.data$seurat_clusters)
  if (n_clusts > 2) {
    cur_DE_genes <- run_lm_stats_limma_group(
      t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
      seu_obj@meta.data %>% dplyr::select(seurat_clusters),
      limma_trend = TRUE) %>%
      dplyr::select(Gene, gene_stat = F_stat)
  } else if (n_clusts == 2) {
    cur_DE_genes <- run_lm_stats_limma(t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
                                               seu_obj@meta.data$cluster,
                                               limma_trend = TRUE) %>%
      dplyr::mutate(gene_stat = abs(t_stat)) %>%
      dplyr::select(Gene, gene_stat)
  } else {
    cur_DE_genes <- data.frame(Gene = colnames(seu_obj), gene_stat = NA)
  }

  return(cur_DE_genes)

}

#' @param TCGA_obj: seurat object containing expression data and sample annotations, usually the tumor data
#' @param CCLE_obj: seurat object containing expression data and sample annotations, usually the cell line data
#' @param pc_dims: the number of cPCs calculated. If set to null then all cPCs will be calculated (this is quite slow), but if set to
#' some value >=4 then an approximate cPCA will be calculated, which just calculates the input number of contrastive principle components,
#' which is quicker.
#' @description run contrastive principal components analysis.
#' Set pc_dims to a value >= 4 to run fast cPCA by just calculating the top contrastive principle components
#'
#' @export run_cPCA
run_cPCA <- function(TCGA_obj, CCLE_obj, pc_dims = NULL) {
  if(nrow(Seurat::GetAssayData(TCGA_obj, assay='RNA', slot='scale.data'))==0) {
    stop("TCGA seurat object doesn't have expression data at scale.data, run 'create_Seurat_object' first")
  }
  if(nrow(Seurat::GetAssayData(CCLE_obj, assay='RNA', slot='scale.data'))==0) {
    stop("CCLE seurat object doesn't have expression data at scale.data, run 'create_Seurat_object' first")
  }
  cov_diff_eig <- run_cPCA_analysis(t(Seurat::GetAssayData(TCGA_obj, assay='RNA', slot='scale.data')),
                                    t(Seurat::GetAssayData(CCLE_obj, assay='RNA', slot='scale.data')),
                                    TCGA_obj@meta.data, CCLE_obj@meta.data, pc_dims=pc_dims)
 return(cov_diff_eig)
}

#' @param CCLE_cor: matrix of samples by genes of cPC corrected data that serves as the reference data in the MNN alignment.
#' In the default Celligner pipeline this the cell line data.
#' @param TCGA_cor: matrix of samples by genes of cPC corrected data that is corrected in the MNN alignment and projected onto the reference data.
#' In the default Celligner pipeline this the tumor data.
#' @param k1: the number of neighbors within the data being corrected (by default the tumor data). By default this
#' pulls from the global paramter mnn_k_tumor, which by default is 50.
#' @param k2: the number of neighbors within the reference data (by default the cell line data). By default this
#' pulls from the global parameter mnn_k_CL, which by default is 5.
#' @param ndist: A numeric scalar specifying the threshold beyond which neighbors are to be ingnored when computing correction vectors.
#' By default this pulls from the global parameter mnn_ndist, which by default is 3.
#' @param subset_genes: the subset of genes used for identifying mutual nearest neighbors within the datasets. The set of differentially
#' expressed genes is usually passed here.
#' @description run MNN batch correction to align data to a reference dataset
#'
#' @export run_MNN
run_MNN <- function(CCLE_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist,
                    subset_genes) {
  mnn_res <- modified_mnnCorrect(CCLE_cor, TCGA_cor, k1 = k1, k2 = k2, ndist = ndist,
                             subset_genes = subset_genes)

  return(mnn_res)
}

#' @param Celligner_aligned_data: Celligner aligned data matrix of samples (cells line and tumors) by genes
#' @param Celligner_info: annotation file of cell line and tumor samples with a column 'type' marking samples as either
#' cell lines or tumors and a column 'sampleID' that matches the row names of Celligner_aligned_data
#' @description calculate the correlation between cell line and tumor samples in the Celligner aligned data
#'
#' @export calc_tumor_CL_cor
calc_tumor_CL_cor <- function(Celligner_aligned_data, Celligner_info) {
  tumors_samples <- dplyr::filter(Celligner_info, type=='tumor')$sampleID
  cl_samples <- dplyr::filter(Celligner_info, type=='CL')$sampleID
  tumor_CL_cor <- cor(t(Celligner_aligned_data[tumor_samples,]), t(Celligner_aligned_data[cl_samples,]),
                      use='pairwise')


  return(tumor_CL_cor)
}


#' @param cell_line_data_name: if cell_line_taiga = TRUE, then the data.name of the taiga file containing the cell line expression data,
#' if cell_line_taiga=FALSE, then the file path to the local folder containing the cell line expression data
#' @param cell_line_data_file: if cell_line_taiga = TRUE, then the data.file of the taiga file containing the cell line expression data,
#' if cell_line_taiga=FALSE, then the name of the file of cell line expression data
#' @param cell_line_version: parameter to specify the version to pull from taiga for the cell line expression data, default set to NULL
#' @param cell_line_taiga: if TRUE then pulls the cell line expression data from taiga, if FALSE then finds cell line expression data in local folder
#' @param cell_line_ann_name: if cell_line_taiga = TRUE, then the data.name of the taiga file containing the cell line annotations,
#' if cell_line_taiga=FALSE, then the file path to the local folder containing the cell line annotations
#' @param cell_line_ann_file: if cell_line_taiga = TRUE, then the data.file of the taiga file containing the cell line annotations,
#' if cell_line_taiga=FALSE, then the name of the file of cell line annotations. If pulling from taiga, assumes that the file is the arxspan
#' file (could also use virtual release sample_info file), if not then assumes it is a local file containing the columns sampleID, lineage, subtype, and type=='CL'.
#' @param cell_line_ann_version: parameter to specify the version to pull from taiga for the cell line annotations, default set to NULL
#' @param cell_line_ann_taiga: if TRUE then pulls the cell line annotations from taiga, if FALSE then finds cell line annotations in local folder
#' @param tumor_data_name: if tumor_taiga = TRUE, then the data.name of the taiga file containing the tumor expression data,
#' if tumor_taiga=FALSE, then the file path to the local folder containing the tumor expression data.
#' If pulling from taiga, assumes that the file is the already create Celligner info file used in the Celligner manuscript,
#'  if not then assumes it is a local file containing the columns sampleID, lineage, subtype, and type=='tumor'.
#' @param tumor_data_file: if tumor_taiga = TRUE, then the data.file of the taiga file containing the tumor expression data,
#' if tumor_taiga=FALSE, then the name of the file the tumor expression data
#' @param tumor_version: parameter to specify the version to pull from taiga for the tumor expression data, default set to NULL
#' @param tumor_taiga: if TRUE then pulls the tumor expression data from taiga, if FALSE then finds tumor expression data in local folder
#' @param tumor_ann_name: if tumor_taiga = TRUE, then the data.name of the taiga file containing the tumor annotations,
#' if tumor_taiga=FALSE, then the file path to the local folder containing the tumor annotations
#' @param tumor_ann_file: if tumor_ann_taiga = TRUE, then the data.file of the taiga file containing the tumor annotations,
#' if tumor_ann_taiga=FALSE, then the name of the file the tumor annotations
#' @param tumor_version: parameter to specify the version to pull from taiga for the tumor annotations, default set to NULL
#' @param tumor_ann_taiga: if TRUE then pulls the tumor annotations from taiga, if FALSE then finds tumor annotations in local folder
#' @param additional_annotations_name: if additional_annotations_taiga = TRUE, then the data.name of the taiga file containing the additional annotations,
#' if additional_annotations_taiga=FALSE, then the file path to the local folder containing the additional annotations. Used to add more fine-grain subtype annotations
#' for the cell lines. If null, assumes there are no additional annotations.
#' @param additional_annotations_file: if additional_annotations_taiga = TRUE, then the data.file of the taiga file containing the additional annotations,
#' if additional_annotations_taiga=FALSE, then the name of the file the additional annotations. If null, assumes there are
#' no additional annotations.
#' @param additional_annotations_version: parameter to specify the version to pull from taiga for the tumor annotations, default set to NULL
#' @param additional_annotations_taiga: if TRUE then pulls the tumor annotations from taiga, if FALSE then finds tumor annotations in local folder
#' @param hgnc_data_name: if hgnc_taiga = TRUE, then the data.name of the taiga file containing the HGNC gene annotations,
#' if hgnc_taiga=FALSE, then the file path to the local folder containing the HGNC gene annotations
#' @param hgnc_data_file: if hgnc_taiga = TRUE, then the data.file of the taiga file containing the HGNC gene annotations,
#' if hgnc_taiga=FALSE, then the name of the file the HGNC gne annotations
#' @param hgnc_version: parameter to specify the version to pull from taiga for the HGNC gene annotations, default set to NULL
#' @param hgnc_taiga: if TRUE then pulls the HGNC gene annotations from taiga, if FALSE then finds HGNC gene annotations in local folder
#' @param save_output: by default is NULL and won't save output, to save output pass in a filepath of where to save the output
#'
#' @importFrom magrittr "%>%"
#'
#' @description run all parts of the Celligner pipeline
#'
#' @export run_Celligner
run_Celligner <- function(cell_line_data_name='public-20q4-a4b3', cell_line_data_file = 'CCLE_expression_full', cell_line_version = NULL, cell_line_taiga=TRUE,
                          cell_line_ann_name='arxspan-cell-line-export-f808', cell_line_ann_file = 'ACH',cell_line_ann_version = NULL, cell_line_ann_taiga=TRUE,
                          tumor_data_name = 'celligner-input-9827', tumor_data_file = 'tumor_expression', tumor_version = NULL, tumor_taiga = TRUE,
                          tumor_ann_name = 'celligner-input-9827', tumor_ann_file = 'tumor_annotations', tumor_ann_version = NULL, tumor_ann_taiga = TRUE,
                          additional_annotations_name = 'celligner-input-9827', additional_annotations_file = 'CCLE_annotations', additional_annotations_version = NULL, additional_annotations_taiga = TRUE,
                          hgnc_data_name = 'hgnc-87ab', hgnc_data_file='hgnc_complete_set', hgnc_version= NULL, hgnc_taiga = TRUE,
                          save_output = NULL) {

  require(magrittr)
  require(tidyverse)

  dat <- load_data(cell_line_data_name, cell_line_data_file, cell_line_version, cell_line_taiga,
                   cell_line_ann_name, cell_line_ann_file,cell_line_ann_version, cell_line_ann_taiga,
                   tumor_data_name, tumor_data_file, tumor_version, tumor_taiga,
                   tumor_ann_name, tumor_ann_file, tumor_ann_version, tumor_ann_taiga,
                   additional_annotations_name, additional_annotations_file, additional_annotations_version, additional_annotations_taiga,
                   hgnc_data_name, hgnc_data_file, hgnc_version, hgnc_taiga)

   gene_stats <- calc_gene_stats(dat, hgnc_data_name, hgnc_data_file, hgnc_version, hgnc_taiga)

  comb_ann <- rbind(
    dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype) %>%
      dplyr::mutate(type = 'tumor'),
    dat$CCLE_ann %>% dplyr::select(sampleID, lineage, subtype) %>%
      dplyr::mutate(type = 'CL')
  )

  TCGA_obj <- create_Seurat_object(dat$TCGA_mat, dat$TCGA_ann, type='tumor')
  CCLE_obj <- create_Seurat_object(dat$CCLE_mat, dat$CCLE_ann, type='CL')

  TCGA_obj <- cluster_data(TCGA_obj)
  CCLE_obj <- cluster_data(CCLE_obj)

  tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj)
  CL_DE_genes <- find_differentially_expressed_genes(CCLE_obj)

  DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
    mutate(
      tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
      CL_rank = dplyr::dense_rank(-gene_stat_CL),
      best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
    dplyr::left_join(gene_stats, by = 'Gene')

  # take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
  DE_gene_set <- DE_genes %>%
    dplyr::filter(best_rank < global$top_DE_genes_per) %>%
    .[['Gene']]


  cov_diff_eig <- run_cPCA(TCGA_obj, CCLE_obj, global$fast_cPCA)

  if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, global$remove_cPCA_dims, drop = FALSE]
  } else {
    cur_vecs <- cov_diff_eig$rotation[, global$remove_cPCA_dims, drop = FALSE]
  }

  rownames(cur_vecs) <- colnames(dat$TCGA_mat)
  TCGA_cor <- resid(lm(t(dat$TCGA_mat) ~ 0 + cur_vecs)) %>% t()
  CCLE_cor <- resid(lm(t(dat$CCLE_mat) ~ 0 + cur_vecs)) %>% t()

  mnn_res <- run_MNN(CCLE_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist,
                      subset_genes = DE_gene_set)

  combined_mat <- rbind(mnn_res$corrected, CCLE_cor)

  comb_obj <- create_Seurat_object(combined_mat, comb_ann)
  comb_obj <- cluster_data(comb_obj)

  if(!is.null(save_output)) {
    if(file.exists(save_output)) {
      print('calculating tumor/cell line correlation')
      tumor_CL_cor <- calc_tumor_CL_cor(combined_mat, comb_obj@meta.data)

      Celligner_ouput <- Seurat::Embeddings(comb_obj, reduction = 'umap') %>%
        as.data.frame() %>%
        set_colnames(c('UMAP_1', 'UMAP_2')) %>%
        rownames_to_column(var = 'sampleID') %>%
        left_join(comb_obj@meta.data, by = 'sampleID')

      print('saving files')
      write.csv(tumor_CL_cor, file.path(save_output, 'tumor_CL_cor.csv'))
      write.csv(combined_mat, file.path(save_output, 'Celligner_aligned_data.csv'))
      write_csv(Celligner_output, file.path(save_output, 'Celligner_info.csv'))


    } else{
      warning("can't save output, folder does not exist")
    }
  }

  return(comb_obj)
}







