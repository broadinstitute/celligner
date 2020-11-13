
#' Load additional expression and annotation data
#' @name load_additional_data
#'
#' @param data_name: if data_taiga = TRUE, then the data.name of the taiga file containing the expression data,
#' if data_taiga=FALSE, then the file path to the local folder containing the expression data. Assumes that genes
#' are labeled using ensembl IDs and that there are fewer samples than genes in the matrix, will transpose the matrix
#' so that rows are samples and columns are genes.
#' @param data_file: if data_taiga = TRUE, then the data.file of the taiga file containing the expression data,
#' if data_taiga = FALSE, then the name of the file of expression data
#' @param data_version: (optional) parameter to specify the version to pull from taiga for the expression data, default set to NULL
#' @param data_taiga: if TRUE then pulls the expression data from taiga, if FALSE then finds expression data in local folder
#' @param ann_name: if ann_taiga = TRUE, then the data.name of the taiga file containing the data annotations,
#' if ann_taiga=FALSE, then the file path to the local folder containing the annotations
#' @param ann_file: if ann_taiga = TRUE, then the data.file of the taiga file containing the data annotations,
#' if ann_taiga=FALSE, then the name of the file of data annotations
#' @param ann_version: (optional) parameter to specify the version to pull from taiga for the annotations, default set to NULL
#' @param ann_taiga: if TRUE (default) then pulls the annotations from taiga, if FALSE then finds cell line annotations in local folder
#' @param data_type: string added to the annotation file under the column type to specify the data, default is ""
#' @description load additional expression and annotation files
#'
#' @return object containing expression matrix and annotations table
#' @export
#'
load_additional_data <- function(data_name, data_file, data_version = NULL, data_taiga = TRUE,
                     ann_name, ann_file, ann_version = NULL, ann_taiga = TRUE, data_type = "") {

  if(data_taiga) {
    data_mat <- taigr::load.from.taiga(data.name = data_name, data.version = data_version, data.file = data_file)
    if(is.null(data_mat)) {
      stop("expression data file input does not exist on taiga")
    }
  } else {
    if(file.exists(file.path(data_name, data_file))) {
      data_mat <-  readr::read_csv(file.path(data_name, data_file)) %>%
        as.data.frame() %>%
        tibble::column_to_rownames('X1') %>%
        as.matrix()
    } else {
      stop('expression data file input does not exist')
    }
  }


  # transpose matrix, if needed, so rownames are samples and column names are genes
  if(nrow(data_mat) > ncol(data_mat)) {
    warning('more rows than columns, taking transpose of expression matrix')
    data_mat <- t(data_mat)
  }


  if(ann_taiga) {
    ann <- taigr::load.from.taiga(data.name = ann_name, data.version = ann_version, data.file = ann_file)
    column_names <- c('sampleID', 'lineage', 'subtype')
    if(is.null(ann)) {
      warning('annotation file does not exist on taiga, creating default annotations')
      ann <- data.frame(sampleID =  rownames(data_mat),
                             lineage = NA,
                             subtype = NA,
                             type = data_type)
    }
    if(!all(column_names %in% colnames(ann))) {
      warning('annotation file does not contain expected columns (sampleID, lineage, & subtype), creating tumor annotations')
      ann <- data.frame(sampleID =  rownames(data_mat),
                             lineage = NA,
                             subtype = NA,
                             type = data_type)
    } else {
      ann <- ann[,column_names]
      ann$type <- data_type
    }
  } else {
    if(file.exists(file.path(ann_name, ann_file))) {
      ann <- data.table::fread(file.path(ann_name, ann_file)) %>%
        as.data.frame()
    } else {
      warning('annotation file does not exist, creating default annotations')
      ann <- data.frame(sampleID =  rownames(data_mat),
                             lineage = NA,
                             subtype = NA,
                             type = data_type)
    }
    if(!all(c('sampleID', 'lineage', 'subtype', 'type') %in% colnames(ann))) {
      warning('annotation file does not contain expected columns (sampleID, lineage, subtype & type), creating default annotations')
      ann <- data.frame(sampleID =  rownames(data_mat),
                             lineage = NA,
                             subtype = NA,
                             type = data_type)
    }
  }

  return(list(mat = data_mat, ann = ann))

}



#' All methods to run Celligner, with additional alignment of Met500 and PDX data, and save the output, if desired
#' @name run_multidataset_alignment
#'
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
#' @param cell_line_ann_taiga: if TRUE (default) then pulls the cell line annotations from taiga, if FALSE then finds cell line annotations in local folder
#' @param tumor_data_name: if tumor_taiga = TRUE, then the data.name of the taiga file containing the tumor expression data,
#' if tumor_taiga=FALSE, then the file path to the local folder containing the tumor expression data.
#' @param tumor_data_file: if tumor_taiga = TRUE, then the data.file of the taiga file containing the tumor expression data,
#' if tumor_taiga=FALSE, then the name of the file the tumor expression data
#' @param tumor_version: parameter to specify the version to pull from taiga for the tumor expression data, default set to NULL
#' @param tumor_taiga: if TRUE (default) then pulls the tumor expression data from taiga, if FALSE then finds tumor expression data in local folder
#' @param tumor_ann_name: if tumor_taiga = TRUE, then the data.name of the taiga file containing the tumor annotations,
#' if tumor_taiga=FALSE, then the file path to the local folder containing the tumor annotations
#' @param tumor_ann_file: if tumor_ann_taiga = TRUE, then the data.file of the taiga file containing the tumor annotations,
#' if tumor_ann_taiga=FALSE, then the name of the file the tumor annotations. If pulling from taiga, assumes that the file is the already create Celligner info file used in the Celligner manuscript,
#'  if not then assumes it is a local file containing the columns sampleID, lineage, subtype, and type=='tumor'.
#' @param tumor_ann_version: parameter to specify the version to pull from taiga for the tumor annotations, default set to NULL
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
#' @param met500_data_name: Met500 expression, default pulls from taiga, this the data_name of the taiga dataset, or path to folder if using met500_taiga=F
#' @param met500_data_file: default pulls from taiga, this the data_file of the taiga dataset, or name of local file if using met500_taiga=F
#' @param met500_version: default NULL, used to specify version of taiga dataset
#' @param met500_taiga: if TRUE (default) pulls Met500 expression from taiga dataset, if FALSE reads from local
#' @param met500_ann_name: Met500 annotations, default pulls from taiga, this the data_name of the taiga dataset, or path to folder is using met500_ann_taiga=F
#' @param met500_ann_file: Met500 annotations, default pulls from taiga, this the data_file of the taiga dataset, or name of local file is using met500_ann_taiga=F
#' @param met500_ann_version: default NULL, used to specify version of taiga dataset
#' @param met500_ann_taiga: if TRUE (default) pulls met500 annotations from taiga dataset, if FALSE reads from local
#' @param Novartis_PDX_data_name: Novartis PDX expression, default pulls from taiga, this the data_name of the taiga dataset, or path to folder if using Novartis_PDX_taiga=F
#' @param Novartis_PDX_data_file: default pulls from taiga, this the data_file of the taiga dataset, or name of local file if using Novartis_PDX_taiga=F
#' @param Novartis_PDX_version: default NULL, used to specify version of taiga dataset
#' @param Novartis_PDX_taiga: if TRUE (default) pulls Novartis PDX expression from taiga dataset, if FALSE reads from local
#' @param Novartis_PDX_ann_name: Novartis PDX annotations, default pulls from taiga, this the data_file of the taiga dataset, or path to folder is using met500_ann_taiga=F
#' @param Novartis_PDX_ann_file: Novartis PDX annotations, default pulls from taiga, this the data_name of the taiga dataset, or name of local file is using Novartis_PDX_ann_taiga=F
#' @param Novartis_PDX_ann_version: default NULL, used to specify version of taiga dataset
#' @param Novartis_PDX_ann_taiga: if TRUE (default) pulls Novartis PDX annotations from taiga dataset, if FALSE reads from local
#' @param pediatric_PDX_data_name: pediatric PDX expression, default pulls from taiga, this the data_name of the taiga dataset, or path to folder if using pediatric_PDX_taiga=F
#' @param pediatric_PDX_data_file: default pulls from taiga, this the data_file of the taiga dataset, or name of local file if using pediatric_PDX_taiga=F
#' @param pediatric_PDX_version: default NULL, used to specify version of taiga dataset
#' @param pediatric_PDX_taiga: if TRUE (default) pulls pediatric PDX expression from taiga dataset, if FALSE reads from local
#' @param pediatric_PDX_ann_name: Pediatric PDX annotations, default pulls from taiga, this the data_name of the taiga dataset, or path to folder is using pediatric_PDX_ann_taiga=F
#' @param pediatric_PDX_ann_file: Pediatric PDX annotations, default pulls from taiga, this the data_file of the taiga dataset, or name of local file is using pediatric_PDX_ann_taiga=F
#' @param pediatric_PDX_ann_version: default NULL, used to specify version of taiga dataset
#' @param pediatric_PDX_ann_taiga: if TRUE (default) pulls pediatric PDX annotations from taiga dataset, if FALSE reads from local
#' @param save_output: by default is NULL and won't save output, to save output pass in a filepath of where to save the output
#'
#' @importFrom magrittr "%>%"
#'
#' @description run all parts of the Celligner pipeline, with alignment of additional datasets
#'
#' @return seurat object of the Celligner-aligned data
#' @export
#'
run_multidataset_alignment <- function(cell_line_data_name='public-20q4-a4b3', cell_line_data_file = 'CCLE_expression_full', cell_line_version = NULL, cell_line_taiga=TRUE,
                          cell_line_ann_name='arxspan-cell-line-export-f808', cell_line_ann_file = 'ACH',cell_line_ann_version = NULL, cell_line_ann_taiga=TRUE,
                          tumor_data_name = 'celligner-input-9827', tumor_data_file = 'tumor_expression', tumor_version = NULL, tumor_taiga = TRUE,
                          tumor_ann_name = 'celligner-input-9827', tumor_ann_file = 'tumor_annotations', tumor_ann_version = NULL, tumor_ann_taiga = TRUE,
                          additional_annotations_name = 'celligner-input-9827', additional_annotations_file = 'CCLE_annotations', additional_annotations_version = NULL, additional_annotations_taiga = TRUE,
                          hgnc_data_name = 'hgnc-87ab', hgnc_data_file='hgnc_complete_set', hgnc_version= NULL, hgnc_taiga = TRUE,
                          met500_data_name = 'met500-fc3c', met500_data_file = 'met500_TPM', met500_version = NULL, met500_taiga = TRUE,
                          met500_ann_name = 'met500-fc3c', met500_ann_file = 'met500_ann', met500_ann_version = NULL, met500_ann_taiga = TRUE,
                          Novartis_PDX_data_name = 'pdx-data-3d29', Novartis_PDX_data_file = 'Novartis_PDX_TPM', Novartis_PDX_version = NULL, Novartis_PDX_taiga = TRUE,
                          Novartis_PDX_ann_name = 'pdx-data-3d29', Novartis_PDX_ann_file = 'Novartis_PDX_ann', Novartis_PDX_ann_version = NULL, Novartis_PDX_ann_taiga = TRUE,
                          pediatric_PDX_data_name = 'pdx-data-3d29', pediatric_PDX_data_file = 'pediatric_PDX_TPM', pediatric_PDX_version = NULL, pediatric_PDX_taiga = TRUE,
                          pediatric_PDX_ann_name = 'pdx-data-3d29', pediatric_PDX_ann_file = 'pediatric_PDX_ann', pediatric_PDX_ann_version = NULL, pediatric_PDX_ann_taiga = TRUE,
                          save_output = NULL) {

  require(magrittr)
  require(tidyverse)

  dat <- load_data(cell_line_data_name, cell_line_data_file, cell_line_version, cell_line_taiga,
                   cell_line_ann_name, cell_line_ann_file,cell_line_ann_version, cell_line_ann_taiga,
                   tumor_data_name, tumor_data_file, tumor_version, tumor_taiga,
                   tumor_ann_name, tumor_ann_file, tumor_ann_version, tumor_ann_taiga,
                   additional_annotations_name, additional_annotations_file, additional_annotations_version, additional_annotations_taiga,
                   hgnc_data_name, hgnc_data_file, hgnc_version, hgnc_taiga)

  met500 <- load_additional_data(met500_data_name, met500_data_file, met500_data_version, met500_data_taiga,
                                 met500_ann_name, met500_ann_file, met500_ann_version, met500_ann_taiga)
  Novartis_PDX <- load_additional_data(Novartis_PDX_data_name, Novartis_PDX_data_file, Novartis_PDX_data_version, Novartis_PDX_data_taiga,
                                       Novartis_PDX_ann_name, Novartis_PDX_ann_file, Novartis_PDX_ann_version, Novartis_PDX_ann_taiga)

  pediatric_PDX <- load_additional_data(pediatric_PDX_data_name, pediatric_PDX_data_file, pediatric_PDX_data_version, pediatric_PDX_data_taiga,
                                       pediatric_PDX_ann_name, pediatric_PDX_ann_file, pediatric_PDX_ann_version, pediatric_PDX_ann_taiga)

  shared_genes <- intersect(colnames(dat$TCGA_mat), colnames(dat$CCLE_mat)) %>%
    intersect(colnames(met500$mat)) %>%
    intersect(colnames(Novartis_PDX$mat)) %>%
    intersect(colnames(pediatric_PDX$mat))

  dat$TCGA_mat <- dat$TCGA_mat[,shared_genes]
  dat$CCLE_mat <- dat$CCLE_mat[,shared_genes]
  met500$mat <- met500$mat[,shared_genes]
  Novartis_PDX$mat <- Novartis_PDX$mat[,shared_genes]
  pediatric_PDX$mat <- pediatric_PDX$mat[,shared_genes]

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
    dplyr::filter(best_rank < celligner_global$top_DE_genes_per) %>%
    .[['Gene']]


  cov_diff_eig <- run_cPCA(TCGA_obj, CCLE_obj, celligner_global$fast_cPCA)

  if(is.null(celligner_global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, celligner_global$remove_cPCA_dims, drop = FALSE]
  } else {
    cur_vecs <- cov_diff_eig$rotation[, celligner_global$remove_cPCA_dims, drop = FALSE]
  }

  rownames(cur_vecs) <- colnames(dat$TCGA_mat)
  TCGA_cor <- resid(lm(t(dat$TCGA_mat) ~ 0 + cur_vecs)) %>% t()
  CCLE_cor <- resid(lm(t(dat$CCLE_mat) ~ 0 + cur_vecs)) %>% t()

  mnn_res <- run_MNN(CCLE_cor, TCGA_cor,  k1 = celligner_global$mnn_k_tumor, k2 = celligner_global$mnn_k_CL, ndist = celligner_global$mnn_ndist,
                     subset_genes = DE_gene_set)

  combined_mat <- rbind(mnn_res$corrected, CCLE_cor)

  # clear unused objects
  rm(TCGA_obj); rm(CCLE_obj); rm(cov_diff_eig); rm(TCGA_cor); rm(CCLE_cor); gc()

  # Met500 alignment
  met500_cor <- resid(lm(t(met500$mat) ~ 0 + cur_vecs)) %>% t()

  mnn_res <- run_MNN(combined_mat, met500_cor,  k1 = 20, k2 = 50, ndist = celligner_global$mnn_ndist,
                     subset_genes = DE_gene_set)
  combined_mat <- rbind(combined_mat, mnn_res$corrected)

  ## align PDX datasets

  ### PDX - Novartis
  Novartis_PDX_cor <- resid(lm(t(Novartis_PDX$mat) ~ 0 + cur_vecs)) %>% t()

  mnn_res_Novartis_PDX <- run_MNN(combined_mat, Novartis_PDX_cor, k1 = 10, k2 = 50, ndist = 3,
                         subset_genes = DE_gene_set)

  combined_mat <- rbind(combined_mat, mnn_res_Novartis_PDX$corrected)

  ### PDX - pediatric
  pediatric_PDX_cor <- resid(lm(t(pediatric_PDX$mat) ~ 0 + cur_vecs)) %>% t()

  mnn_res_pediatric_PDX <- run_MNN(combined_mat[-which(rownames(combined_mat) %in% rownames(pediatric_PDX_cor)),],
                               pediatric_PDX_cor, k1 = 10, k2 = 50, ndist = 3,
                               subset_genes = DE_gene_set)

  combined_mat <- t(rbind(combined_mat, mnn_res_pediatric_PDX$corrected))

  # combine all output
  comb_ann <- rbind.data.frame(comb_ann[,c('sampleID', 'lineage', 'subtype', 'type')],
                               met500$ann[,c('sampleID', 'lineage', 'subtype', 'type')],
                               Novartis_PDX$ann[,c('sampleID', 'lineage', 'subtype', 'type')],
                               pediatric_PDX$ann[,c('sampleID', 'lineage', 'subtype', 'type')])
  rownames(comb_ann) <- comb_ann$sampleID
  comb_ann <- comb_ann[colnames(combined_mat),]

  # clear unused object
  rm(met500); rm(Novartis_PDX); rm(pediatric_PDX); rm(met500_cor); rm(Novartis_PDX_cor); rm(pediatric_PDX_cor); gc()

   # create seurat object
  comb_obj <- create_Seurat_object(combined_mat, comb_ann)
  comb_obj <- cluster_data(comb_obj)

  Celligner_res <- Seurat::Embeddings(comb_obj, reduction = 'umap') %>%
    as.data.frame() %>%
    magrittr::set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    tibble::rownames_to_column(var = 'sampleID') %>%
    dplyr::left_join(comb_obj@meta.data, by = 'sampleID')

  lineage_averages <- Celligner_res %>%
    dplyr::filter(!lineage %in% c('adrenal_cortex', 'embryo', 'endocrine', 'engineered', 'engineered_blood',
                                  'engineered_breast', 'engineered_central_nervous_system', 'engineered_kidney',
                                  'engineered_lung', 'engineered_ovary', 'engineered_prostate', 'epidermoid_carcinoma',
                                  'nasopharynx', 'nerve','pineal', 'teratoma', 'unknown')) %>%
    dplyr::group_by(lineage) %>%
    dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm=T),
                     UMAP_2 = median(UMAP_2, na.rm=T))
  lineage_averages$lineage <- gsub("_", " ", lineage_averages$lineage)
  lineage_lab_aes <- ggplot2::geom_text(data = lineage_averages, mapping = aes(x = UMAP_1, y = UMAP_2, label = lineage), size = 3, color="#000000")


  if('type' %in% colnames(Celligner_res) & 'tumor' %in% Celligner_res$type & 'CL' %in% Celligner_res$type) {
    celligner_plot <- ggplot2::ggplot(Celligner_res,  ggplot2::aes(UMAP_1, UMAP_2)) +
      ggplot2::geom_point(alpha=0.7, pch=21,  ggplot2::aes(color = type, fill = lineage, size = type)) +
      ggplot2::scale_color_manual(values = c(tumor = 'white', CL = 'black')) +
      ggplot2::scale_size_manual(values=c(tumor=0.75, CL=1.5)) +
      ggplot2::xlab('UMAP 1') + ggplot2::ylab('UMAP 2') +
      ggplot2::guides(fill=FALSE,
                      color = ggplot2::guide_legend(override.aes = list(color=c('black', 'white'), fill = c('white','black')))) +
      ggplot2::theme_classic()
  } else {
    celligner_plot <-  ggplot2::ggplot(Celligner_res, ggplot2::aes(UMAP_1, UMAP_2)) +
      ggplot2::geom_point(alpha=0.7, pch=21, size = 1, ggplot2::aes(fill = lineage)) +
      ggplot2::xlab('UMAP 1') + ggplot2::ylab('UMAP 2') +
      ggplot2::theme_classic() + ggplot2::theme(legend.position = 'none')
  }

  print(celligner_plot)
  print(celligner_plot + lineage_lab_aes)


  if(!is.null(save_output)) {
    if(file.exists(save_output)) {
      print('saving files')
      write.csv(combined_mat, file.path(save_output, 'Celligner_multidataset_aligned_data.csv'))
      readr::write_csv(Celligner_res, file.path(save_output, 'Celligner_multidataset_info.csv'))
      ggplot2::ggsave(file.path(save_output, 'Celligner_multidataset_plot.png'), celligner_plot, device='png', width = 8, height = 6)
      ggplot2::ggsave(file.path(save_output, 'labeled_Celligner_multidataset_plot.png'), celligner_plot + lineage_lab_aes, device='png', width = 8, height = 6)

    } else{
      warning("can't save output, folder does not exist")
    }
  }

  return(comb_obj)
}

