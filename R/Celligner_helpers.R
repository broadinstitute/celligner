library(magrittr)
library(tidyverse)

#' check for NAs in the expression data and remove samples with NAs
#' @name check_NAs
#'
#' @param mat: matrix of gene expression data that is samples by genes
#' @return matrix of gene expression data, removing samples that have NAs
#' @export
#'
check_NAs <- function(mat) {
  if(length(which(is.na(rowSums(mat))==T))>0) {
    warning("Removing sample(s) due to NAs in the data")
    mat <- mat[!is.na(rowSums(mat)),]
  }

  return(mat)
}

#'
#' Differentially expressed genes
#' @name run_lm_stats_limma_group
#'
#' @param mat: Nxp data matrix of N cell lines and p genes
#' @param phenos: N vector of independent variables. Can be two-group labels as factors, bools, or can be numeric
#' @param covars: optional Nxk matrix of sample covariates
#' @param weights: optional N vector of precision weights for each data point
#' @param target_type: name of the column variable in the data (default 'Gene')
#' @return table of gene level stata
#' @description  Estimate linear-model stats for a matrix of data with respect to a group of phenotype variables
# using limma with empirical Bayes moderated F-stats for p-values
#' @export
#'
run_lm_stats_limma_group <- function (mat, phenos, covars = NULL, weights = NULL, target_type = "Gene",
          limma_trend = FALSE)
{
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  udata <- rownames(mat) %>% intersect(rownames(phenos))
  if (!is.null(covars)) {
    udata %<>% intersect(rownames(covars))
  }
  form <- as.formula(paste("~", paste0(colnames(phenos), collapse = " + ")))
  design <- model.matrix(form, data = phenos[udata, , drop = F])
  if (!is.null(covars)) {
    covars <- data.frame(covars)
    form <- as.formula(paste("~", paste0(colnames(covars),
                                         collapse = " + ")))
    Cdesign <- model.matrix(form, data = covars[udata, ,
                                                drop = F])
    Cdesign <- Cdesign[, setdiff(colnames(Cdesign), "(Intercept)"),
                       drop = FALSE]
    stopifnot(length(intersect(colnames(Cdesign), colnames(design))) ==
                0)
    design %<>% cbind(Cdesign)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata, ])
    }
    else {
      weights <- weights[udata]
    }
  }
  design <- design[, colSums(design) > 2, drop = FALSE]
  targ_coefs <- setdiff(colnames(design), "(Intercept)")
  fit <- limma::lmFit(t(mat[udata, ]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- which(colnames(design) %in% targ_coefs)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf,
                             sort.by = "F", genelist = colnames(mat))
  results %<>% tibble::rownames_to_column(var = target_type)
  results %<>% magrittr::set_colnames(revalue(colnames(.), c(AveExpr = "Avg",
                                                   F = "F_stat", P.Value = "p.value", adj.P.Val = "q.value"))) %>%
    na.omit() %>% dplyr::select(-ProbeID)
  return(results)
}

#'
#' cPCA
#' @name run_cPCA_analysis
#'
#' @param TCGA_dat: sample by genes matrix of scaled expression data
#' @param CCLE_dat: sample by genes matrix of scaled expression data
#' @param tumor_cluster_df: table of sample metadata that includes a column 'seurat_clusters',
#'  containing transcriptional clusters in the TCGA data
#' @param CL_cluster_df: table of sample metadata that includes a column 'seurat_clusters',
#'  containing transcriptional clusters in the CCLE data
#' @param pc_dims: numbers of cPCs calculated. If set to NULL (default) all cPCs will be calculated, if set to a value
#' then that number of cPCs will be approximated. Values input should be >= 4.
#' @return contrastive principal component object containing cPC vectors and values
#' @description Run contrastive principal components analysis, first removing average cluster expression, to
# estimate the average intra-cluster covariance. If pc_dims = NULL, all cPCs are calculated. Faster cPCA can be run by setting pc_dims to a
# value >=4 and approximating just those cPCs.
#' @export
#'
run_cPCA_analysis <- function(TCGA_dat, CCLE_dat, tumor_cluster_df, CL_cluster_df, pc_dims=NULL) {
  tumor_clust_avgs <- get_cluster_averages(TCGA_dat, tumor_cluster_df)
  CL_clust_avgs <- get_cluster_averages(CCLE_dat, CL_cluster_df)

  TCGA_subtype_ms <- TCGA_dat - tumor_clust_avgs[tumor_cluster_df$seurat_clusters,]
  CCLE_subtype_ms <- CCLE_dat - CL_clust_avgs[CL_cluster_df$seurat_clusters,]

  TCGA_cov <- cov(TCGA_subtype_ms)
  CCLE_cov <- cov(CCLE_subtype_ms)

  if(!is.null(pc_dims)) {
    cov_diff_eig <- irlba::prcomp_irlba(TCGA_cov - CCLE_cov, n = pc_dims)
  } else {
    cov_diff_eig <- eigen(TCGA_cov - CCLE_cov)
  }
  return(cov_diff_eig)
}

#'
#' calculate the average expression per cluster
#' @name get_cluster_averages
#'
#' @param mat: sample by genes matrix of expression data
#' @param cluster_df: table of sample metadata that includes a column 'seurat_clusters',
#' containing transcriptional clusters
#' @return average cluster expression
#' @description calculate the average expression per cluster
#' @export
#'
get_cluster_averages <- function(mat, cluster_df) {
  n_clusts <- nlevels(cluster_df$seurat_clusters)
  clust_avgs <- matrix(NA, nrow = n_clusts, ncol = ncol(mat)) %>%
    magrittr::set_colnames(colnames(mat)) %>%
    magrittr::set_rownames(levels(cluster_df$seurat_clusters))
  for (ii in levels(cluster_df$seurat_clusters)) {
    clust_avgs[ii,] <- colMeans(mat[cluster_df$seurat_clusters == ii,], na.rm=T)
  }
  return(clust_avgs)
}

# MNN --------------------------------------------------------------------

#'
#' MNN
#' @name modified_mnnCorrect
#'
#' @param ref_mat: matrix of samples by genes of cPC corrected data that serves as the reference data in the MNN alignment.
#' In the standard Celligner pipeline this the cell line data.
#' @param targ_mat: matrix of samples by genes of cPC corrected data that is corrected in the MNN alignment and projected onto the reference data.
#' In the standard Celligner pipeline this the tumor data.
#' @param k1: the number of neighbors within the data being corrected (in standard pipeline the tumor data). By default this is 20.
#' @param k2: the number of neighbors within the reference data (in standard pipeline the cell line data). By default this is 20.
#' @param ndist: A numeric scalar specifying the threshold beyond which neighbors are to be ignored when computing correction vectors.
#' By default is 3.
#' @param subset_genes: the subset of genes used for identifying mutual nearest neighbors within the datasets. The set of differentially
#' expressed genes is usually passed here. By default is NULL, meaning all genes are used
#' @return MNN object, containing the targ_mat corrected data and the mutual nearest neighbor pairs.
#' @description Mutual nearest neighbors correction. Modification of the scran::fastMNN (https://github.com/MarioniLab/scran).
#' Allows for separate k values per dataset, and simplifies some of the IO and doesn't use PCA reduction
#' @export
#'
modified_mnnCorrect <- function(ref_mat, targ_mat, k1 = 20, k2 = 20,
                            ndist = 3, subset_genes = NULL) {
  if (is.null(subset_genes)) {
    subset_genes <- colnames(ref_mat)
  }

  sets <- batchelor::findMutualNN(ref_mat[, subset_genes],
                                  targ_mat[, subset_genes],
                                  k1 = k2, k2 = k1,
                                  BPPARAM = BiocParallel::SerialParam())
  mnn_pairs <- as.data.frame(sets) %>%
    dplyr::mutate(ref_ID = rownames(ref_mat)[first],
           targ_ID = rownames(targ_mat)[second],
           pair = seq(nrow(.))) %>%
    dplyr::select(-first, -second)

  # Estimate the overall batch vector.
  ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  overall.batch <- colMeans(ave.out$averaged)

  #remove variation along the overall batch vector
  ref_mat <- .center_along_batch_vector(ref_mat, overall.batch)
  targ_mat <- .center_along_batch_vector(targ_mat, overall.batch)

  # Recompute correction vectors and apply them.
  re.ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  targ_mat <- .tricube_weighted_correction(targ_mat, re.ave.out$averaged, re.ave.out$second, k=k2, ndist=ndist, subset_genes, BPPARAM=BiocParallel::SerialParam())

  final <- list(corrected = targ_mat,
                pairs = mnn_pairs)
  return(final)
}

#'
#' calculate the average correction vector
#' @name .average_correction
#'
#' @param refdata: matrix of samples by genes of cPC corrected data that serves as the reference data in the MNN alignment.
#' In the standard Celligner pipeline this the cell line data.
#' @param mnn1: mnn1 pairs
#' @param curdata: matrix of samples by genes of cPC corrected data that is corrected in the MNN alignment and projected onto the reference data.
#' In the standard Celligner pipeline this the tumor data.
#' @param mnn2: mnn2 pairs
#' @return correction vector and pairs
#' @description Computes correction vectors for each MNN pair, and then averages them for each MNN-involved cell in the second batch.
#' Copied from dev version of scran (2018-10-28), with slight modifications as noted https://github.com/MarioniLab/scran
#' @export
#'
.average_correction <- function(refdata, mnn1, curdata, mnn2)
  # Computes correction vectors for each MNN pair, and then
  # averages them for each MNN-involved cell in the second batch.
{
  corvec <- refdata[mnn1,,drop=FALSE] - curdata[mnn2,,drop=FALSE]
  corvec <- rowsum(corvec, mnn2)
  npairs <- table(mnn2)
  stopifnot(identical(names(npairs), rownames(corvec)))
  corvec <- unname(corvec)/as.vector(npairs)
  list(averaged=corvec, second=as.integer(names(npairs)))
}


#'
#' centers samples within each batch
#' @name .center_along_batch_vector
#'
#' @param mat: matrix of samples by genes
#' @param batch.vec: batch vector
#' @return correction vector and pairs
#' @description Projecting along the batch vector, and shifting all samples to the center within each batch.
#' This removes any variation along the overall batch vector within each matrix.
#' @export
#'
.center_along_batch_vector <- function(mat, batch.vec)
  # Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
  # This removes any variation along the overall batch vector within each matrix.
{
  batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
  batch.loc <- as.vector(mat %*% batch.vec)
  central.loc <- mean(batch.loc)
  mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
  return(mat)
}


#' tricube-weighted correction
#' @name .tricube_weighted_correction
#'
#' @param curdata: target matrix of samples by genes
#' @param correction: corrected vector
#' @param in.mnn: mnn pairs
#' @param k: k values, default 20
#' @param ndist: A numeric scalar specifying the threshold beyond which neighbors are to be ignored when computing correction vectors.
#' By default is 3.
#' @param subset_genes: genes used to identify mutual nearest neighbors
#' @param BNPARAM: default NULL
#' @param BPPARAM: default BiocParallel::SerialParam()
#' @return MNN corrected data
#' @description Computing tricube-weighted correction vectors for individual samples,
#' using the nearest neighbouring samples involved in MNN pairs.
#' Modified to use FNN rather than queryKNN for nearest neighbor finding
#' @export
#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocParallel SerialParam
#'
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, subset_genes, BNPARAM=NULL, BPPARAM=BiocParallel::SerialParam())
  # Computing tricube-weighted correction vectors for individual cells,
  # using the nearest neighbouring cells _involved in MNN pairs_.
  # Modified to use FNN rather than queryKNN for nearest neighbor finding
{
  cur.uniq <- curdata[in.mnn,,drop=FALSE]
  safe.k <- min(k, nrow(cur.uniq))
  # closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
  closest <- FNN::get.knnx(cur.uniq[, subset_genes], query=curdata[, subset_genes], k=safe.k)
  # weighted.correction <- .compute_tricube_average(correction, closest$index, closest$distance, ndist=ndist)
  weighted.correction <- .compute_tricube_average(correction, closest$nn.index, closest$nn.dist, ndist=ndist)
  curdata + weighted.correction
}

#'
#' compute tricube averages
#' @name .compute_tricube_average
#'
#' @param values: correction vector
#' @param indices: nxk matrix for the nearest neighbor indice
#' @param distances: nxk matrix for the nearest neighbor Euclidea distances
#' @param bandwidth: Is set at 'ndist' times the median distance, if not specified.
#' @param ndist: By default is 3.
#' @description Centralized function to compute tricube averages.
#' @export
#'
.compute_tricube_average <- function(vals, indices, distances, bandwidth=NULL, ndist=3)
  # Centralized function to compute tricube averages.
  # Bandwidth is set at 'ndist' times the median distance, if not specified.
{
  if (is.null(bandwidth)) {
    middle <- ceiling(ncol(indices)/2L)
    mid.dist <- distances[,middle]
    bandwidth <- mid.dist * ndist
  }
  bandwidth <- pmax(1e-8, bandwidth)

  rel.dist <- distances/bandwidth
  rel.dist[rel.dist > 1] <- 1 # don't use pmin(), as this destroys dimensions.
  tricube <- (1 - rel.dist^3)^3
  weight <- tricube/rowSums(tricube)

  output <- 0
  for (kdx in seq_len(ncol(indices))) {
    output <- output + vals[indices[,kdx],,drop=FALSE] * weight[,kdx]
  }

  if (is.null(dim(output))) {
    matrix(0, nrow(vals), ncol(vals))
  } else {
    output
  }
}
