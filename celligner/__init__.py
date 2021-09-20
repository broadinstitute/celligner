# Celligner
from celligner.params import *
from genepy.utils import helper as h
from genepy.utils import plot
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.linear_model import LinearRegression
import mnnpy
from scanpy.tl import louvain
from scanpy.pp import neighbors
from anndata import AnnData
# import louvain
# import pynndescent
import os
import sys

import pandas as pd
import numpy as np
import umap.umap_ as umap
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/snn/")
from SNN import snn
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import limma
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/contrastive/")
from contrastive import CPCA

def check_Xpression(X_pression, gene_file):
  """

  Args:
      X_pression (pd.Dataframe): [description]
      gene_file (pd.Dataframe): [description]

  Raises:
      ValueError: if the number of genes in the expression matrix and the gene file do not match
      ValueError: if the expression matrix contains nan values

  Returns:
      pd.Dataframe: 
  """
  common_genes = set(X_pression.columns) & set(gene_file.ensembl_gene_id)
  if len(common_genes) < MIN_GENES:
    raise ValueError("X_pression columns do not match gene_file enough only " +
                    str(len(common_genes))+" common genes")
  else:
    print("found "+str(len(common_genes))+" common genes")
  # drop genes not in gene_file
  X_pression = X_pression.loc[:, common_genes].astype(float)
  # raise issue if there are any NaNs
  if X_pression.isnull().values.any():
    raise ValueError("X_pression contains NaNs")
  return X_pression


def runDiffExprOnCluster(expression, clustered, clust_covariates=None, pvalue_threshold=DESEQ_MAX_PVAL,):
  """
  Runs DESEQ2 on the clustered data.

  Args:
    expression: pandas.DataFrameThe expression data.
    clustered: list 

  """
  n_clusts = len(set(clustered))
  print('running differential expression on '+str(n_clusts)+' clusters')
  clusts = set(clustered)-set([-1])
  #TODO: add covariates
  if clust_covariates:
    if len(clust_covariates) != n_clusts:
      raise ValueError("number of covariates does not match number of clusters")
    design_matrix = clust_covariates
  # make a design matrix
  design_matrix = pd.DataFrame(index=expression.index,
                              data=np.array([clustered==i for i in clusts]).T,
                               columns=["C"+ str(i)+"C" for i in clusts])
  design_matrix.index = design_matrix.index.astype(str).str.replace('-', '.')
  design_matrix = design_matrix[design_matrix.sum(1) > 0]
  # creating the matrix
  data = expression.T
  data = data[data.columns[clustered!=-1].tolist()]
  # running limmapy
  print("running limmapy on the samples")
  res = limma.limmapy().lmFit(data, design_matrix).eBayes(trend=False).topTable(
      number=len(data)).iloc[:,len(clusts):]
  return res.sort_values(by='F')

def marioniCorrect(ref_mat, targ_mat, k, ndist, subset_genes, **mnn_kwargs):
  """marioniCorrect is a function that corrects for batch effects using the Marioni method.

  Args:
    ref_mat (pd.Dataframe): matrix of samples by genes of cPC corrected data that serves as the reference data in the MNN alignment.
      In the standard Celligner pipeline this the cell line data.
    targ_mat matrix of samples by genes of cPC corrected data that is corrected in the MNN alignment and projected onto the reference data.
      In the standard Celligner pipeline this the tumor data.
    mnn_kwargs (dict): args to mnnCorrect

  Returns:
      pd.Dataframe: corrected dataframe
  """
  if "var_subset" in mnn_kwargs:
    ref_mat = ref_mat.loc[:, mnn_kwargs["var_subset"]]
    targ_mat = targ_mat.loc[:, mnn_kwargs["var_subset"]]
  # find mutual nearest neighbors
  corrected, mnn_pairs, _  = mnnpy.mnn_correct(ref_mat, targ_mat, **mnn_kwargs)
  # make a dataframe
  mnn_pairs = pd.DataFrame(mnn_pairs).rename(columns={'first':'ref_ID', 'second':'targ_ID', 'pair':'pair'})
  # compute the overall batch vector
  ave_out = _averageCorrection(ref_mat, sets['first'], targ_mat, sets['second'])
  overall_batch = ave_out['averaged'].mean(axis=0)
  # remove variation along the overall batch vector
  ref_mat = _centerAlongBatchVector(ref_mat, overall_batch)
  targ_mat = _centerAlongBatchVector(targ_mat, overall_batch)
  # recompute correction vectors and apply them
  re_ave_out = _averageCorrection(ref_mat, sets['first'], targ_mat, sets['second'])
  targ_mat = _tricubeWeightedCorrection(targ_mat, re_ave_out['averaged'], re_ave_out['second'], k=k, ndist=ndist)
  final = {'corrected':targ_mat, 'pairs':mnn_pairs}
  return final

def _averageCorrection(refdata, mnn1, curdata, mnn2):
  """_averageCorrection computes correction vectors for each MNN pair, and then averages them for each MNN-involved cell in the second batch.

  Args:
      refdata (pandas.DataFrame): matrix of samples by genes of cPC corrected data that serves as the reference data in the MNN alignment.
      mnn1 (list): mnn1 pairs
      curdata (pandas.DataFrame): matrix of samples by genes of cPC corrected data that is corrected in the MNN alignment and projected onto the reference data.
      mnn2 (list): mnn2 pairs

  Returns:
      dict: correction vector and pairs
  """
  corvec = refdata.iloc[mnn1] - curdata.iloc[mnn2]
  corvec = corvec.sum(axis=0)
  npairs = pd.Series(mnn2).value_counts()
  #stopifnot(npairs.index.equals(corvec.index))
  corvec = corvec/npairs
  return {'averaged':corvec, 'second':npairs.index.astype(int)}

def _centerAlongBatchVector(mat, batch_vec):
  """_centerAlongBatchVector - Projecting along the batch vector, and shifting all samples to the center within each batch.

  Args:
      mat (pandas.DataFrame): matrix of samples by genes
      batch_vec (pandas.Series): batch vector

  Returns:
      pandas.DataFrame: corrected matrix
  """
  batch_vec = batch_vec/np.sqrt(np.sum(batch_vec**2))
  batch_loc = np.dot(mat, batch_vec)
  central_loc = np.mean(batch_loc)
  mat = mat + np.outer(central_loc - batch_loc, batch_vec)
  return mat

def _tricubeWeightedCorrection(curdata, correction, in_mnn, k=20, ndist=3):
  """_tricubeWeightedCorrection computes tricube-weighted correction vectors for individual cells,

  Args:
      curdata (pandas.DataFrame): target matrix of samples by genes
      correction (pandas.Series): corrected vector
      in_mnn (pandas.DataFrame): mnn pairs
      k (int, optional): k values, default 20
      ndist (int, optional): A numeric scalar specifying the threshold beyond which neighbors are to be ignored when computing correction vectors.
      subset_genes (list, optional): genes used to identify mutual nearest neighbors
      BNPARAM (None, optional): default None
      BPPARAM (None, optional): default BiocParallel::SerialParam()
  """
  cur_uniq = curdata.iloc[in_mnn]
  safe_k = min(k, cur_uniq.shape[0])
  # closest = queryKNN(query=curdata, X=cur_uniq, k=safe_k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
  closest = []#FNN.get_knnx(cur_uniq[:, subset_genes], query=curdata[:, subset_genes], k=safe_k)
  # weighted_correction = compute_tricube_average(correction, closest['index'], closest['distance'], ndist=ndist)
  weighted_correction = _compute_tricube_average(correction, closest['nn_index'], closest['nn_dist'], ndist=ndist)
  curdata + weighted_correction


def _compute_tricube_average(vals, indices, distances, bandwidth=None, ndist=3):
  """_compute_tricube_average - Centralized function to compute tricube averages.

  Args:
      vals (pandas.DataFrame): correction vector
      indices (pandas.DataFrame): nxk matrix for the nearest neighbor indice
      distances (pandas.DataFrame): nxk matrix for the nearest neighbor Euclidea distances
      bandwidth (float): Is set at 'ndist' times the median distance, if not specified.
      ndist (int, optional): By default is 3.

  Returns:
      [type]: [description]
  """
  if bandwidth is None:
    middle = int(np.ceil(indices.shape[1]/2))
    mid_dist = distances[:,middle]
    bandwidth = mid_dist * ndist
  bandwidth = np.maximum(1e-8, bandwidth)

  rel_dist = distances/bandwidth
  rel_dist[rel_dist > 1] = 1 # don't use pmin(), as this destroys dimensions.
  tricube = (1 - rel_dist**3)**3
  weight = tricube/np.sum(tricube, axis=0)

  output = 0
  for kdx in range(indices.shape[1]):
    output = output + vals[indices[:,kdx]] * weight[:,kdx]

  if output.shape[0] == 0:
    output = np.zeros((vals.shape[0], vals.shape[1]))
  return output


class Celligner(object):
  def __init__(self, args={}, gene_file=None, onlyGenes=GENE_TYPE,
               ensemble_server="http://nov2020.archive.ensembl.org/biomart",
               umap_kwargs=UMAP_PARAMS, pca_kwargs=PCA_PARAMS,
               snn_kwargs=SNN_PARAMS, topKGenes=TOP_K_GENES, cpca_kwargs=CPCA_PARAMS, 
               cpca_ncomp=CPCA_NCOMP, mnn_kwargs=MNN_PARAMS, make_plots=False,
               low_mem=False,): #scneigh_kwargs=SCneigh_PARAMS
    """initialize Celligner object

    Args:
        args ([type]): [description]
        onlyGenes (str, optional): one of 'usefull', 'all', 'protein_coding'. Defaults to "usefull".
        gene_file (pd.Dataframe, optional): Needs to contain at least 15000 genes 
          and an "ensembl_gene_id", columns. Defaults to None.
        ensemble_server (str, optional): [description]. Defaults to "http://nov2020.archive.ensembl.org/biomart".
        umap_kwargs (dict, optional): see umap_pamarameters.md or . Defaults to {}.
    """
    self.args = args
    if gene_file:
      self.gene_file = gene_file
    else:
      self.gene_file = h.generateGeneNames(ensemble_server=ensemble_server,
      useCache=True)
      if onlyGenes == "protein_coding":
        print('using only protein coding genes')
        self.gene_file[self.gene_file.gene_biotype == "protein_coding"]
      elif onlyGenes == "usefull":
        print("using only usefull genes")
        self.gene_file[self.gene_file.gene_biotype.isin(USEFUL_GENE_BIOTYPES)]
      else:
        print('using all genes')
      self.gene_file.ensembl_gene_id.drop_duplicates(
          keep='first', inplace=True)
    self.umap_kwargs = umap_kwargs
    self.pca_kwargs = pca_kwargs
    self.snn_kwargs = snn_kwargs
    self.topKGenes = topKGenes
    self.cpca_kwargs = cpca_kwargs
    self.cpca_ncomp = cpca_ncomp
    self.mnn_kwargs = mnn_kwargs
    self.number_of_datasets = 0
    self.make_plots = make_plots
    self.low_mem = low_mem
    #self.scneigh_kwargs = scneigh_kwargs

    self.fit_input=None
    self.fit_reduced=None
    self.fit_clusters=None
    self.differential_genes_input=None
    self.differential_genes_names=None
    self.fit_annotations = None
    self.transform_annotations = None
    self.transform_input=None
    self.transform_clusters=None
    self.transform_reduced=None
    self.corrected=None
    self.common_genes = None


  def addToFit(self, X_pression, annotations=None, dofit=True, doAdd=True):
    """adds expression data to the fit dataframe

    Args:
        X_pression (pd.Dataframe): 
        annotations (pd.Dataframe, optional): [description]. Defaults to None.

    Raises:
        ValueError: if the expression matrix and annotations matrix do not have the same index
        ValueError: if the new expression matrix has different gene names than the current one
    """
    count = X_pression.shape[0]+(self.fit_input.shape[0] if self.fit_input is not None else 0)
    print('looking at '+str(count)+' samples.')
    fit_input = check_Xpression(X_pression, self.gene_file)
    if annotations is not None:
      if len(annotations) != len(fit_input) or list(fit_input.index) != list(annotations.index):
        raise ValueError("annotations do not match X_pression")
    else:
      # create fake annotations
      annotations = pd.DataFrame(index=X_pression.index,
                                 columns=['cell_type',
                                          'disease_type', 'tissue_type'],
                                 data=np.zeros((len(X_pression), 3))+self.number_of_datasets)
      
    if self.common_genes is None or not doAdd:
      # it is the first time we run it.
      print("creating a fit dataset..")
      self.common_genes = fit_input.columns
      self.fit_input = fit_input
      self.fit_annotations = annotations
    else:
      print("adding to fit dataset..")
      # check if same genes
      if set(fit_input.columns) != set(self.common_genes):
        raise ValueError("fit_input's genes do not match common_genes")
      # add annotations together
      self.fit_annotations = self.fit_annotations.append(annotations)
      # add fit_input together
      self.fit_input = self.fit_input.append(fit_input)
    self.number_of_datasets +=1
    if dofit:
      return self.fit()


  def fit(self, X_pression=None, annotations=None):
    """fit the model using X_pression

    Args:
        X_pression (pd.Dataframe): contains the expression data as RSEM expected counts with 
          ensembl_gene_id as columns and samplenames as index.
        annotations (pd.Dataframe, optional): sample annotations, for each sample, 
          needs to contain ['cell_type', 'disease_type', 'tissue_type']. 
          Defaults to None (will create an empty dataframe).
    """
    # check if X_pression is compatible with the model
    if X_pression is not None:
      self.addToFit(X_pression, annotations, dofit=False, doAdd=False)
    elif self.fit_input is None:
      raise ValueError("no input provided")
  
    # mean center the dataframe
    fit_input = self.fit_input.sub(self.fit_input.mean(axis=1), axis=0)
    # dimensionality reduction
    print('reducing dimensionality...')
    pca = PCA(**self.pca_kwargs) if not self.low_mem else IncrementalPCA(**self.pca_kwargs)
    self.fit_reduced = pca.fit_transform(fit_input)
    # clustering: doing SNN on the reduced data
    print('clustering...')
    #anndata from df
    adata = AnnData(self.fit_reduced)
    neighbors(adata) # **scneigh_kwargs
    louvain(adata)
    self.fit_clusters = adata.obs['louvain'].values
    del adata
    #self.fit_clusters = snn.SNN(**self.snn_kwargs).fit_predict(self.fit_reduced,
    #    sample_weight=None)
    # do differential expression between clusters and getting the top K most expressed genes
    if self.make_plots:
      # plotting
      plot.scatter(umap.UMAP(
        **self.umap_kwargs).fit_transform(self.fit_reduced), 
        xname="UMAP1", yname="UMAP2", colors=self.fit_clusters,
        labels=["C"+str(i) for i in self.fit_clusters],
        title="SNN clusters", radi=.1)
    
    print('doing differential expression analysis on the clusters')
    if len(set(self.fit_clusters)) < 2:
      raise ValueError("only one cluster found, no differential expression possible\
        try to change your parameters...")
    self.differential_genes_input = runDiffExprOnCluster(
        self.fit_input, self.fit_clusters)
    # need enough genes to be significant
    if len(self.differential_genes_input[self.differential_genes_input.F>10]) < self.topKGenes:
      raise ValueError("not enough differentially expressed genes found..")
    print('done')
    return self


  def addTotransform(self, X_pression, annotations=None, dotransform=True, doAdd=True):
    """adds expression data to the transform dataframe

    Args:
        X_pression (pd.Dataframe): 
        annotations (pd.Dataframe, optional): [description]. Defaults to None.

    Raises:
        ValueError: if the expression matrix and annotations matrix do not have the same index
        ValueError: if the new expression matrix has different gene names than the current one
        ValueError: if the model has not been fitted yet
    """
    count = X_pression.shape[0]+(self.transform_input.shape[0]
                                 if self.transform_input is not None and doAdd else 0)
    print('looking at '+str(count)+' samples.')
    if self.fit_input is None:
      raise ValueError("no fit data available, need to run fit or addToFit first")
    transform_input = check_Xpression(X_pression, self.gene_file)
    if annotations is not None:
      if len(annotations) != len(transform_input) or list(transform_input.index) != list(annotations.index):
        raise ValueError("annotations do not match X_pression")
    else:
      # create fake annotations
      annotations = pd.DataFrame(index=X_pression.index,
                                 columns=['cell_type',
                                          'disease_type', 'tissue_type'],
                                 data=np.zeros((len(X_pression), 3))+self.number_of_datasets)
    if self.transform_input is None or not doAdd:
      # this is the first time we run it.
      print('creating a transform input..')
      self.common_genes = transform_input.columns
      self.fit_input = self.fit_input[self.common_genes]
      self.transform_input = transform_input
      self.transform_annotations = annotations
    else:
      # check if same genes
      print('adding to transform..')
      if set(transform_input.columns) != set(self.common_genes):
        raise ValueError("transform_input's genes do not match common_genes")
      # add annotations together
      self.transform_annotations = self.transform_annotations.append(annotations)
      # add transform_input together
      self.transform_input = self.transform_input.append(transform_input)
    self.number_of_datasets +=1
    if dotransform:
      return self.transform(only_transform=True)


  def transform(self, X_pression=None, annotations=None, only_transform=False):
    """transform the cell type for each sample in X_pression

    Args:
        X_pression ([type], optional): [description]. Defaults to None.
        annotations ([type], optional): [description]. Defaults to None.

    Raises:
        ValueError: [description]
    """
    if X_pression is not None:
      self.addTotransform(X_pression, annotations, dotransform=False, doAdd=False)
    elif self.transform_input is None:
      raise ValueError("no transform Expression data provided")
    # mean center the dataframe
    transform_input = self.transform_input.sub(self.transform_input.mean(axis=1), axis=0)
    # dimensionality reduction
    print('reducing dimensionality...')
    pca = PCA(**self.pca_kwargs) if not self.low_mem else IncrementalPCA(**self.pca_kwargs)
    pca_reduced = pca.fit_transform(transform_input)
    # clustering: doing SNN on the reduced data
    print('clustering..')
    self.transform_clusters = snn.SNN(**self.snn_kwargs).fit_predict(pca_reduced,
                                              sample_weight=None)
    if self.make_plots:
      # plotting
      plot.scatter(umap.UMAP(
        **self.umap_kwargs).fit_transform(np.vstack([self.fit_reduced, pca_reduced])), 
        xname="UMAP1", yname="UMAP2", colors=list(self.fit_clusters)+list(self.transform_clusters+len(set(self.fit_clusters))),
        labels=["fit_C"+str(i) for i in self.fit_clusters]+["transform_C"+str(i) for i in self.transform_clusters],
        title="SNN clusters", radi=.1, importance=[0]*len(self.fit_reduced) + [1]*len(pca_reduced))
    # do differential expression between clusters and getting the top K most expressed genes
    print('doing differential expression analysis on the clusters..')
    if len(set(self.transform_clusters)) < 2:
      raise ValueError("only one cluster found, no differential expression, try changing the parameters...")
    differential_genes = runDiffExprOnCluster(self.transform_input, self.transform_clusters)
    # need enough genes to be significant
    if len(differential_genes) < self.topKGenes:
      raise ValueError("not enough differentially expressed genes found, try changing the parameters..")
    # combining both ranks
    overlap = len(set(differential_genes.index[:self.topKGenes]) & 
      set(self.differential_genes_input.index[:self.topKGenes]))/self.topKGenes
    print("there is "+str(overlap)+" overlap between the fit and transform dataset in their most variable genes")
    # merge ranks
    self.differential_genes_names = []
    for i in range(self.topKGenes*2):
      if i%2==0:
        self.differential_genes_names.append(self.differential_genes_input.index[i//2])
      else:
        self.differential_genes_names.append(differential_genes.index[i//2])
    # removing cluster averages to samples clusters
    # TODO: take care of outlier cluster when outlier is authorized
    centered_fit_input = pd.concat([self.fit_input.loc[self.fit_clusters==val]\
      - self.fit_input.loc[self.fit_clusters==val].mean(axis=0) for val in set(self.fit_clusters)])
    centered_transform_input = pd.concat([self.transform_input.loc[self.transform_clusters == val]\
      - self.transform_input.loc[self.transform_clusters==val].mean(axis=0) for val in set(self.transform_clusters)])
    # doing cPCA on the dataset
    print('doing cPCA..')
    # TODO: try the automated version, (select the best alpha above 1?)
    cpca_loadings = CPCA(standardize=False, n_components=self.cpca_ncomp, low_memory=self.low_mem).fit(
      background=centered_transform_input, foreground=centered_fit_input, preprocess_with_pca_dim=centered_fit_input.shape[1]
      ).transform(only_loadings=True, return_alphas=False, alpha_selection = 'manual', **self.cpca_kwargs).T
    # regress out the cPCA components from the data
    print('regressing out the cPCA components..')
    # take the residuals of the linear regression of fit_input with the cpca_loadings
    del centered_transform_input, centered_fit_input
    transformed_fit = self.fit_input - LinearRegression(fit_intercept=False).fit(
      cpca_loadings, self.fit_input.T).predict(cpca_loadings).T
    transformed_transform = self.transform_input - LinearRegression(fit_intercept=False).fit(
      cpca_loadings, self.transform_input.T).predict(cpca_loadings).T
    
    # TODO?: how come we take the top K genes from a transformed matrix where we removed key axes of variances?
    # TODO?: in Allie's version, it was using different Ks for the two datasets. this tool only uses one
    varsubset = np.array([1 if i in self.differential_genes_names else 0 for i in self.transform_input.columns]).astype(bool)
    self.corrected, self.mnn_pairs, self.other  = mnnpy.mnn_correct(transformed_fit.values, 
                      transformed_transform.values,
                      var_index=list(range(len(transformed_fit.columns))),
                      **self.mnn_kwargs)
    import pdb; pdb.set_trace()
    del transformed_fit, transformed_transform
    self.corrected = pd.DataFrame(self.corrected, index=list(self.fit_input.index)+list(self.transform_input.index),
      columns=self.fit_input.columns)

    corrected_fit = self.corrected.iloc[:len(self.fit_input)]
    corrected_transform = self.corrected.iloc[len(self.fit_input):]
    self.mnn_pairs = mnn_pairs[-1]
    print("done")
    if only_transform:
      return corrected_transform
    else:
      return corrected_transform, corrected_fit, mnn_pairs


  def fit_transform(self, fit_X_pression=None, fit_annotations=None, 
    transform_X_pression=None, transform_annotations=None, only_transform=False):
    """fit_transform the data and transform the data.

    Args:
        fit_X_pression (pandas.DataFrame): the expression data to fit the model.
        fit_annotations (pandas.DataFrame): the annotations to fit the model.
        transform_X_pression (pandas.DataFrame): the expression data to transform.
        transform_annotations (pandas.DataFrame): the annotations to transform.
        only_transform (bool): if True, only transform the data.

    Returns:
        pandas.DataFrame: the transformed data.
    """
    self.fit(fit_X_pression, fit_annotations)
    return self.transform(transform_X_pression, transform_annotations)


  def save(self, folder, asData=False):
    """save the model to a folder

    Args:
      folder (str): folder to save the model
      asData (bool): if True, save the model as a dataframe, otherwise save it as a pickle file

    """
    # save the model
    if not os.path.exists(folder):
      os.makedirs(folder)
    if not asData:
      with open(os.path.join(folder, 'model.pkl'), 'wb') as f:
        pickle.dump(self, f)
      # save the data
    else:
      if not os.path.exists(os.path.join(folder, 'data')):
        os.makedirs(os.path.join(folder, 'data'))
      if self.fit_input is not None:
        self.fit_input.to_csv(os.path.join(folder, 'data', 'fit_input.csv'), index=None)
        self.fit_annotations.to_csv(os.path.join(folder, 'data', 'fit_annotations.csv'), index=None)
        self.fit_reduced.to_csv(os.path.join(folder, 'data', 'fit_reduced.csv'), index=None)
        h.listToFile(self.fit_clusters, os.path.join(folder, 'data', 'fit_clusters.csv'), index=None)
        self.differential_genes_input.to_csv(os.path.join(
          folder, 'data', 'differential_genes_input.csv'))
        h.listToFile(self.common_genes, os.path.join(folder, 'data', 'common_genes.csv'))
      if self.transform_input is not None:
        self.transform_input.to_csv(os.path.join(folder, 'data', 'transform_input.csv'), index=None)
        self.transform_annotations.to_csv(os.path.join(folder, 'data', 'transform_annotations.csv'), index=None)
        h.listToFile(self.transform_clusters, os.path.join(folder, 'data', 'transform_clusters.csv'))
        h.listToFile(self.differential_genes_names, os.path.join(folder, 'data', 'differential_genes_names.csv'))
        self.corrected.to_csv(os.path.join(folder, 'data', 'corrected.csv'), index=None)
        self.transform_reduced.to_csv(os.path.join(folder, 'data', 'transform_reduced.csv'), index=None)
        self.mnn_pairs.to_csv(os.path.join(folder, 'data', 'mnn_pairs.csv'), index=None)


  def load(self, folder):
    """load the model from a folder

    Args:
      folder (str): folder to load the model from
    """
    # if folder contains data folder
    if os.path.exists(os.path.join(folder, 'data')):
      # load the data
      if os.path.exists(os.path.join(folder, 'data', 'fit_input.csv')):
        self.fit_input = pd.read_csv(os.path.join(folder, 'data', 'fit_input.csv'))
        self.fit_reduced = pd.read_csv(os.path.join(folder, 'data', 'fit_reduced.csv'))
        self.fit_annotations = pd.read_csv(os.path.join(folder, 'data', 'fit_annotations.csv'))
        self.fit_clusters = h.fileToList(os.path.join(folder, 'data', 'fit_clusters.csv'))
        self.differential_genes_input = pd.read_csv(os.path.join(folder, 'data', 'differential_genes_input.csv'))
        self.common_genes = h.fileToList(os.path.join(folder, 'data', 'common_genes.csv'))
      if os.path.exists(os.path.join(folder, 'data', 'transform_input.csv')):
        self.transform_input = pd.read_csv(os.path.join(folder, 'data', 'transform_input.csv'))
        self.transform_reduced = pd.read_csv(os.path.join(folder, 'data', 'transform_reduced.csv'))
        self.transform_annotations = pd.read_csv(os.path.join(folder, 'data', 'transform_annotations.csv'))
        self.transform_clusters = h.fileToList(os.path.join(folder, 'data', 'transform_clusters.csv'))
        self.differential_genes_names = h.fileToList(os.path.join(folder, 'data', 'differential_genes_names.csv'))
        self.corrected = pd.read_csv(os.path.join(folder, 'data', 'corrected.csv'))
        self.mnn_pairs = pd.read_csv(os.path.join(folder, 'data', 'mnn_pairs.csv'))
    else:
      # load the model
      with open(os.path.join(folder, 'model.pkl'), 'rb') as f:
        model=pickle.load(f)
      self.__dict__.update(model.__dict__)

  def plot(self, onlyfit=False, onlytransform=False, corrected=True, umap_kwargs={},
           plot_kwargs={}, color_column="cell_type", show_clusts=True,annotations = None,
           smaller="fit"):
    """plot the model

    Args:
        onlyfit (bool, optional): if True, only plot the fit data. Defaults to False.
        onlytransform (bool, optional): if True, only plot the transform data. Defaults to False.
        corrected (bool, optional): if True, plot the corrected data. Defaults to True.
        umap_kwargs (dict, optional): kwargs for the umap plot. Defaults to {}.
        plot_kwargs (dict, optional): kwargs for the plot. Defaults to {}.
        color_column (str, optional): column to use for color. Defaults to "cell_type".
        show_clusts (bool, optional): if True, show the clusters. Defaults to True.
        annotations (pd.DataFrame, optional): annotations to use for the plot if none passed before. Defaults to None.
        smaller (str, optional): if "fit", plot the fit data smaller. If "transform", plot the transform data smaller.

    Raises:
        ValueError: model not fitted
    """
    # load the data based on availability
    if self.fit_input is None:
      raise ValueError('model not fitted yet')
    if onlyfit:
      if corrected:
        if self.corrected is None:
          print('no corrected fit data')
          data = self.fit_reduced
        else:
          data= self.corrected[:len(self.fit_input)]
      else:
        data = self.fit_reduced
      ann = self.fit_annotations
      clusts = ["fit_C"+str(i) for i in self.fit_clusters]
    elif onlytransform:
      if corrected:
        if self.corrected is None:
          print('no corrected transform data')
          data=self.transform_reduced
        else:
          self.corrected[len(self.fit_input):]
      else:
        data = self.transform_reduced
      ann = self.transform_annotations
      clusts = ["transform_C"+str(i) for i in self.transform_clusters]
    else:
      ann = self.fit_annotations
      clusts = ["fit_C"+str(i) for i in self.fit_clusters]
      if corrected:
        if self.corrected is None:
          print('no corrected data')
          data = self.fit_reduced
        else:
          data=self.corrected
          ann = ann.append(self.transform_annotations)
          clusts.extend(["transform_C"+str(i) for i in self.transform_clusters])
      else:
        if self.transform_reduced is None:
          data = self.fit_reduced
        else:
          data = self.fit_reduced.append(self.transform_reduced)
          ann = ann.append(self.transform_annotations)
          clusts.extend(["transform_C"+str(i) for i in self.transform_clusters])
    # doing UMAP
    self.umap_kwargs.update(umap_kwargs)
    umap_reduced=umap.UMAP(
        **umap_kwargs).fit_transform(data)
    if annotations is None:
      annotations = ann
    # plotting
    if 'labels' not in plot_kwargs and annotations is not None:
      # annotations to dict
      plot_kwargs['labels'] = {k: list(v) for k, v in annotations.T.iterrows()}
      plot_kwargs['labels'].update({'clusters': clusts})
    if 'colors' not in plot_kwargs:
      if show_clusts:
        col = { l: i for i, l in enumerate(set(clusts))}
        plot_kwargs.update({'colors':[col[x] for x in clusts]})
      else:
        col = { l: i for i, l in enumerate(set(annotations[color_column]))}
        plot_kwargs.update({'colors':[col[x] for x in annotations[color_column].tolist()]})
    # managing size
    if "importance" not in plot_kwargs:
      # 1 for all fit and 0 for all predict
      imp = np.zeros(len(data))
      if smaller == "fit":
        imp[:len(self.fit_input)]=1
      else:
        imp[len(self.fit_input):]=1
      plot_kwargs.update({'importance':imp})
    if 'xname' not in plot_kwargs:
      plot_kwargs.update({'xname':'UMAP1'})
    if 'yname' not in plot_kwargs:
      plot_kwargs.update({'yname':'UMAP2'})
    if 'title' not in plot_kwargs:
      plot_kwargs.update({'title':'Celligner plot'})
    if 'radi' not in plot_kwargs:
      plot_kwargs.update({'radi':0.1})
    plot.scatter(umap_reduced, **plot_kwargs)