# Celligner
from celligner.params import *
from genepy.utils import helper as h
from genepy.utils import plot
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.linear_model import LinearRegression
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
from contrastive import CPCA
import mnnpy
from celligner import limma

def check_Xpression(X_pression, gene_file):
  """

  Args:
    X_pression (pd.Dataframe): expression data
    gene_file (pd.Dataframe): gene file with an ensembl_gene_id column

  Raises:
    ValueError: if the number of genes in the expression matrix and the gene file do not match
    ValueError: if the expression matrix contains nan values

  Returns:
    (pd.Dataframe): the expression matrix
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


def runDiffExprOnCluster(expression, clustered, clust_covariates=None,):
  """
  Runs DESEQ2 on the clustered data.

  Args:
    expression (pd.Dataframe): expression data
    clustered (list): the clusters
    clust_covariates (pd.Dataframe, optional): covariates for the clustering. Defaults to None.

  Returns:
    (pd.Dataframe): limmapy results

  Raises:
    ValueError: if the number of genes in the expression matrix and the gene file do not match
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

class Celligner(object):
  def __init__(self, gene_file=None, onlyGenes=GENE_TYPE,
               ensemble_server="http://nov2020.archive.ensembl.org/biomart",
               umap_kwargs=UMAP_PARAMS, pca_kwargs=PCA_PARAMS,
               snn_kwargs=SNN_PARAMS, topKGenes=TOP_K_GENES, cpca_kwargs=CPCA_PARAMS, 
               cpca_ncomp=CPCA_NCOMP, mnn_kwargs=MNN_PARAMS, make_plots=False,
               low_mem=False, louvain_kwargs=LOUVAIN_PARAMS, method="mnn_marioni"): #scneigh_kwargs=SC_NEIGH_PARAMS
    """initialize Celligner object

    Args:
      onlyGenes (str, optional): one of 'usefull', 'all', 'protein_coding'. Defaults to "usefull".
      gene_file (pd.Dataframe, optional): Needs to contain at least 15000 genes 
        and an "ensembl_gene_id", columns. Defaults to None.
      ensemble_server (str, optional): the ensembl biomart server to map genes to. 
        Defaults to "http://nov2020.archive.ensembl.org/biomart".
      umap_kwargs (dict, optional): see params.py . Defaults to {}.
      pca_kwargs (dict, optional): see see params.py . Defaults to {}.
      snn_kwargs (dict, optional): see see params.py . Defaults to {}.
      topKGenes (int, optional): see params.py. Defaults to 1000.
      cpca_kwargs (dict, optional): see see params.py . Defaults to {}.
      cpca_ncomp (int, optional): see params.py. Defaults to 10.
      mnn_kwargs (dict, optional): see params.py . Defaults to {}.
      make_plots (bool, optional): whether to log multiple plots along the way. Defaults to False.
      low_mem (bool, optional): adviced if you have less than 32Gb of RAM. Defaults to False.
      louvain_kwargs (dict, optional): see params.py . Defaults to {}.
      method (str, optional): either "mnn_marioni" or "mnn". Defaults to "mnn_marioni".
    """
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
    self.louvain_kwargs = louvain_kwargs
    #self.scneigh_kwargs = scneigh_kwargs
    self.method = method

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
    self.pca_fit=None
    self.pca_transform=None
    self.common_genes = None
    self.cpca_loadings=None


  def addToFit(self, X_pression, annotations=None, dofit=True, doAdd=True):
    """adds expression data to the fit dataframe

    Args:
      X_pression (pd.Dataframe): expression data
      annotations (pd.Dataframe, optional): sample annotations. Defaults to None.

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
    elif doAdd: 
      return self.fit(_rerun=False)


  def fit(self, X_pression=None, annotations=None, _rerun=True):
    """fit the model using X_pression

    Args:
      X_pression (pd.Dataframe): contains the expression data as RSEM expected counts with 
        ensembl_gene_id as columns and samplenames as index.
      annotations (pd.Dataframe, optional): sample annotations, for each sample, 
        needs to contain ['cell_type', 'disease_type', 'tissue_type']. 
        Defaults to None (will create an empty dataframe).
      _rerun (bool, optional): whether to rerun the function entirely or not. Defaults to True.

    Raises:
        ValueError: if the expression matrix and annotations matrix do not have the same index
        ValueError: if the new expression matrix has different gene names than the current one
    """
    # check if X_pression is compatible with the model
    if X_pression is not None:
      self.addToFit(X_pression, annotations, dofit=False, doAdd=False)
    elif self.fit_input is None:
      raise ValueError("no input provided")
  
    # mean center the dataframe
    self.fit_input = self.fit_input.sub(self.fit_input.mean(axis=1), axis=0)
    # dimensionality reduction
    print('reducing dimensionality...')
    if _rerun:
      self.pca_fit = PCA(**self.pca_kwargs) if not self.low_mem else IncrementalPCA(**self.pca_kwargs)
      self.fit_reduced = self.pca_fit.fit_transform(self.fit_input)
    else:
      self.fit_reduced = self.pca_fit.transform(self.fit_input)
    # clustering: doing SNN on the reduced data
    print('clustering...')
    #anndata from df
    adata = AnnData(self.fit_reduced)
    neighbors(adata) # **scneigh_kwargs
    louvain(adata, **self.louvain_kwargs)
    self.fit_clusters = adata.obs['louvain'].values.astype(int)
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
    
    if len(set(self.fit_clusters)) < 2:
      raise ValueError("only one cluster found, no differential expression possible\
        try to change your parameters...")
    if _rerun:
      print('doing differential expression analysis on the clusters')
      self.differential_genes_input = runDiffExprOnCluster(
          self.fit_input, self.fit_clusters)
      # need enough genes to be significant
    if len(self.differential_genes_input[self.differential_genes_input.F>10]) < self.topKGenes:
      raise ValueError("not enough differentially expressed genes found..")
    print('done')
    return self

  def putAllToFit(self):
    """puts all the data to the fit dataframe"""
    self.fit_annotations = self.fit_annotations.append(self.transform_annotations)
    self.fit_input = self.fit_input.append(self.transform_input)
    return self.fit()


  def addTotransform(self, X_pression, annotations=None, dotransform=True, doAdd=True):
    """adds expression data to the transform dataframe

    Args:
      X_pression (pd.Dataframe): the expression data as RSEM expected counts
        with ensembl_gene_id as columns and samplenames as index.
      annotations (pd.Dataframe, optional): sample annotations, for each sample,
      dotransform (bool, optional): if True, will transform the data. Defaults to True.
      doAdd (bool, optional): if True, will add the data to the transform dataframe.

    Returns:
      (, optional): transform()'s output

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
    elif doAdd:
      return self.transform(_rerun=False)


  def transform(self, X_pression=None, annotations=None, only_transform=False, _rerun=True):
    """transform the cell type for each sample in X_pression

    Args:
      X_pression (pd.Dataframe, optional): expression dataframe. Defaults to None.
      annotations (pd.Dataframe, optional): annotations dataframe. Defaults to None.
      only_transform (bool, optional): if True, will only transform the dataframe.
      _rerun (bool, optional): if True, will rerun the PCA and SNN. Defaults to True.

    Raises:
      ValueError: if the model has not been fitted yet
      ValueError: if the expression matrix and annotations matrix do not have the same index
      ValueError: if the new expression matrix has different gene names than the current one
    """
    if X_pression is not None:
      self.addTotransform(X_pression, annotations, dotransform=False, doAdd=False)
    elif self.transform_input is None:
      raise ValueError("no transform Expression data provided")
    # mean center the dataframe
    self.transform_input = self.transform_input.sub(self.transform_input.mean(axis=1), axis=0)
    
    # dimensionality reduction
    print('reducing dimensionality...')
    if _rerun:
      self.pca_transform = PCA(**self.pca_kwargs) if not self.low_mem else IncrementalPCA(**self.pca_kwargs)
      self.pca_transform = self.pca_transform.fit_transform(self.transform_input)
    else:
      self.transform_reduced = self.pca_transform.transform(self.transform_input)
    
    if _rerun:
      # clustering: doing SNN on the reduced data
      print('clustering..')
      #anndata from df
      adata = AnnData(self.transform_reduced)
      neighbors(adata) # **scneigh_kwargs
      louvain(adata, **self.louvain_kwargs)
      self.transform_clusters = adata.obs['louvain'].values.astype(int)
      del adata
      #self.transform_clusters = snn.SNN(**self.snn_kwargs).fit_predict(self.transform_reduced,
      #                                          sample_weight=None)

    if self.make_plots:
      # plotting
      plot.scatter(umap.UMAP(
        **self.umap_kwargs).fit_transform(np.vstack([self.fit_reduced, self.transform_reduced])), 
        xname="UMAP1", yname="UMAP2", colors=list(self.fit_clusters)+list(self.transform_clusters+len(set(self.fit_clusters))),
        labels=["fit_C"+str(i) for i in self.fit_clusters]+["transform_C"+str(i) for i in self.transform_clusters],
        title="SNN clusters", radi=.1, importance=[0]*len(self.fit_reduced) + [1]*len(self.transform_reduced))
    
    # do differential expression between clusters and getting the top K most expressed genes
    print('doing differential expression analysis on the clusters..')
    if len(set(self.transform_clusters)) < 2:
      raise ValueError("only one cluster found, no differential expression, try changing the parameters...")
    
    if _rerun:
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
      self.cpca_loadings = CPCA(standardize=False, n_components=self.cpca_ncomp, low_memory=self.low_mem).fit(
        background=centered_transform_input, foreground=centered_fit_input, preprocess_with_pca_dim=centered_fit_input.shape[1]
        ).transform(only_loadings=True, return_alphas=False, alpha_selection = 'manual', **self.cpca_kwargs).T
      del centered_transform_input, centered_fit_input
    
    # regress out the cPCA components from the data
    print('regressing out the cPCA components..')
    # take the residuals of the linear regression of fit_input with the cpca_loadings
    transformed_fit = self.fit_input - LinearRegression(fit_intercept=False).fit(
      self.cpca_loadings, self.fit_input.T).predict(self.cpca_loadings).T
    transformed_transform = self.transform_input - LinearRegression(fit_intercept=False).fit(
      self.cpca_loadings, self.transform_input.T).predict(self.cpca_loadings).T
    
    # TODO?: how come we take the top K genes from a transformed matrix where we removed key axes of variances?
    # TODO?: in Allie's version, it was using different Ks for the two datasets. this tool only uses one
    varsubset = np.array([1 if i in self.differential_genes_names else 0 for i in self.transform_input.columns]).astype(bool)
    if self.method == 'mnn_marioni':
      print('doing the MNN analysis using Marioni et al. method..')
      self.corrected, self.mnn_pairs = mnnpy.marioniCorrect(transformed_fit,
        transformed_transform, var_index = list(range(len(transformed_fit.columns))), 
        var_subset=varsubset, **self.mnn_kwargs)
      #marioniCorrect(transformed_fit.values, 
      #  transformed_transform.values, var_index = list(range(len(transformed_fit.columns))),
      #  var_subset=varsubset, **self.mnn_kwargs)
    elif self.method == "mnn":
      print('doing the MNN analysis using scanPy MNN...')
      self.corrected, mnn_pairs, self.other  = mnnpy.mnn_correct(transformed_fit.values,
                        transformed_transform.values, 
                        var_index=list(range(len(transformed_fit.columns))),
                        varsubset=varsubset,
                        **self.mnn_kwargs)
      self.mnn_pairs = mnn_pairs[-1]
      self.corrected = pd.DataFrame(self.corrected, index=list(self.fit_input.index)+list(self.transform_input.index),
        columns=self.fit_input.columns)
      self.fit_input = self.corrected.iloc[:len(self.fit_input)]
      self.corrected = self.corrected.iloc[len(self.fit_input):]
    del transformed_fit, transformed_transform

    print("done")
    if self.make_plots:
      self.plot()
    if only_transform:
      return self.corrected
    else:
      return self.corrected, self.fit_input, self.mnn_pairs


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
        self.pca_fit.to_csv(os.path.join(folder, 'data', 'pca_transform.csv'), index=None)
      if self.transform_input is not None:
        self.transform_input.to_csv(os.path.join(folder, 'data', 'transform_input.csv'), index=None)
        self.transform_annotations.to_csv(os.path.join(folder, 'data', 'transform_annotations.csv'), index=None)
        h.listToFile(self.transform_clusters, os.path.join(folder, 'data', 'transform_clusters.csv'))
        h.listToFile(self.differential_genes_names, os.path.join(folder, 'data', 'differential_genes_names.csv'))
        self.corrected.to_csv(os.path.join(folder, 'data', 'corrected.csv'), index=None)
        self.transform_reduced.to_csv(os.path.join(folder, 'data', 'transform_reduced.csv'), index=None)
        self.mnn_pairs.to_csv(os.path.join(folder, 'data', 'mnn_pairs.csv'), index=None)
        self.pca_transform.to_csv(os.path.join(folder, 'data', 'pca_transform.csv'), index=None)
        self.cpca_loadings.to_csv(os.path.join(folder, 'data', 'cpca_loadings.csv'), index=None)


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
        self.pca_fit = pd.read_csv(os.path.join(folder, 'data', 'pca_transform.csv'))
      if os.path.exists(os.path.join(folder, 'data', 'transform_input.csv')):
        self.transform_input = pd.read_csv(os.path.join(folder, 'data', 'transform_input.csv'))
        self.transform_reduced = pd.read_csv(os.path.join(folder, 'data', 'transform_reduced.csv'))
        self.transform_annotations = pd.read_csv(os.path.join(folder, 'data', 'transform_annotations.csv'))
        self.transform_clusters = h.fileToList(os.path.join(folder, 'data', 'transform_clusters.csv'))
        self.differential_genes_names = h.fileToList(os.path.join(folder, 'data', 'differential_genes_names.csv'))
        self.corrected = pd.read_csv(os.path.join(folder, 'data', 'corrected.csv'))
        self.mnn_pairs = pd.read_csv(os.path.join(folder, 'data', 'mnn_pairs.csv'))
        self.pca_transform = pd.read_csv(os.path.join(folder, 'data', 'pca_transform.csv'))
        self.cpca_loadings = pd.read_csv(os.path.join(folder, 'data', 'cpca_loadings.csv'))
    else:
      # load the model
      with open(os.path.join(folder, 'model.pkl'), 'rb') as f:
        model=pickle.load(f)
      self.__dict__.update(model.__dict__)

  def plot(self, onlyfit=False, onlytransform=False, corrected=True, umap_kwargs={},
           plot_kwargs={}, color_column="cell_type", show_clusts=False,annotations = None,
           smaller="predict", rerun=True, colortable=None,):
    """plot the model

    Args:
      onlyfit (bool, optional): if True, only plot the fit data. Defaults to False.
      onlytransform (bool, optional): if True, only plot the transform data. Defaults to False.
      corrected (bool, optional): if True, plot the corrected data. Defaults to True.
      umap_kwargs (dict, optional): kwargs for the umap plot. Defaults to {}.
      plot_kwargs (dict, optional): kwargs for the plot. Defaults to {}.
      color_column (str, optional): column to use for color. Defaults to "cell_type".
      show_clusts (bool, optional): if True, show the clusters. Defaults to True.
      annotations (pd.DataFrame, optional): annotations to use for the plot if none passed before. Defaults to None. If None, use the fit annotations.
      smaller (str, optional): if "fit", plot the fit data smaller. If "transform", plot the transform data smaller. Defaults to "fit".
      rerun (bool, optional): if True, rerun the umap and plot. Defaults to True.
      colortable (str, optional): if not None, use this colortable else chooses viridis. Defaults to None.

    Raises:
      ValueError: model not fitted
    """
    # load the data based on availability
    if rerun or self.umap_reduced is None:
      if self.fit_input is None:
        raise ValueError('model not fitted yet')
      if onlyfit:
        data = self.fit_input
        ann = self.fit_annotations
        clusts = ["fit_C"+str(i) for i in self.fit_clusters]
      elif onlytransform:
        if corrected:
          if self.corrected is None:
            print('no corrected transform data')
            data=self.transform_input
          else:
            data=self.corrected
        else:
          data = self.transform_input
        ann = self.transform_annotations
        clusts = ["transform_C"+str(i) for i in self.transform_clusters]
      else:
        ann = self.fit_annotations
        clusts = ["fit_C"+str(i) for i in self.fit_clusters]
        if corrected:
          if self.corrected is None:
            print('no corrected data')
            data = self.fit_input
          else:
            data=self.fit_input.append(self.corrected)
            ann = ann.append(self.transform_annotations)
            clusts.extend(["transform_C"+str(i) for i in self.transform_clusters])
        else:
          if self.transform_input is None:
            data = self.fit_input
          else:
            data = self.fit_input.append(self.transform_input)
            ann = ann.append(self.transform_annotations)
            clusts.extend(["transform_C"+str(i) for i in self.transform_clusters])
      # doing UMAP
      self.umap_kwargs.update(umap_kwargs)
      print('reducing dimensionality...')
      pca = PCA(**self.pca_kwargs) if not self.low_mem else IncrementalPCA(**self.pca_kwargs)
      data = pca.fit_transform(data)
      umap_reduced=umap.UMAP(
          **umap_kwargs).fit_transform(data)
      if annotations is None:
        annotations = ann
      self.umap_reduced = umap_reduced
      self.annotations = annotations
      self.clusts = clusts
    # plotting
    if 'labels' not in plot_kwargs and self.annotations is not None:
      # annotations to dict
      plot_kwargs['labels'] = {k: list(v) for k, v in self.annotations.T.iterrows()}
      plot_kwargs['labels'].update({'clusters': self.clusts})
    if 'colors' not in plot_kwargs:
      if show_clusts:
        col = { l: i for i, l in enumerate(set(self.clusts))}
        plot_kwargs.update({'colors':[col[x] for x in self.clusts]})
      else:
        if colortable is None:
          col = { l: i for i, l in enumerate(set(self.annotations[color_column]))}
        else:
          col = colortable
          plot_kwargs.update({"colprovided": True})
        plot_kwargs.update({'colors':[col[x] for x in self.annotations[color_column].tolist()]})
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
    print('making plot...')
    plot.scatter(self.umap_reduced, **plot_kwargs)