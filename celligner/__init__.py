# Celligner
from celligner.params import *
from genepy.utils import helper as h
from genepy.utils import plot
from sklearn.decomposition import PCA
from genepy.rna import pyDESeq2
from snn import SNN
from contrastive import CPCA
import mnnpy

import pandas as pd
import numpy as np
import umap
import os
import pickle

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
  common_genes = set(X_pression.columns) & set(gene_file.ensembl_id)
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
  #TODO: add covariates
  if clust_covariates:
    if len(clust_covariates) != n_clusts:
      raise ValueError("number of covariates does not match number of clusters")
    design_matrix = clust_covariates
  # make a design matrix
  design_matrix = pd.DataFrame(index=expression.index,
                              data=[clustered==i for i in range(n_clusts)],
                              columns=['C'+str(i) for i in range(n_clusts)])
  formula = '~'+'+'.join(['C'+str(i) for i in range(n_clusts)])
  # creating the matrix
  data = expression.T
  data['gene_id'] = data.index
  data = data.reset_index(drop=True)
  # running DESeq2
  print("running DESeq2 on the samples")
  # Note: the differentially expressed genes don't make much sense themselves
  # they are computed by taking the difference across each cluster compared to the 
  # cluster 1
  # TODO: use the genes with the greatest variance?
  deseq = pyDESeq2.pyDESeq2(count_matrix=data, design_matrix=design_matrix,
                            design_formula=formula, gene_column="gene_id")
  deseq.run_deseq()
  deseq.get_deseq_result()
  r = deseq.deseq_result
  r.pvalue = np.nan_to_num(np.array(r.pvalue), 1)
  r.log2FoldChange = np.nan_to_num(np.array(r.log2FoldChange), 0)
  r = r.set_index("gene_id", drop=True)
  # TODO: why do we do it across each dataset separatly?
  # filtereing low p-values
  return r[(r.p_adj < pvalue_threshold) & (r.log2FoldChange > 0)]



class Celligner(object):
  def __init__(self, args, gene_file=None, onlyGenes=GENE_TYPE,
    ensemble_server="http://nov2020.archive.ensembl.org/biomart",
               umap_kwargs=UMAP_PARAMS, pca_kwargs=PCA_PARAMS,
               snn_kwargs=SNN_PARAMS, topKGenes=TOP_K_GENES, cpca_kwargs=CPCA_PARAMS, 
               cpca_ncomp=CPCA_NCOMP):
    """initialize Celligner object

    Args:
        args ([type]): [description]
        onlyGenes (str, optional): one of 'usefull', 'all', 'protein_coding'. Defaults to "usefull".
        gene_file (pd.Dataframe, optional): Needs to contain at least 15000 genes 
          and an "ensembl_id", columns. Defaults to None.
        ensemble_server (str, optional): [description]. Defaults to "http://nov2020.archive.ensembl.org/biomart".
        umap_kwargs (dict, optional): see umap_pamarameters.md or . Defaults to {}.
    """
    self.args = args
    if gene_file:
      self.gene_file = gene_file
    else:
      self.gene_file = h.generateGeneNames(ensemble_server=ensemble_server,
      useCache=True, cache_folder=".genenames/")
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
    self.number_of_datasets = 0

    self.clustered = None
    self.fit_input=None
    self.fit_reduced=None
    self.fit_clusters=None
    self.differential_genes_input=None
    self.differential_genes_names=None
    self.fit_annotations = None
    self.fit_umap=None
    self.fit_n_clusts=None
    self.predict_annotations = None
    self.predict_input=None
    self.predict_umap=None
    self.predict_clusters=None
    self.predict_reduced=None
    self.predict_n_clusts=None
    self.common_genes = None


      
  def addToFit(self, X_pression, annotations=None, dofit=True):
    """adds expression data to the fit dataframe

    Args:
        X_pression (pd.Dataframe): 
        annotations (pd.Dataframe, optional): [description]. Defaults to None.

    Raises:
        ValueError: if the expression matrix and annotations matrix do not have the same index
        ValueError: if the new expression matrix has different gene names than the current one
    """
    fit_input = check_Xpression(X_pression, self.gene_file)
    if annotations:
      if len(annotations) != len(fit_input) or fit_input.index.values != annotations.index.values:
        raise ValueError("annotations do not match X_pression")
    else:
      # create fake annotations
      annotations = pd.DataFrame(index=self.fit_input.index,
                                 columns=['cell_type',
                                          'disease_type', 'tissue_type'],
                                 data=np.zeros((len(self.fit_input), 3))+self.number_of_datasets)
      
    if self.common_genes is None:
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
          ensembl_id as columns and samplenames as index.
        annotations (pd.Dataframe, optional): sample annotations, for each sample, 
          needs to contain ['cell_type', 'disease type', 'tissue type']. 
          Defaults to None (will create an empty dataframe).
    """
    # check if X_pression is compatible with the model
    if X_pression:
      self.addToFit(X_pression, annotations, dofit=False)
    else:
      if self.fit_input is None:
        raise ValueError("no input provided")
  
    # mean center the dataframe
    fit_input = self.fit_input.sub(self.fit_input.mean(axis=1), axis=0)
    # dimensionality reduction
    print('reducing dimensionality...')
    self.fit_reduced = PCA(**self.pca_kwargs).fit_transform(fit_input) if self.doPCA else fit_input
    # clustering: doing SNN on the reduced data
    # TODO: try on the full data?
    # TODO: try DBSCAN or spectral clustering instead? (we might want more outliers)
    print('clustering...')
    self.fit_clusters = SNN(**self.snn_kwargs).fit_predict(self.fit_reduced,
        sample_weight=None)
    # do differential expression between clusters and getting the top K most expressed genes
    self.fit_n_clusts = len(set(self.fit_clusters))
    print('doing differential expression analysis on the clusters')
    if self.fit_n_clusts < 2:
      print("only one cluster found, no differential expression, using variance...")
      differential_genes = pd.DataFrame(index=self.fit_input.columns,
                                             columns=['var'],
                                             data=self.fit_input.var(axis=0))
    else:
      differential_genes = runDiffExprOnCluster(self.fit_input, self.fit_clusters)
      # need enough genes to be significant
      if len(differential_genes) < self.topKGenes:
        print("not enough differentially expressed genes found, using variance..")
        differential_genes = pd.DataFrame(index=self.fit_input.columns,
                                               columns=['var'],
                                               data=self.fit_input.var(axis=0))
      else:
        differential_genes = differential_genes.rename(columns={
          "log2FoldChange": "var"})
    self.differential_genes_input = differential_genes.sort_values(by='var', ascending=True)
    print('done')
    return self


  def addToPredict(self, X_pression, annotations=None, dopredict=True):
    """adds expression data to the predict dataframe

    Args:
        X_pression ([type]): [description]
        annotations ([type], optional): [description]. Defaults to None.

    Raises:
        ValueError: [description]
        ValueError: [description]
        ValueError: [description]
    """
    if self.fit_input is None:
      raise ValueError("no fit data available, need to run fit or addToFit first")
    predict_input = check_Xpression(X_pression, self.gene_file)
    if annotations:
      if len(annotations) != len(predict_input) or predict_input.index.values != annotations.index.values:
        raise ValueError("annotations do not match X_pression")
    else:
      # create fake annotations
      annotations = pd.DataFrame(index=self.predict_input.index,
                                 columns=['cell_type',
                                          'disease type', 'tissue type'],
                                 data=np.zeros((len(self.predict_input), 3))+self.number_of_datasets)
    if self.predict_input is None:
      print('creating a predict input..')
      self.predict_input = predict_input
      self.predict_annotations = annotations
    else:
      # check if same genes
      print('adding to predict..')
      if set(predict_input.columns) != set(self.common_genes):
        raise ValueError("predict_input's genes do not match common_genes")
      # add annotations together
      self.predict_annotations = self.predict_annotations.append(annotations)
      # add predict_input together
      self.predict_input = self.predict_input.append(predict_input)
    self.number_of_datasets +=1
    if dopredict:
      return self.predict(return_predict=True)


  def predict(self, X_pression=None, annotations=None, return_predict=False):
    """predict the cell type for each sample in X_pression

    Args:
        X_pression ([type], optional): [description]. Defaults to None.
        annotations ([type], optional): [description]. Defaults to None.

    Raises:
        ValueError: [description]
    """
    if X_pression:
      self.addToPredict(X_pression, annotations, dopredict=False)
    else:
      if self.predict_input is None:
        raise ValueError("no predict Expression data provided")
    # mean center the dataframe
    predict_input = self.predict_input.sub(self.predict_input.mean(axis=1), axis=0)
    # dimensionality reduction
    print('reducing dimensionality...')
    pca_reduced = PCA(**self.pca_kwargs).fit_transform(predict_input)
    # clustering: doing SNN on the reduced data
    print('clustering..')
    self.predict_clusters = SNN(**self.snn_kwargs).fit_predict(pca_reduced,
                                              sample_weight=None)
    # do differential expression between clusters and getting the top K most expressed genes
    self.predict_n_clusts = len(set(self.predict_clusters))
    print('doing differential expression analysis on the clusters..')
    if self.predict_n_clusts < 2:
      print("only one cluster found, no differential expression, using variance...")
      differential_genes = pd.DataFrame(index=self.predict_input.columns,
                                             columns=['var'],
                                             data=self.predict_input.var(axis=1))
    else:
      differential_genes = runDiffExprOnCluster(self.predict_input, self.predict_clusters)
      differential_genes.log2FoldChange = differential_genes.log2FoldChange.abs()
      # need enough genes to be significant
      if len(differential_genes) < self.topKGenes:
        print("not enough differentially expressed genes found, using variance..")
        differential_genes = pd.DataFrame(index=self.predict_input.columns,
                                               columns=['var'],
                                               data=self.predict_input.var(axis=1))
      else:
        differential_genes = differential_genes.rename(columns={
          "log2FoldChange": "var"})
    differential_genes = differential_genes.sort_values(by='var', ascending=True)
    # combining both ranks
    overlap = len(set(differential_genes.index[:self.topKGenes]) & 
      set(self.differential_genes_input.index[:self.topKGenes]))/self.topKGenes
    print("there is "+str(overlap)+" overlap between the fit and predict dataset in their most variable genes")
    # merge ranks
    self.differential_genes_names = []
    for i in range(self.topKGenes*2):
      if i%2==0:
        self.differential_genes_names.append(self.differential_genes_input.index[i//2])
      else:
        self.differential_genes_names.append(differential_genes.index[i//2])
    # removing cluster averages to samples clusters
    centered_fit_input = pd.DataFrame()
    #TODO move to fit
    centered_predict_input = pd.DataFrame()
    for val in range(self.fit_n_clusts):
      centered_fit_input.loc[self.fit_clusters==val] = self.fit_input.loc[self.fit_clusters==val]\
        - self.fit_input.loc[self.fit_clusters==val].mean(axis=0)
    for val in range(self.predict_n_clusts):
      centered_predict_input.loc[self.predict_clusters == val] = self.predict_input.loc[self.predict_clusters == val]\
        - self.predict_input.loc[self.predict_clusters==val].mean(axis=0)
    
    # doing cPCA on the dataset
    print('doing cPCA..')
    # TODO: in the examples of cPCA only the first 1000 highest dispersed genes in the foreground dataset are used
    # TODO: try to increase values of alpha so that we do not penalize too much the background:
    # TODO: try the automated version, (select the best alpha above 1?)
    cpca_loadings = CPCA(standardize=False, n_components=self.cpca_ncomp).fit(
        background=centered_predict_input, foreground=centered_fit_input).predict(
        only_loading=True, return_alphas=False, alpha_selection = 'manual', **self.cpca_kwargs)

    # we now remove the first n_components from the predict_input and fit_input
    # TODO: in Allie's version, we use the top eigenvectors to regress them from the matrix
    # take the residuals of the linear regression of fit_input with the cpca_loadings

     
    
    # TODO: in Allie's version, it was using different Ks for the two datasets. this tool only uses one
    
    # TODO: can we just do mnn and do almost as good?
    # TODO: how come we take the top K genes from a transformed matrix where we removed key axes of variances?
    # TODO: in Allie's version MNN vectors are applied to the entire cluster

    # TODO: in Allie's version each of the two dataset gets aligned by taking half of the vector

    # TODO: in Allie's version, the MNN vectors are computed with different values of K for the two datasets
    self.corrected = mnnpy.mnn_correct(transformed_fit,
                                  transformed_predict,
                                  **self.mnn_kwargs)
    corrected_fit = self.corrected.loc[:len(transformed_fit)]
    corrected_predict = self.corrected.loc[len(transformed_fit):]
    # TODO: recompute MNN vectors with the corrected data and apply a tricube weighting for the neighboring lines
    print("done")
    if return_predict:
      return corrected_predict
    else:
      return corrected_predict, corrected_fit


  def fit_predict(self, fit_X_pression=None, fit_annotations=None, 
    predict_X_pression=None, predict_annotations=None, return_predict=False):
    self.fit(fit_X_pression, fit_annotations)
    return self.predict(predict_X_pression, predict_annotations)

  def score(self): 
    # do GSEA on the cPCAs
    # 

  
  def save(self, folder, asData=False):
    """save the model to a folder
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
        h.listToFile(self.fit_clusters, os.path.join(folder, 'data', 'fit_clusters.csv'), index=None)
        self.differential_genes_input.to_csv(os.path.join(
          folder, 'data', 'differential_genes_input.csv'))
      if self.predict_input is not None:
        self.predict_input.to_csv(os.path.join(folder, 'data', 'predict_input.csv'), index=None)
        self.predict_annotations.to_csv(os.path.join(folder, 'data', 'predict_annotations.csv'), index=None)
        h.listToFile(self.predict_clusters, os.path.join(folder, 'data', 'predict_clusters.csv'))
        h.listToFile(self.differential_genes_names, os.path.join(folder, 'data', 'differential_genes_names.csv'))
        self.corrected.to_csv(os.path.join(folder, 'data', 'corrected.csv'), index=None)


  def load(self, folder):
    """load the model from a folder
    """
    # if folder contains data folder
    if os.path.exists(os.path.join(folder, 'data')):
      # load the data
      if os.path.exists(os.path.join(folder, 'data', 'fit_input.csv')):
        self.fit_input = pd.read_csv(os.path.join(folder, 'data', 'fit_input.csv'))
        self.fit_annotations = pd.read_csv(os.path.join(folder, 'data', 'fit_annotations.csv'))
        self.fit_clusters = h.fileToList(os.path.join(folder, 'data', 'fit_clusters.csv'))
        self.differential_genes_input = pd.read_csv(os.path.join(folder, 'data', 'differential_genes_input.csv'))
      if os.path.exists(os.path.join(folder, 'data', 'predict_input.csv')):
        self.predict_input = pd.read_csv(os.path.join(folder, 'data', 'predict_input.csv'))
        self.predict_annotations = pd.read_csv(os.path.join(folder, 'data', 'predict_annotations.csv'))
        self.predict_clusters = h.fileToList(os.path.join(folder, 'data', 'predict_clusters.csv'))
        self.differential_genes_names = h.fileToList(os.path.join(folder, 'data', 'differential_genes_names.csv'))
        self.corrected = pd.read_csv(os.path.join(folder, 'data', 'corrected.csv'))
    else:
      # load the model
      with open(os.path.join(folder, 'model.pkl'), 'rb') as f:
        model=pickle.load(f)
      self.__dict__.update(model.__dict__)

  def plot(self, onlyfit=False, onlypredict=False, corrected=True, umap_kwargs={},
           plot_kwargs={}, color_column="cell_type"):
    annotations = None
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
      annotations = self.fit_annotations
    elif onlypredict:
      if corrected:
        if self.corrected is None:
          print('no corrected predict data')
          data=self.predict_reduced
        else:
          self.corrected[len(self.fit_input):]
      else:
        data = self.predict_reduced
      annotations = self.predict_annotations
    else:
      annotations = self.fit_annotations
      if corrected:
        if self.corrected is None:
          print('no corrected data')
          data = self.fit_reduced
        else:
          data=self.corrected
          annotations = annotations.append(self.predict_annotations)
      else:
        if self.predict_reduced is None:
          data = self.fit_reduced
        else:
          data = self.fit_reduced.append(self.predict_reduced)
          annotations = annotations.append(self.predict_annotations)
    # doing UMAP
    umap_reduced=umap.UMAP(
        **umap_kwargs).fit_transform(data)
    # plotting
    if 'labels' not in plot_kwargs and annotations is not None:
      plot_kwargs.update({'labels':annotations})
    if 'colors' not in plot_kwargs:
      l = list(set(annotations[color_column]))
      col = { l[i]:i for i in range(len(l))}
      plot_kwargs.update({'colors':[col[x] for x in annotations[color_column].tolist()]})
    if 'xname' not in plot_kwargs:
      plot_kwargs.update({'xname':'UMAP1'})
    if 'yname' not in plot_kwargs:
      plot_kwargs.update({'yname':'UMAP2'})
    if 'title' not in plot_kwargs:
      plot_kwargs.update({'title':'Celligner plot'})
    plot.scatter(umap_reduced, **plot_kwargs)

