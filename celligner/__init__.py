from tabnanny import check
from celligner.params import *
from celligner import limma

from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.linear_model import LinearRegression
import sklearn.metrics as metrics
import umap.umap_ as umap

import scanpy as sc
from anndata import AnnData

import os
import pickle
import gc

import pandas as pd
import numpy as np

#from contrastive import CPCA
import mnnpy


class Celligner(object):
    def __init__(
        self,
        umap_kwargs=UMAP_PARAMS,
        pca_kwargs=PCA_PARAMS,
        neighbors_kwargs=SC_NEIGH_PARAMS,
        topKGenes=TOP_K_GENES,
        cpca_kwargs=CPCA_PARAMS,
        cpca_ncomp=CPCA_NCOMP,
        mnn_kwargs=MNN_PARAMS,
        low_mem=False,
        louvain_kwargs=LOUVAIN_PARAMS,
        mnn_method="mnn_marioni"
    ):
        """Initialize Celligner object

        Args:
            umap_kwargs (dict, optional): see params.py . 
            pca_kwargs (dict, optional): see see params.py . Defaults to {}.
            topKGenes (int, optional): see params.py. Defaults to 1000.
            cpca_kwargs (dict, optional): see see params.py . Currently unused.
            cpca_ncomp (int, optional): see params.py. Defaults to 4.
            mnn_kwargs (dict, optional): see params.py . 
            low_mem (bool, optional): adviced if you have less than 32Gb of RAM. Defaults to False.
            louvain_kwargs (dict, optional): see params.py . Defaults to {}.
            neighbors_kwargs (dict, optional): see params.py . Defaults to {}.
            mnn_method (str, optional): Only default "mnn_marioni" supported right now.
        """
        self.umap_kwargs = umap_kwargs
        self.pca_kwargs = pca_kwargs
        self.topKGenes = topKGenes
        self.cpca_kwargs = cpca_kwargs
        self.cpca_ncomp = cpca_ncomp
        self.mnn_kwargs = mnn_kwargs
        self.low_mem = low_mem
        self.louvain_kwargs = louvain_kwargs
        self.neighbors_kwargs = neighbors_kwargs
        self.mnn_method = mnn_method

        self.ref_input = None
        self.ref_annotations = None
        self.ref_clusters = None
        self.ref_de_genes = None
        
        self.target_input = None
        self.target_annotations = None
        self.target_clusters = None
        self.target_de_genes = None

        self.de_genes = None
        self.cpca_loadings = None
        self.target_corrected = None
        self.combined_output = None
        
        self.umap_reduced = None
        self.output_clusters = None


    def __checkExpression(self, expression, annot=None, check_genes=False):
        """
        Checks gene overlap for reference and target, checks for NaNs, then mean center the expression dataframe

        Args:
            expression (pd.Dataframe): expression data as samples x genes
            annot (pd.DataFrame): annotations for expression data samples (rows)
            gene_table (pd.Dataframe): gene dataframe with an ensembl_gene_id column

        Raises:
            
            ValueError: if some genes from the reference dataset are missing in the target dataset
            ValueError: if the expression matrix contains nan values
            ValueError: if the annotation matrix does not match the expression dataset samples

        Returns:
            (pd.Dataframe): the expression matrix
            (pd.Dataframe): the annotation matrix
        """
        
        # Check gene overlap if this is the target dataset (called from transform)
        if check_genes:
            if expression.loc[:, expression.columns.isin(self.common_genes)].shape[1] < len(self.common_genes):
                raise ValueError("Some genes from reference dataset not found in target dataset")
            expression = expression.loc[:, self.common_genes].astype(float)
        
        # Raise issue if there are any NaNs in the expression dataframe
        if expression.isnull().values.any():
            raise ValueError("Expression contains NaNs")

        # Mean center the expression dataframe
        expression = expression.sub(expression.mean(0), 1)

        # Check annotations
        if annot is not None:
            if len(annot) != expression.shape[0] or list(expression.index) != list(annot.index):
                raise ValueError("Annotations do not match expression dataframe")
        else:
            # Create empty annotations dataframe
            annot = pd.DataFrame(
                index=expression.index,
                columns=["cell_type", "disease_type", "tissue_type"],
                data=np.zeros((len(expression), 3)),
            )
        
        return expression, annot


    def __cluster(self, expression):
        """
        Cluster expression in 70-dimensional PCA space using a shared nearest neighbor based method

        Args:
            expression (pd.Dataframe): expression data as samples x genes

        Returns:
            (list): cluster label for each sample
        """
        # Create anndata object
        adata = AnnData(expression, dtype='float64')

        # Get PCs
        print("Doing PCA...")
        sc.tl.pca(adata, n_comps=70, zero_center=True)

        print("Computing neighbors...")
        # Find shared nearest neighbors (SNN) in PCA space
        # TODO? a bit different from R's version. ScanPy and Seurat differ in their implementation.
        sc.pp.neighbors(adata, knn=True, use_rep='X_pca', **self.neighbors_kwargs)
        
        print("Clustering...")
        sc.tl.louvain(adata, **self.louvain_kwargs)
        fit_clusters = adata.obs["louvain"].values.astype(int)
        
        del adata
        gc.collect()

        return fit_clusters


    def __runDiffExprOnClusters(self, expression, clusters):
        """
        Runs limma (R) on the clustered data.

        Args:
            expression (pd.Dataframe): expression data
            clusters (list): the cluster labels (per sample)

        Returns:
            (pd.Dataframe): limmapy results
        """

        n_clusts = len(set(clusters))
        print("Running differential expression on " + str(n_clusts) + " clusters")
        clusts = set(clusters) - set([-1])
        
        # make a design matrix
        design_matrix = pd.DataFrame(
            index=expression.index,
            data=np.array([clusters == i for i in clusts]).T,
            columns=["C" + str(i) + "C" for i in clusts],
        )
        design_matrix.index = design_matrix.index.astype(str).str.replace("-", ".")
        design_matrix = design_matrix[design_matrix.sum(1) > 0]
        
        # creating the matrix
        data = expression.T
        data = data[data.columns[clusters != -1].tolist()]
        
        # running limmapy
        print("Running limmapy on the samples")
        res = (
            limma.limmapy()
            .lmFit(data, design_matrix)
            .eBayes(trend=False)
            .topTable(number=len(data))
            .iloc[:, len(clusts) :]
        )
        return res.sort_values(by="F", ascending=False)
    

    def __runCPCA(self, centered_ref_input, centered_target_input):
        """
        Do cPCA on the centered reference and target expression datasets

        Args:
            centered_ref_input (pd.DataFrame): reference expression matrix where the cluster mean has been subtracted
            centered_target_input (pd.DataFrame): target expression matrix where the cluster mean has been subtracted

        """
        target_cov = centered_target_input.cov()
        ref_cov = centered_ref_input.cov()
        if not self.low_mem:
            pca = PCA(self.cpca_ncomp, svd_solver="randomized", copy=False)
        else: 
            pca = IncrementalPCA(self.cpca_ncomp, copy=False, batch_size=1000)
        
        return pca.fit(target_cov - ref_cov).components_


    def fit(self, ref_expr, ref_annot=None):
        """
        Fit the model to the reference expression dataset

        Args:
            ref_expr (pd.Dataframe): reference expression matrix of samples (rows) by genes (columns), 
                where genes are ensembl gene IDs. Data should be log2(X+1) TPM data. 
                In the standard Celligner pipeline this the cell line data.
            ref_annot (pd.Dataframe, optional): sample annotations for the reference dataframe,
                needs to contain ['cell_type', 'disease_type', 'tissue_type'].
                Defaults to None (will create an empty dataframe).

        Raises:
                ValueError: if only 1 cluster is found in the PCs of the expression
        """
        
        self.common_genes = list(ref_expr.columns)

        print("Using " + str(len(self.common_genes)) + " genes in reference expression matrix")

        self.ref_input, self.ref_annotations = self.__checkExpression(ref_expr, ref_annot)
        
        # Cluster and find differential expression for reference data
        self.ref_clusters = self.__cluster(self.ref_input)
        if len(set(self.ref_clusters)) < 2:
            raise ValueError("Only one cluster found in reference data, no differential expression possible")
        self.ref_de_genes = self.__runDiffExprOnClusters(self.ref_input, self.ref_clusters)

        return self


    def transform(self, target_expr=None, target_annot=None, compute_cPCs=True):
        """Align samples in the target dataset to samples in the reference dataset

        Args:
            target_expr (pd.Dataframe, optional): target expression matrix of samples (rows) by genes (columns), 
                where genes are ensembl gene IDs. Data should be log2(X+1) TPM data.
                In the standard Celligner pipeline this the tumor data (TCGA). 
                Set to None if re-running transform with new reference data.
            target_annot (pd.Dataframe, optional): sample annotations for the target dataframe,
                needs to contain ['cell_type', 'disease_type', 'tissue_type'].
                Defaults to None (will create an empty dataframe).
            compute_cPCs (bool, optional): if True, compute cPCs from the fitted reference and target expression. Defaults to True.

        Raises:
            ValueError: if compute_cPCs is True but there is no reference input (fit has not been run)
            ValueError: if compute_cPCs is False but there are no previously computed cPCs available
            ValueError: if compute_target_fit is False but there are no previously computed DE genes available for the target dataset
            ValueError: if compute_target_fit is True but compute_cPCs is false (there is no use case for this)
            ValueError: if there are not enough clusters to compute DE genes for the target dataset
        """

        if self.ref_input is None and compute_cPCs:
            raise ValueError("Need previously fitted reference dataset to compute cPCs, run fit() first")

        if not compute_cPCs and self.cpca_loadings is None:
            raise ValueError("Transform needs to be run with compute_cPCs==True at least once")

        if not target_expr is None and self.target_de_genes is None:
            raise ValueError("Transform needs to be run with a target expression dataset at least once")

        if target_expr is not None:

            self.target_input, self.target_annotations = self.__checkExpression(target_expr, target_annot, check_genes=True)

            # Cluster and find differential expression for target data
            self.target_clusters = self.__cluster(self.target_input)
            if len(set(self.target_clusters)) < 2:
                raise ValueError("Only one cluster found in reference data, no differential expression possible")
            self.target_de_genes = self.__runDiffExprOnClusters(self.target_input, self.target_clusters)

            # Union of the top 1000 differentially expressed genes in each dataset
            self.de_genes = list(set(self.ref_de_genes[:self.topKGenes].index) | 
                                set(self.target_de_genes[:self.topKGenes].index))

        else:
            print("No target expression provided, using previously fit target dataset")

        if compute_cPCs:

            # Subtract cluster average from cluster samples
            centered_ref_input = pd.concat(
                [
                    self.ref_input.loc[self.ref_clusters == val] - self.ref_input.loc[self.ref_clusters == val].mean(axis=0)
                    for val in set(self.ref_clusters)
                ]
            ).loc[self.ref_input.index]
            
            centered_target_input = pd.concat(
                [
                    self.target_input.loc[self.target_clusters == val] - self.target_input.loc[self.target_clusters == val].mean(axis=0)
                    for val in set(self.target_clusters)
                ]
            ).loc[self.target_input.index]
            
            # Compute contrastive PCs
            print("doing cPCA..")
            self.cpca_loadings = self.__runCPCA(centered_ref_input, centered_target_input)

            del centered_ref_input, centered_target_input
            gc.collect()

            print("Regressing top cPCs out of reference dataset..")
             # Take the residuals of the linear regression of ref_input with the cpca_loadings
            transformed_ref = (self.ref_input - 
                LinearRegression(fit_intercept=False)
                    .fit(self.cpca_loadings.T, self.ref_input.T)
                    .predict(self.cpca_loadings.T)
                    .T
            )

        else:
            self.target_input, self.target_annotations = self.__checkExpression(target_expr, target_annot, check_genes=True)
            transformed_ref = self.ref_input
        
        # Only need to regress out of target dataset if using previously computed cPCs
        print("Regressing top cPCs out of target dataset..")
        transformed_target= (self.target_input - 
            LinearRegression(fit_intercept=False)
                .fit(self.cpca_loadings.T, self.target_input.T)
                .predict(self.cpca_loadings.T)
                .T
        )

        # Do MNN (dropped Marioni version for now for simplicity/getting this running)
        varsubset = np.array([1 if i in self.de_genes else 0 for i in self.target_input.columns]).astype(bool)

        # if self.mnn_method == "mnn_marioni":
        print("Doing the MNN analysis using Marioni et al. method..")
        self.target_corrected, self.mnn_pairs = mnnpy.marioniCorrect(
            transformed_ref,
            transformed_target,
            var_index=list(range(len(self.ref_input.columns))),
            var_subset=varsubset,
            **self.mnn_kwargs,
        )
        
        # elif self.mnn_method == "mnn":
        #     print("doing the MNN analysis using scanPy MNN...")
        #     self.corrected, mnn_pairs, self.other = mnnpy.mnn_correct(
        #         transformed_ref.values,
        #         transformed_target.values,
        #         var_index=list(range(len(transformed_ref.columns))),
        #         varsubset=varsubset,
        #         **self.mnn_kwargs,
        #     )
        #     self.mnn_pairs = mnn_pairs[-1]
        #     self.corrected = pd.DataFrame(
        #         self.corrected[len(self.fit_input) :],
        #         index=list(self.transform_input.index),
        #         columns=self.transform_input.columns,
        #     )

        if compute_cPCs:
            self.combined_output =  pd.concat([self.target_corrected, transformed_ref])
            self.combined_annotations = pd.concat([self.target_annotations, self.ref_annotations])
        else: # Append at the end
            self.combined_output =  pd.concat([self.ref_input, self.target_corrected]) 
            self.combined_annotations = pd.concat([self.ref_annotations, self.target_annotations])
        
        print('Done')

        return self

    def computeMetricsForOutput(self, cl_ids=None, tumor_ids=None):
        """Compute UMAP embedding, clusters and tumor-cell line distance for combined output

        Args:
            cl_ids (list): cell line IDs for computing tumor-CL distance. Defaults to None, in which case the reference index is used.
            tumor_ids (list): tumor IDs for computing tumor-CL distance. Defaults to None, in which case the target index is used.

        Raises:
            ValueError: if there are no corrected expression matrices
        """
        if self.combined_output is None:
            raise ValueError("No corrected expression matrix found, run this function after transform()")

        print("Computing UMAP embedding...")
        # Compute UMAP embedding for results
        pca = PCA(**self.pca_kwargs)
        pcs = pca.fit_transform(self.combined_output)
        umap_reduced = umap.UMAP(**self.umap_kwargs, random_state=0).fit_transform(pcs)
        self.umap_reduced = pd.DataFrame(umap_reduced, index=self.combined_output.index, columns=['umap1','umap2'])
        # Add annotations to UMAP output
        self.umap_reduced = pd.concat([self.umap_reduced, pd.concat([self.target_annotations, self.ref_annotations])], axis=1)

        # Compute clusters in Celligner space
        self.output_clusters = self.__cluster(self.combined_output)

        print("Computing tumor-CL distance...")
        # Compute Euclidean distances between tumor samples and cell line samples in the combined 70 principal components
        pcs = pd.DataFrame(pcs, index=self.combined_output.index)
        if cl_ids is None: cl_ids = self.ref_input.index
        if tumor_ids is None: tumor_ids = self.target_input.index
        cl_pcs = pcs[pcs.index.isin(cl_ids)]
        tumor_pcs = pcs[pcs.index.isin(tumor_ids)]
        self.tumor_CL_dist = pd.DataFrame(metrics.pairwise_distances(tumor_pcs, cl_pcs), index=tumor_pcs.index, columns=cl_pcs.index)

        return self


    def makeNewReference(self):
        """Make a new reference dataset from the previously transformed fit+target datasets. Use for multi-dataset algignment.
        
        """
        self.ref_input = self.combined_output
        self.ref_annotations = self.combined_annotations

        return self
    
    
    def save(self, file_name):
        """save the model as a pickle file

        Args:
            file_name (str): name of file in which to save the model
        """
        # save the model
        with open(os.path.normpath(file_name), "wb") as f:
            pickle.dump(self, f)


    def load(self, file_name):
        """load the model from a pickle file

        Args:
            file_name (str): pickle file to load the model from
        """
        with open(os.path.normpath(file_name), "rb") as f:
            model = pickle.load(f)
            self.__dict__.update(model.__dict__)
        return self