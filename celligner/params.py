TISSUE_COLOR = {
    "engineered": "#bcdfbd",
    "fibroblast": "#9eAeAe",
    "other": "#A3969d",
    "skin": "#969696",
    "soft_tissue": "#cedb9c",  # put it closer to bone
    "sarcomatoid": "#cdcdbd",
    "unknown": "#bdbdbd",
    "NS": "#becdbd",
    "teratoma": "#252525",
    "germ_cell": "#c7c7c7",
    "embryo": "#7f7f7f",
    "bone": "#aec7e8",
    "lymphocyte": "#17becf",
    "plasma_cell": "#9edae5",
    "blood": "#1f77b4",
    "engineered_blood": "#2f87b4",
    "central_nervous_system": "#ff7f0e",
    "engineered_central_nervous_system": "#ff8f3f",
    "peripheral_nervous_system": "#ffbb78",
    "nerve": "#dbdb8d",
    "autonomic_ganglia": "#ebcb8d",
    "eye": "#bcbd22",
    "lung": "#d62728",
    "engineered_lung": "#ee2e3e",
    "upper_aerodigestive": "#ff9896",
    "esophagus": "#e7969c",
    "nasopharynx": "#f7b6d2",
    "oral": "#feceee",
    "parotid": "#fdbf6f",
    "stomach": "#e377c2",
    "gall_bladder": "#ff7f0e",
    "bile_duct": "#a55194",
    "engineered_bile_duct": "#a55194",
    "ampulla_of_vater": "#ad3184",
    "pancreas": "#e377c2",
    "liver": "#9467bd",
    "gastric": "#c49c94",
    "small_intestine": "#9e5e6e",
    "colon": "#8c564b",
    "ovary": "#2ca02c",
    "engineered_ovary": "#4eae4e",
    "uterus": "#98df8a",
    "cervix": "#5ab172",
    "breast": "#393b79",
    "engineered_breast": "#4e3e7e",
    "kidney": "#386cb0",
    "engineered_kidney": "#386cb0",
    "bladder": "#397cb9",
    "urinary_tract": "#b644dc",
    "prostate": "#637939",
    "engineered_prostate": "#6e7e3e",
    "testis": "#8c6d31",
    "thyroid": "#8f7e3e",
    "endocrine": "#bd9e39",
    "pineal": "#e7ba52",
    "adrenal": "#8ca252",
}

TISSUE_COLOR_R = {
    "central_nervous_system": "#f5899e",
    "engineered_central_nervous_system": "#f5899e",
    "teratoma": "#f5899e",
    "bone": "#9f55bb",
    "pancreas": "#b644dc",
    "soft_tissue": "#5fdb69",
    "skin": "#6c55e2",
    "liver": "#9c5e2b",
    "blood": "#da45bb",
    "lymphocyte": "#abd23f",
    "peripheral_nervous_system": "#73e03d",
    "ovary": "#56e79d",
    "engineered_ovary": "#56e79d",
    "adrenal": "#e13978",
    "adrenal_cortex": "#e13978",
    "upper_aerodigestive": "#5da134",
    "kidney": "#1f8fff",
    "engineered_kidney": "#1f8fff",
    "gastric": "#dfbc3a",
    "eye": "#349077",
    "nasopharynx": "#a9e082",
    "nerve": "#c44c90",
    "unknown": "#999999",
    "cervix": "#5ab172",
    "thyroid": "#d74829",
    "lung": "#51d5e0",
    "engineered_lung": "#51d5e0",
    "rhabdoid": "#d04850",
    "germ_cell": "#75dfbb",
    "embryo": "#75dfbb",
    "colorectal": "#96568e",
    "endocrine": "#d1d684",
    "bile_duct": "#c091e3",
    "pineal": "#949031",
    "thymus": "#659fd9",
    "mesothelioma": "#dc882d",
    "prostate": "#3870c9",
    "engineered_prostate": "#3870c9",
    "uterus": "#e491c1",
    "breast": "#45a132",
    "engineered_breast": "#45a132",
    "urinary_tract": "#e08571",
    "esophagus": "#6a6c2c",
    "fibroblast": "#d8ab6a",
    "plasma_cell": "#e6c241",
}

USEFUL_GENE_BIOTYPES = {
    "rRNA",
    "Mt_rRNA",
    "TR_V_gene",
    "IG_J_gene",
    "TEC",
    "TR_D_gene",
    "snoRNA",
    "snRNA",
    "IG_V_gene",
    "scRNA",
    "vault_RNA",
    "lncRNA",
    "miRNA",
    "scaRNA",
    "ribozyme",
    "protein_coding",
    "TR_C_gene",
    "sRNA",
    "TR_J_gene",
    "IG_D_gene",
    "IG_C_gene",
}

MIN_GENES = 5000

# differentially expressed genes with a rank better than this is in the cell line
# or tumor data are used to identify mutual nearest neighbors in the MNN alignment step
TOP_K_GENES = 1000


"""
@see https://umap-learn.readthedocs.io/en/latest/parameters.html

UMAP Parameters
----------
    n_neighbors: float(optional, default 15)
        The size of local neighborhood (in terms of number of neighboring
        sample points) used for manifold approximation. Larger values
        result in more global views of the manifold, while smaller
        values result in more local data being preserved. In general
        values should be in the range 2 to 100.
    n_components: int(optional, default 2)
        The dimension of the space to embed into. This defaults to 2 to
        provide easy visualization, but can reasonably be set to any
        integer value in the range 2 to 100.
    metric: string or function(optional, default 'euclidean')
        The metric to use to compute distances in high dimensional space.
        If a string is passed it must match a valid predefined metric. If
        a general metric is required a function that takes two 1d arrays and
        returns a float can be provided. For performance purposes it is
        required that this be a numba jit'd function. Valid string metrics
        include:
                * euclidean
                * manhattan
                * chebyshev
                * minkowski
                * canberra
                * braycurtis
                * mahalanobis
                * wminkowski
                * seuclidean
                * cosine
                * correlation
                * haversine
                * hamming
                * jaccard
                * dice
                * russelrao
                * kulsinski
                * ll_dirichlet
                * hellinger
                * rogerstanimoto
                * sokalmichener
                * sokalsneath
                * yule
        Metrics that take arguments(such as minkowski, mahalanobis etc.)
        can have arguments passed via the metric_kwds dictionary. At this
        time care must be taken and dictionary elements must be ordered
        appropriately; this will hopefully be fixed in the future.
    min_dist: float(optional, default 0.1)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points
        on the manifold are drawn closer together, while larger values will
        result on a more even dispersal of points. The value should be set
        relative to the ``spread`` value, which determines the scale at which
        embedded points will be spread out.
    low_memory: bool(optional, default True)
        For some datasets the nearest neighbor computation can consume a lot of
        memory. If you find that UMAP is failing due to memory constraints
        consider setting this option to True. This approach is more
        computationally expensive, but avoids excessive memory use.
"""
UMAP_PARAMS = {
    "n_neighbors": 10,
    "min_dist": 0.5,
    "metric": "euclidean",
    "n_components": 2,
}

# @see https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
PCA_PARAMS = {
    "n_components": 70,
}

# @see https://github.com/abidlabs/contrastive README
CPCA_PARAMS = {"alpha_value": 1}
# number of dimensions to reduce to in the embedding
CPCA_NCOMP = 4


"""
    :param datas: `numpy.ndarray` or class:`anndata.AnnData`
        Expression matrices or AnnData objects. Matrices should be shaped like n_obs * n_vars
        (n_cell * n_gene) and have consistent number of columns. AnnData objects should have same
        number of vars.
    :param var_index: `list` or `None`, optional (default: None)
        The index (list of str) of vars (genes). Necessary when using only a subset of vars to
        perform MNN correction, and should be supplied with var_subset. When datas are AnnData
        objects, var_index is ignored.
    :param var_subset: `list` or `None`, optional (default: None)
        The subset of vars (list of str) to be used when performing MNN correction. Typically, a
        list of highly variable genes (HVGs). When set to None, uses all vars.
    :param batch_key: `str`, optional (default: 'batch')
        The batch_key for AnnData.concatenate. Only valid when do_concatenate and supplying AnnData
        objects.
    :param index_unique: `str`, optional (default: '-')
        The index_unique for AnnData.concatenate. Only valid when do_concatenate and supplying
        AnnData objects.
    :param batch_categories: `list` or `None`, optional (default: None)
        The batch_categories for AnnData.concatenate. Only valid when do_concatenate and supplying
        AnnData objects.
    :param k1: `int`, optional (default: 20)
        Number of mutual nearest neighbors of 2 in 1.
    :param k2: `int`, optional (default: 20)
        Number of mutual nearest neighbors of 1 in 2.
    :param sigma: `float`, optional (default: 1)
        The bandwidth of the Gaussian smoothing kernel used to compute the correction vectors.
    :param cos_norm_in: `bool`, optional (default: True)
        Whether cosine normalization should be performed on the input data prior to calculating
        distances between cells.
    :param cos_norm_out: `bool`, optional (default: True)
        Whether cosine normalization should be performed prior to computing corrected expression
        values.
    :param svd_dim: `int` or `None`, optional (default: None)
        The number of dimensions to use for summarizing biological substructure within each batch.
        If set to None, biological components will not be removed from the correction vectors.
    :param var_adj: `bool`, optional (default: True)
        Whether to adjust variance of the correction vectors. Note this step takes most computing
        time.
    :param compute_angle: `bool`, optional (default: False)
        Whether to compute the angle between each cell’s correction vector and the biological
        subspace of the reference batch.
    :param mnn_order: `list` or `None`, optional (default: None)
        The order in which batches are to be corrected. When set to None, datas are corrected
        sequentially.
    :param svd_mode: `str`, optional (default: 'rsvd')
        One of 'svd', 'rsvd', and 'irlb'. 'svd' computes SVD using a non-randomized SVD-via-ID
        algorithm, while 'rsvd' uses a randomized version. 'irlb' performes truncated SVD by
        implicitly restarted Lanczos bidiagonalization (forked from https://github.com/airysen/irlbpy).
    :param do_concatenate: `bool`, optional (default: True)
        Whether to concatenate the corrected matrices or AnnData objects. Default is True.
    :param save_raw: `bool`, optional (default: False)
        Whether to save the original expression data in the .raw attribute of AnnData objects.
    :param n_jobs: `int` or `None`, optional (default: None)
        The number of jobs. When set to None, automatically uses the number of cores.
    :param kwargs: `dict` or `None`, optional (default: None)
        optional keyword arguments for irlb.
"""
# @see https://github.com/chriscainx/mnnpy/blob/master/mnnpy/mnn.py
# MNN_PARAMS = {
#  "var_adj": True,
#  "svd_mode": "rsvd",
#  "svd_dim": None,
# }

# For Mariona method (default)
MNN_PARAMS = {
    "k1": 5,
    "k2": 50,
    "cosine_norm": False,
    "fk": 5
}

# @see https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.louvain.html
LOUVAIN_PARAMS = {
    "resolution": 5,
}

# @see https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html
SC_NEIGH_PARAMS = {
    "n_neighbors": 20,
    "n_pcs": 70,
}
