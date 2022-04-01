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

# @see Celligner package
GENE_TYPE = "usefull"

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
    n_epochs: int(optional, default None)
        The number of training epochs to be used in optimizing the
        low dimensional embedding. Larger values result in more accurate
        embeddings. If None is specified a value will be selected based on
        the size of the input dataset(200 for large datasets, 500 for small).
    learning_rate: float(optional, default 1.0)
        The initial learning rate for the embedding optimization.
    init: string(optional, default 'spectral')
        How to initialize the low dimensional embedding. Options are:
                * 'spectral': use a spectral embedding of the fuzzy 1-skeleton
                * 'random': assign initial embedding positions at random.
                * A numpy array of initial embedding positions.
    min_dist: float(optional, default 0.1)
        The effective minimum distance between embedded points. Smaller values
        will result in a more clustered/clumped embedding where nearby points
        on the manifold are drawn closer together, while larger values will
        result on a more even dispersal of points. The value should be set
        relative to the ``spread`` value, which determines the scale at which
        embedded points will be spread out.
    spread: float(optional, default 1.0)
        The effective scale of embedded points. In combination with ``min_dist``
        this determines how clustered/clumped the embedded points are.
    low_memory: bool(optional, default True)
        For some datasets the nearest neighbor computation can consume a lot of
        memory. If you find that UMAP is failing due to memory constraints
        consider setting this option to True. This approach is more
        computationally expensive, but avoids excessive memory use.
    set_op_mix_ratio: float(optional, default 1.0)
        Interpolate between(fuzzy) union and intersection as the set operation
        used to combine local fuzzy simplicial sets to obtain a global fuzzy
        simplicial sets. Both fuzzy set operations use the product t-norm.
        The value of this parameter should be between 0.0 and 1.0; a value of
        1.0 will use a pure fuzzy union, while 0.0 will use a pure fuzzy
        intersection.
    local_connectivity: int(optional, default 1)
        The local connectivity required - - i.e. the number of nearest
        neighbors that should be assumed to be connected at a local level.
        The higher this value the more connected the manifold becomes
        locally. In practice this should be not more than the local intrinsic
        dimension of the manifold.
    repulsion_strength: float(optional, default 1.0)
        Weighting applied to negative samples in low dimensional embedding
        optimization. Values higher than one will result in greater weight
        being given to negative samples.
    negative_sample_rate: int(optional, default 5)
        The number of negative samples to select per positive sample
        in the optimization process. Increasing this value will result
        in greater repulsive force being applied, greater optimization
        cost, but slightly more accuracy.
    transform_queue_size: float(optional, default 4.0)
        For transform operations(embedding new points using a trained model_
        this will control how aggressively to search for nearest neighbors.
        Larger values will result in slower performance but more accurate
        nearest neighbor evaluation.
    a: float(optional, default None)
        More specific parameters controlling the embedding. If None these
        values are set automatically as determined by ``min_dist`` and ``spread``.
    b: float(optional, default None)
        More specific parameters controlling the embedding. If None these
        values are set automatically as determined by ``min_dist`` and ``spread``.
    random_state: int, RandomState instance or None, optional(default: None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.
    metric_kwds: dict(optional, default None)
        Arguments to pass on to the metric, such as the ``p`` value for
        Minkowski distance. If None then no arguments are passed on.
    angular_rp_forest: bool(optional, default False)
        Whether to use an angular random projection forest to initialise
        the approximate nearest neighbor search. This can be faster, but is
        mostly on useful for metric that use an angular style distance such
        as cosine, correlation etc. In the case of those metrics angular forests
        will be chosen automatically.
    target_n_neighbors: int(optional, default - 1)
        The number of nearest neighbors to use to construct the target simplcial
        set. If set to - 1 use the ``n_neighbors`` value.
    target_metric: string or callable(optional, default 'categorical')
        The metric used to measure distance for a target array is using supervised
        dimension reduction. By default this is 'categorical' which will measure
        distance in terms of whether categories match or are different. Furthermore,
        if semi-supervised is required target values of - 1 will be trated as
        unlabelled under the 'categorical' metric. If the target array takes
        continuous values(e.g. for a regression problem) then metric of 'l1'
        or 'l2' is probably more appropriate.
    target_metric_kwds: dict(optional, default None)
        Keyword argument to pass to the target metric when performing
        supervised dimension reduction. If None then no arguments are passed on.
    target_weight: float(optional, default 0.5)
        weighting factor between data topology and target topology. A value of
        0.0 weights entirely on data, a value of 1.0 weights entirely on target.
        The default of 0.5 balances the weighting equally between data and target.
    transform_seed: int(optional, default 42)
        Random seed used for the stochastic aspects of the transform operation.
        This ensures consistency in transform operations.
    verbose: bool(optional, default False)
        Controls verbosity of logging.
    tqdm_kwds: dict(optional, defaul None)
        Key word arguments to be used by the tqdm progress bar.
    unique: bool(optional, default False)
        Controls if the rows of your data should be uniqued before being
        embedded.  If you have more duplicates than you have n_neighbour
        you can have the identical data points lying in different regions of
        your space.  It also violates the definition of a metric.
        For to map from internal structures back to your data use the variable
        _unique_inverse_.
    densmap: bool(optional, default False)
        Specifies whether the density-augmented objective of densMAP
        should be used for optimization. Turning on this option generates
        an embedding where the local densities are encouraged to be correlated
        with those in the original space. Parameters below with the prefix 'dens'
        further control the behavior of this extension.
    dens_lambda: float(optional, default 2.0)
        Controls the regularization weight of the density correlation term
        in densMAP. Higher values prioritize density preservation over the
        UMAP objective, and vice versa for values closer to zero. Setting this
        parameter to zero is equivalent to running the original UMAP algorithm.
    dens_frac: float(optional, default 0.3)
        Controls the fraction of epochs(between 0 and 1) where the
        density-augmented objective is used in densMAP. The first
        (1 - dens_frac) fraction of epochs optimize the original UMAP objective
        before introducing the density correlation term.
    dens_var_shift: float(optional, default 0.1)
        A small constant added to the variance of local radii in the
        embedding when calculating the density correlation objective to
        prevent numerical instability from dividing by a small number
    output_dens: float(optional, default False)
        Determines whether the local radii of the final embedding(an inverse
        measure of local density) are computed and returned in addition to
        the embedding. If set to True, local radii of the original data
        are also included in the output for comparison; the output is a tuple
        (embedding, original local radii, embedding local radii). This option
        can also be used when densmap=False to calculate the densities for
        UMAP embeddings.
    disconnection_distance: float(optional, default np.inf or maximal value for bounded distances)
        Disconnect any vertices of distance greater than or equal to disconnection_distance when approximating the
        manifold via our k-nn graph. This is particularly useful in the case that you have a bounded metric.  The
        UMAP assumption that we have a connected manifold can be problematic when you have points that are maximally
        different from all the rest of your data.  The connected manifold assumption will make such points have perfect
        similarity to a random set of other points.  Too many such points will artificially connect your space.
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

"""
SNN parameters:
---------------
    neighbor_num : int
        K number of neighbors to consider for shared nearest neighbor similarity
    min_shared_neighbor_proportion : float [0, 1]
        Proportion of the K nearest neighbors that need to share two data points to be considered part of the same cluster
    Note: Naming conventions for attributes are based on the analogous ones of DBSCAN
"""
SNN_PARAMS = {
    "neighbor_num": 20,
    "min_shared_neighbor_proportion": 1 / 15,
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
MNN_PARAMS = {
    "k1": 5,
    "k2": 50,
}

# MNN_PARAMS = {
#  "var_adj": True,
#  "svd_mode": "rsvd",
#  "svd_dim": None,
# }
# @see https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.louvain.html
LOUVAIN_PARAMS = {
    "resolution": 5,
}

# @see https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html
SC_NEIGH_PARAMS = {
    "n_neighbors": 20,
    "n_pcs": 70,
}
