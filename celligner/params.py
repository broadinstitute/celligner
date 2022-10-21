# Oncotree tissue colors
TISSUE_COLOR_OT = {
    "Adrenal Gland": "#E13978",
    "Ampulla of Vater": "#F5899E",
    "Biliary Tract": "#C091E3",
    "Bladder/Urinary Tract":"#E08571",
    "Bone": "#9F55BB",
    "Breast":"#45A132",
    "Bowel":"#96568E",
    "CNS/Brain": "#F5899E",
    "Cervix":"#5AB172",
    "Esophagus/Stomach": "#DFBC3A",
    "Eye": "#349077",
    "Fibroblast": "#D8AB6A",
    "Embryonal":"#75DFBB",
    "Head and Neck": "#5DA134",
    "Kidney": "#1F8FFF",
    "Liver": "#9C5E2B",
    "Lung": "#51D5E0",
    "Lymphoid": "#ABD23F",
    "Myeloid": "#DA45BB",
    "Normal":"#555555",
    "Ovary/Fallopian Tube": "#56E79D",
    "Pancreas": "#B644DC",
    "Peripheral Nervous System": "#73E03D",
    "Pleura": "#F5899E", ###
    "Prostate": "#3870C9",
    "Skin": "#6C55E2",
    "Soft Tissue": "#5FDB69",
    "Testis": "#F5899E", ###
    "Thymus": "#659FD9", 
    "Thyroid": "#D74829",
    "Other/Unknown": "#bdbdbd",
    "Uterus": "#E491C1",
    "Vulva/Vagina":"#E491C1"
}

TISSUE_COLOR = {
    "engineered": "#bcdfbd",
    "fibroblast": "#9eAeAe",
    "other": "#A3969d",
    "skin": "#969696",
    "soft_tissue": "#cedb9c",
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
    "biliary_tract": "#e7ba52",
    "adrenal": "#8ca252",
    "thymus": "#659fd9"
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



#mnn_ndist = 3, # ndist parameter used for MNN

# Differentially expressed genes with a rank better than this is in the cell line
# or tumor data are used to identify mutual nearest neighbors in the MNN alignment step
TOP_K_GENES = 1000

# number of PCs to use for dimensionality reduction
PCA_NCOMP = 70 

# number of cPCA dimensions to regress out of the data
CPCA_NCOMP = 4

# @see https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.louvain.html
LOUVAIN_PARAMS = {
    "resolution": 5, # resolution parameter used for clustering the data
}

# For Mariona method (default)
MNN_PARAMS = {
    "k1": 5, # number of nearest neighbors of tumors in the cell line data
    "k2": 50, # number of nearest neighbors of cell lines in the tumor data
    "cosine_norm": False,
    "fk": 5 
}

UMAP_PARAMS = {
    "n_neighbors": 10, # num nearest neighbors used to create UMAP plot
    "n_components": 2, 
    "metric": "euclidean", # distance metric used for the UMAP projection
    "min_dist": 0.5 # min distance used to create UMAP plot
}
