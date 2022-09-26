import celligner
import pandas as pd
import numpy as np
import re

from taigapy import TaigaClient
tc = TaigaClient()

##### 1. LOAD DATA - expression and annotations for the same samples #####
# CCLE
CCLE_expression = tc.get(name='internal-21q3-fe4c', version=16, file='CCLE_expression_full')
CCLE_annot = tc.get(name='internal-21q3-fe4c', version=16, file='sample_info')

# TCGA
TCGA_expression = tc.get(name='celligner-input-9827', version=1, file='tumor_expression')
TCGA_annot = tc.get(name='celligner-input-9827', version=1, file='tumor_annotations')

# # met500 
# met500_ann = tc.get(name='met500-fc3c', version=1, file='met500_ann')
# met500_meta = tc.get(name='met500-fc3c', version=1, file='met500_meta')
# met500_TPM = tc.get(name='met500-fc3c', version=1, file='met500_TPM') #20,979x868 matrix

# # Novartis_PDX
# Novartis_PDX_ann = tc.get(name='pdx-data-3d29', version=2, file='Novartis_PDX_ann')
# Novartis_PDX_TPM = tc.get(name='pdx-data-3d29', version=2, file='Novartis_PDX_TPM').T # 38,087x445

# # Pediatric_PDX
# pediatric_PDX_ann = tc.get(name='pdx-data-3d29', version=2, file='pediatric_PDX_ann')
# pediatric_PDX_TPM = tc.get(name='pdx-data-3d29', version=2, file='pediatric_PDX_TPM') #80,000x250


###### 2. IDENTIFY COMMON GENES #####
hgnc_complete_set = tc.get(name='hgnc-87ab', version=5, file='hgnc_complete_set')
func_genes = hgnc_complete_set[~hgnc_complete_set.locus_group.isin(["non-coding RNA", "pseudogene"])]

# Filter to columns with ensembl id
CCLE_expression = CCLE_expression.filter(like='ENSG')
CCLE_expression.columns = pd.Series(CCLE_expression.columns).apply(lambda x: re.search('(ENSG\d+)', x).group(1))

# Subset to common genes
common_genes = set(CCLE_expression.columns).intersection(set(TCGA_expression.columns))
# common_genes = common_genes.intersection(set(met500_TPM.columns))
# common_genes = common_genes.intersection(set(Novartis_PDX_TPM.columns))
# common_genes = common_genes.intersection(set(pediatric_PDX_TPM.columns))
common_genes = list(common_genes.intersection(set(func_genes.ensembl_gene_id)))
print('Common genes:', len(common_genes))

CCLE_expression = CCLE_expression[common_genes]
TCGA_expression = TCGA_expression[common_genes]
# met500_TPM = met500_TPM[common_genes]
# Novartis_PDX_TPM = Novartis_PDX_TPM[common_genes]
# pediatric_PDX_TPM = pediatric_PDX_TPM[common_genes]


##### 3. CLEAN UP ANNOTATIONS #####
# Specify tissue_type, disease_type and cell_type columns in both annotation dataframes

# Some tissue names are not consistent between the datasets
renamed_tissues = {np.nan: "unknown", "adrenal_cortex": "adrenal", "colorectal": "colon", 
                   'thymus': 'thyroid', 'meninges':"central_nervous_system", None: "unknown", 
                   'brain': "central_nervous_system"}

# CCLE
CCLE_annot = CCLE_annot.set_index("DepMap_ID")
CCLE_annot = CCLE_annot.loc[CCLE_expression.index, ["lineage", 'Subtype']]\
                       .rename(columns={"lineage": "tissue_type", "Subtype": 'disease_type'})
CCLE_annot["cell_type"] = "CCLE cell line"
CCLE_annot = CCLE_annot.replace({"tissue_type": renamed_tissues})

# TCGA
TCGA_annot = TCGA_annot.set_index("sampleID").loc[TCGA_expression.index, ["lineage","subtype"]]\
                       .rename(columns={"lineage":"tissue_type", "subtype": 'disease_type'})
TCGA_annot['cell_type'] = "TCGA tumor"
TCGA_annot = TCGA_annot.replace({"tissue_type": renamed_tissues})

# # met500
# met500_meta["primary_site"] = met500_ann['primary_site'].values
# met500_ann = met500_meta.rename(columns={"Sample_id": 'sample_id', 'tissue': 'tissue_type',
#                                          'primary_site': "disease_type", "sample_type": "cell_type"})\
#                         .set_index('sample_id', drop=True)[["tissue_type","disease_type","cell_type"]]\
#                         .replace({"tissue_type":renamed_tissues, "cell_type": {"tumor": "met500 tumor"}})

# # Novartis_PDX
# Novartis_PDX_ann = Novartis_PDX_ann.rename(columns={"sampleID": 'sample_id', 'lineage': 'tissue_type', 
#                                                     'subtype': "disease_type", "type": "cell_type"})\
#                                    .set_index('sample_id', drop=True)[['cell_type', 'disease_type', 'tissue_type']]\
#                                    .replace({"tissue_type":renamed_tissues})
# Novartis_PDX_ann = Novartis_PDX_ann.loc[Novartis_PDX_TPM.index,:]

# # Pediatric_PDX
# pediatric_PDX_ann = pediatric_PDX_ann.rename(columns={"sampleID": 'sample_id', 'lineage': 'tissue_type', 
#                                                       'subtype': "disease_type", "type": "cell_type"})\
#                                     .set_index('sample_id', drop=True)[['cell_type', 'disease_type', 'tissue_type']]\
#                                     .replace({"tissue_type":renamed_tissues})
# pediatric_PDX_ann = pediatric_PDX_ann.loc[pediatric_PDX_TPM.index]



##### RUN CELLIGNER #####
# Create Celligner object and fit + transform the expression datasets
my_celligner = celligner.Celligner()
my_celligner.fit(CCLE_expression, CCLE_annot)
my_celligner.transform(TCGA_expression, TCGA_annot)
my_celligner.computeMetricsForOutput()
my_celligner.save("model.pkl")

# Multi-dataset alignment - sequentially aligning additional expression datasets
# my_celligner.makeNewReference()
# my_celligner.mnn_kwargs.update({"k1":20, "k2":50})
# my_celligner.transform(met500_TPM, met500_ann, compute_cPCs=False)

# my_celligner.makeNewReference()
# my_celligner.mnn_kwargs.update({"k1":10, "k2":50})
# my_celligner.transform(Novartis_PDX_TPM, Novartis_PDX_ann, compute_cPCs=False)

# my_celligner.makeNewReference()
# my_celligner.mnn_kwargs.update({"k1":10, "k2":50})
# my_celligner.transform(pediatric_PDX_TPM, pediatric_PDX_ann, compute_cPCs=False)

# my_celligner.computeMetricsForOutput()
# my_celligner.save("model_multi.pkl")
