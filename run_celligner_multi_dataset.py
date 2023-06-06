import celligner
import pandas as pd
import numpy as np
import re
from taigapy import TaigaClient
tc = TaigaClient()

# Load data
CCLE_expression = tc.get(name='dmc-22q2-5e51', version=16, file='CCLE_expression_full')
tumor_expression = tc.get(name='celligner-input-9827', version=2, file='tumor_expression')
met500_TPM = tc.get(name='met500-fc3c', version=1, file='met500_TPM')
Novartis_PDX_TPM = tc.get(name='pdx-data-3d29', version=2, file='Novartis_PDX_TPM').T
pediatric_PDX_TPM = tc.get(name='pdx-data-3d29', version=2, file='pediatric_PDX_TPM')

# Filter to columns with ensembl id
CCLE_expression = CCLE_expression.filter(like='ENSG')
CCLE_expression.columns = pd.Series(CCLE_expression.columns).apply(lambda x: re.search('(ENSG\d+)', x).group(1))

## Load HGNC gene set, filter to functional subset
hgnc_complete_set = tc.get(name='hgnc-87ab', version=5, file='hgnc_complete_set')
func_genes = hgnc_complete_set[~hgnc_complete_set.locus_group.isin(["non-coding RNA", "pseudogene"])]

# Identify common genes - maintaining order from CCLE expression matrix
gene_sets = [set(tumor_expression.columns), set(met500_TPM.columns), set(Novartis_PDX_TPM.columns), set(pediatric_PDX_TPM.columns)]
gene_set = set(func_genes.ensembl_gene_id).intersection(*gene_sets)
common_genes = [x for x in CCLE_expression.columns if x in gene_set]
print('Common genes:', len(common_genes))

# Subset all matrices to common genes
CCLE_expression = CCLE_expression[common_genes]
tumor_expression = tumor_expression[common_genes]
met500_TPM = met500_TPM[common_genes]
Novartis_PDX_TPM = Novartis_PDX_TPM[common_genes]
pediatric_PDX_TPM = pediatric_PDX_TPM[common_genes]

# Create Celligner object and fit + transform the reference (CCLE) and target (TCGA) expression datasets
my_celligner = celligner.Celligner()
my_celligner.fit(CCLE_expression)
my_celligner.transform(tumor_expression)

# Multi-dataset alignment - sequentially aligning additional expression datasets
# Met500
my_celligner.makeNewReference()
my_celligner.mnn_kwargs.update({"k1":20, "k2":50})
my_celligner.transform(met500_TPM, compute_cPCs=False)
# Novartis PDX
my_celligner.makeNewReference()
my_celligner.mnn_kwargs.update({"k1":10, "k2":50})
my_celligner.transform(Novartis_PDX_TPM, compute_cPCs=False)
# Pediatric PDX
my_celligner.makeNewReference()
my_celligner.mnn_kwargs.update({"k1":10, "k2":50})
my_celligner.transform(pediatric_PDX_TPM, compute_cPCs=False)

# Compute UMAP, clusters and tumor - model distance
model_ids = list(CCLE_expression.index)+list(Novartis_PDX_TPM.index)+list(pediatric_PDX_TPM.index)
tumor_ids = list(tumor_expression.index)+list(met500_TPM.index)
my_celligner.computeMetricsForOutput(model_ids=model_ids, tumor_ids=tumor_ids)

my_celligner.save("model_22q2_dmc_multi_dataset.pkl")
