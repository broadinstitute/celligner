import celligner
import pandas as pd
import numpy as np
import re
from taigapy import TaigaClient
tc = TaigaClient()

# Load data
CCLE_expression_full = tc.get(name='public-22q2-de04', version=14, file='CCLE_expression_full')
tumor_expression = tc.get(name='celligner-input-9827', version=2, file='tumor_expression')

# Filter to columns with ensembl id
CCLE_expression = CCLE_expression_full.filter(like='ENSG')
CCLE_expression.columns = pd.Series(CCLE_expression.columns).apply(lambda x: re.search('(ENSG\d+)', x).group(1))

## Load HGNC gene set, filter to functional subset
hgnc_complete_set = tc.get(name='hgnc-87ab', version=5, file='hgnc_complete_set')
func_genes = hgnc_complete_set[~hgnc_complete_set.locus_group.isin(["non-coding RNA", "pseudogene"])]

# Identify common genes - maintaining order from CCLE expression matrix
gene_set = set(func_genes.ensembl_gene_id).intersection(set(tumor_expression.columns))
common_genes = [x for x in CCLE_expression.columns if x in gene_set]
print('Common genes:', len(common_genes))

# Subset all matrices to common genes
CCLE_expression = CCLE_expression[common_genes]
tumor_expression = tumor_expression[common_genes]

# Create Celligner object and fit + transform the reference (CCLE) and target (TCGA) expression datasets
my_celligner = celligner.Celligner()
my_celligner.fit(CCLE_expression)
my_celligner.transform(tumor_expression)

# Compute UMAP, clusters and tumor - model distance
my_celligner.computeMetricsForOutput()

my_celligner.save("model_22q2_public.pkl")
