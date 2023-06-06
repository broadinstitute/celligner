import celligner
import pandas as pd
import re
from taigapy import TaigaClient

tc = TaigaClient()

portal_in_dict = {}

portal_out_dict = {}

depmap_params = {'name': 'depmap_public_23q2',
                 'taiga_name': 'public-23q2-19de',
                 'taiga_file': 'OmicsExpressionProteinCodingGenesTPMLogp1',
                 'dset_type': 'model',
                 'mnn_params': None}

tcga_params = {'name': 'tcga',
               'taiga_name': 'celligner-input-9827',
               'taiga_file': 'tumor_expression',
               'dset_type': 'tumor',
               'mnn_params': None}

met500_params = {'name': 'met500',
                 'taiga_name': 'met500-fc3c',
                 'taiga_file': 'met500_TPM',
                 'dset_type': 'tumor',
                 'mnn_params': {'k1': 20, 'k2': 50}}

pdx_nv_params = {'name': 'pdx_novartis',
                 'taiga_name': 'pdx-data-3d29',
                 'taiga_file': 'Novartis_PDX_TPM',
                 'dset_type': 'model',
                 'mnn_params': {'k1': 10, 'k2': 50}}

pdx_ped_params = {'name': 'pdx_pediatric',
                  'taiga_name': 'pdx-data-3d29',
                  'taiga_file': 'pediatric_PDX_TPM',
                  'dset_type': 'model',
                  'mnn_params': {'k1': 10, 'k2': 50}}

celligner_default_extras = [met500_params, pdx_nv_params, pdx_ped_params]


def process_data(bg_df, contrast_df, extra_dfs):
    # Filter to columns with ensembl id
    # Need to transform from ensembl transcript ids to ensemble gene ids
    # insert below

    # Load HGNC gene set, filter to functional subset
    hgnc_complete_set = tc.get(name='hgnc-87ab', version=7, file='hgnc_complete_set')
    func_genes = hgnc_complete_set[~hgnc_complete_set.locus_group.isin(["non-coding RNA", "pseudogene"])]

    # background datafrrame should have genes labeled by gene symbol
    # need to reindex
    # print('bg_genes going in:', bg_df.keys())
    bg_genes = pd.Series(bg_df.keys()).apply(lambda s: re.search(r'^([\w.-]+) \(', s).group(1)).rename('symbol')
    bg_genes = bg_genes.set_axis(bg_df.keys()).to_frame()
    bg_genes = bg_genes.reset_index().merge(hgnc_complete_set[['symbol', 'ensembl_gene_id']],
                              left_on='symbol', right_on='symbol')[['index','ensembl_gene_id']].set_index('index')
    # print('bg_genes processed:', bg_genes)
    bg_df = bg_df.rename(columns=bg_genes.to_dict()['ensembl_gene_id'])
    # print(bg_df.head())

    gene_sets = [set(bg_df.columns), set(contrast_df.columns)] + [set(_df.columns) for _df in extra_dfs]
    gene_set = set(func_genes.ensembl_gene_id).intersection(*gene_sets)
    common_genes = [x for x in bg_df.columns if x in gene_set]
    print('Common genes:', len(common_genes))

    bg_df = bg_df[common_genes]
    contrast_df = contrast_df[common_genes]
    extra_dfs = [_df[common_genes] for _df in extra_dfs]

    return bg_df, contrast_df, extra_dfs


def run_celligner(bg=depmap_params, contrast=tcga_params, extra_dsets=celligner_default_extras):
    # load data to be used as the background and label the source
    bg_df = tc.get(name=bg['taiga_name'], file=bg['taiga_file'])
    bg_df = pd.concat({bg['name']: bg_df}, names=['source'])

    # load data to be contrasted with the background and label the source
    contrast_df = tc.get(name=contrast['taiga_name'], file=contrast['taiga_file'])
    contrast_df = pd.concat({contrast['name']: contrast_df}, names=['source'])

    # if there are additional datasets to be projected then collect and label them as well
    extra_dfs = []
    for dset in extra_dsets:
        _df = tc.get(name=dset['taiga_name'], file=dset['taiga_file'])
        if dset['name'] == 'pdx_novartis':
            _df = _df.T
        _df = pd.concat({dset['name']: _df}, names=['source'])
        extra_dfs.append(_df)

    '''print(bg_df.head())
    print(contrast_df.head())
    for _df in extra_dfs:
        print(_df.head())'''
    # make sure all datasets are using ensembl gene ids and are restricted to common sets of genes
    bg_df, contrast_df, extra_dfs = process_data(bg_df, contrast_df, extra_dfs)
    '''print(bg_df.head())
    print(contrast_df.head())
    for _df in extra_dfs:
        print(_df.head())'''
    # Create Celligner object and fit + transform the reference (depmap) and target (TCGA) expression datasets
    my_celligner = celligner.Celligner()
    my_celligner.fit(bg_df.droplevel(0,0))
    my_celligner.transform(contrast_df.droplevel(0,0))
    # add in additional datasets to be projected if they are given
    for _df in extra_dfs:
        df_name = _df.index.get_level_values(0).unique()[0]
        for dset in extra_dsets:
            if dset['name'] == df_name:
                break

        my_celligner.makeNewReference()
        if dset['mnn_params']:
            p = dset['mnn_params']
            my_celligner.mnn_kwargs.update({'k1': p['k1'], 'k2': p['k2']})
        my_celligner.transform(_df.droplevel(0,0), compute_cPCs=False)

    # Compute UMAP, clusters and tumor - model distance
    model_ids = list(bg_df.index.get_level_values(1))
    tumor_ids = list(contrast_df.index.get_level_values(1))

    for _df in extra_dfs:
        df_name = _df.index.get_level_values(0).unique()[0]
        for dset in extra_dsets:
            if dset['name'] == df_name:
                break

        if dset['dset_type'] == 'model':
            model_ids += list(_df.index.get_level_values(1))
        elif dset['dset_type'] == 'tumor':
            tumor_ids += list(_df.index.get_level_values(1))

    my_celligner.computeMetricsForOutput(model_ids=model_ids, tumor_ids=tumor_ids)

    outname = 'celligner_output_' \
              + bg['name'] + '_' \
              + contrast['name'] + '_' \
              + '_'.join([d['name'] for d in extra_dsets]) + ".pkl"
    my_celligner.save(outname)
    print('Model saved to: ', outname)
    return outname


if __name__ == "__main__":
    run_celligner()
