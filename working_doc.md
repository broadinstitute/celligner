## what mfmap is not doing:

- no purity estimation?
- not robust to outliers / not able to detect undif cluster
- Maybe it should use gene network from cmap data.

## what we would want to do

min(A - (X_a*Y_a + I_a))
min(B - (X_b*Y_b + I_b)) s.t. min(MIN(dist(X_a, X_b))) ; max(Y_a\*Y_b)

## ways to test:

- distance between known good similar lines.
  - list of close matching lines using other paper's data and Allie's data (tumorComparer)
  - fake tumor data (cell line + purity*normal)
- ability to find out misslabeled lines:
  - use list of putative mislabelled (outliers in the bioarxiv paper)
  - create fake misslabeling
- using HCMI's line
- does known gene dependency of typical cancer lines match with clustering?
- ask sanger for their RNAseq

- celligner that uses gene loadings found by Josh's tool (Webster)

- celligner is already working on subsets:
  - set of genes are droped when cPCA, when mNN realignment

- celligner that finds [a mapping / a cell line for a tumor] given a specific gene grouping (dependency/geneset..)
- celligner that finds the best set of lines for all groupings for a cancer specific expression
