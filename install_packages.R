options(repos = c("https://cran.cnr.berkeley.edu"))

cran_packages <- c('tidyverse', 'reshape2', 'plyr', 'data.table', 'Seurat', 'pdist','FNN', 'irlba')
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

bioconductor_packages <- c('limma', 'batchelor', 'BiocParallel')
new_bioconductor_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
if(length(new_bioconductor_packages)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(new_bioconductor_packages)
}

