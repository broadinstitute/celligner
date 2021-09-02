## celligner
FROM debian:latest
MAINTAINER Jeremie Kalfon

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils && \
  apt-get install -y sudo &&\
  sudo apt-get install -y wget libterm-readline-gnu-perl &&\

  # all nice packages
  ## install the [tools](https://www.datacamp.com/community/tutorials/google-cloud-data-science) sudo apt-get -y install htop parallel curl  tar  vim  nano  bzip2  unzip libssl-dev  make cmake libcurl4-openssl-dev  default-jre  && sudo apt-get -y install dirmngr apt-transport-https  ca-certificates  gnupg2  software-properties-common  zlib1g-dev  libbz2-dev  liblzma-dev  openssh-server  default-libmysqlclient-dev  acl  g++
  ## sudo apt install git libmagickwand-dev libtool libexpat1-dev ghostscript graphviz libgraphviz-dev pkg-config libxml-simple-perl zlib1g-dev
  sudo apt-get -y install \
  htop \
  parallel \
  curl \
  tar \
  vim \
  nano \
  bc \
  bzip2 \
  unzip \
  libssl-dev \
  make \
  cmake \
  libcurl4-openssl-dev \
  default-jre \
  && sudo apt-get -y install dirmngr apt-transport-https \
  ca-certificates \
  gnupg2 \
  software-properties-common \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libxml2-dev \
  openssh-server \
  default-libmysqlclient-dev \
  acl \
  g++ \
  autoconf \
  automake \
  git\
  libmagickwand-dev\
  libtool \
  libexpat1-dev \
  ghostscript \
  graphviz \
  libgraphviz-dev \
  pkg-config \
  libxml-simple-perl \
  zlib1g-dev &&\

# install R	sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && echo "deb http://http.debian.net/debian sid main" | sudo tee -a /etc/apt/sources.list && echo "deb http://ftp.de.debian.org/debian testing main" | sudo tee -a /etc/apt/sources.list && sudo add-apt-repository 'deb http://cran.rstudio.com/bin/linux/debian buster-cran35/' && sudo apt update && sudo apt -y install r-base && sudo apt -y install python3 python3-pip &&
  sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
  echo "deb http://http.debian.net/debian sid main" | sudo tee -a /etc/apt/sources.list && \
  echo "deb http://ftp.de.debian.org/debian testing main" | sudo tee -a /etc/apt/sources.list &&\
  sudo add-apt-repository 'deb http://cran.rstudio.com/bin/linux/debian buster-cran35/' && \
  sudo apt update && \
  sudo apt -y install r-base &&\
  sudo apt -y install python3 python3-pip &&\
  # all python config pip3 install numpy pandas && pip3 install MACS2 && pip3 install dxpy  jupytext  scikit-learn  google-api-core  igv  igv-jupyter firecloud-dalmatian  awscli  seaborn  pipreqs  pysradb  nbstripout  bokeh  matplotlib  deeptools  tensorflow  cutadapt ipykernel jupyter_contrib_nbextensions && jupyter contrib nbextension install && nbstripout --install --global && ipykernel install && nbstripout --install --global
  pip3 install numpy pandas &&\

# search history
  touch ~/.inputrc &&\
  echo "$include /etc/inputrc" > ~/.inputrc &&\
  echo ""\e[A":history-search-backward" > ~/.inputrc &&\
  echo ""\e[B":history-search-forward" > ~/.inputrc &&\

# all R config
  export R_LIBS="~/R/x86_64-pc-linux-gnu-library/3.6" &&\
# install all nice R packages R -e "install.packages(c('plyr','dplyr','tidyverse','magrittr','reshape2','useful','ggplot2','ggthemes','ggrepel','gridExtra','ggridges','GGally','plotly','VennDiagram','RColorBrewer','extrafont','cowplot', 'network','data.table','DT','readr','readxl','clues','mclust','pheatmap','Rtsne','NMF','hash', 'stringr', 'irr', 'zoo', 'devtools', 'scales', 'rlang', 'rmarkdown','lsa','BiocManager'), dependencies=TRUE, repos='http://cran.rstudio.com/');  font_import(); loadfonts(); BiocManager::install(c('GSEABase','limma','org.Hs.eg.db','GenomicRanges','DESeq2'));	  print('if can't use broad intranet, install from source with [R CMD INSTALL .] for 'taigr',	    'cdsr','svacd', 'cell_line_mapping/celllinemapr')"
  R -e "install.packages(c('plyr','dplyr','tidyverse','magrittr','reshape2','useful','ggplot2','ggthemes',\
  'ggrepel','gridExtra','ggridges','GGally','plotly','VennDiagram','RColorBrewer','extrafont','cowplot',\
  'network','data.table','DT','readr','readxl','clues','mclust','pheatmap','Rtsne','NMF','hash',\
  'stringr', 'irr', 'zoo', 'devtools', 'scales', 'rlang', 'rmarkdown','lsa','BiocManager'), \
  dependencies=TRUE, repos='http://cran.rstudio.com/'); \
  font_import();\
  loadfonts();\
  BiocManager::install(c('GSEABase','Seurat','batchelor','limma','org.Hs.eg.db','GenomicRanges','DESeq2'))\
  devtools::install_github('broadinstitute/taigr');\
  devtools::install_github('broadinstitute/celligner');\
  " &&\
  R -e "library(devtools);\
        devtools::install_github('broadinstitute/celligner')"
