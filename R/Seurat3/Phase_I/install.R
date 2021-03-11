conda create --name r3.6.2
conda activate r3.6.2
conda install -c conda-forge r-base rpy2 python matplotlib numpy pandas scipy anndata r-dplyr r-qlcmatrix r-tidyr r-matrix xlrd notebook openssl pyopenssl python-igraph r-cellranger r-data.table r-ggpubr hdf5 ipykernel ipython jpeg jupyterlab r-tidyverse r-matrix r-reshape2 r-devtools leidenalg r-leiden r-kableextra r-hdf5r r-svmisc r-patchwork r-remotes r-heatmaply r-slam r-eulerr r-openxlsx r-shinydashboard r-shinydashboardplus r-biocmanager r-hmisc r-networkd3
conda install -c bioconda r-seurat scanpy bbknn anndata2ri bioconductor-dropletutils bioconductor-scater bioconductor-singler bioconductor-singlecellexperiment bioconductor-mast bioconductor-monocle
conda install -c r rstudio
#pip install harmonyTS
#remotes::install_github("mojaveazure/seurat-disk")
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
