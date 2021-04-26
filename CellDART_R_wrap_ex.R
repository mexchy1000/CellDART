## CellDART in R
library(reticulate)
use_condaenv(condaenv = 'CellDART', required = TRUE)
py_config()

# Appropriate conda environment to run CellDART
# 1. Install conda environment
# conda create -n CellDART python=3.7 scanpy=1.5.1 scanorama=1.6 numba=0.52.0 umap-learn=0.4.6 pandas jupyter numpy seaborn leidenalg ipython python-igraph louvain matplotlib keras=2.3.1 tensorflow=1.14.0 h5py=2.10.0 --channel bioconda --channel conda-forge
# 
# 2. Install keras-gpu
# conda install tensorflow-gpu
# 
# 3. Install bbknn
# conda install -c bioconda bbknn

# CellDART_R_wrap.py should be in the same directory with da_cellfraction.py and utils.py
# set directory for CellDART_R.py, da_cellfraction.py and utils.py
setwd('/home/user/DATA1/Spatial/')

## RUN CellDART
# Command for implementing CellDART
path_os = paste0("python ./CellDART/CellDART_R_wrap.py",
                " --spdir ./Mouse_sp/ --spfilter T",
                " --spfilgene 5 --spfilcell 50",
                " --scdir ./Mouse_sc/ --sc10x_mtx T --sc10x_h5 F --sctranspose F",
                " --celltype cluster --num_markers 20 --nmix 10",
                " --npseudo 20000 --alpha 1 --alpha_lr 5 --batch_size 512",
                " --emb_dim 64 --n_iterations 3000 --init_train_epoch 10",
                " --outdir ./CellDART_output/")
# RUN
system(path_os)

## Explanation of the variables
# spdir: file directory for spatial data 
# -> Visium data should be separated in different folders

# spfilter: check whether to filter the number of cells and genes in spatial data (T: run filter)
# spfilgene: keep genes that are expressed in at least 'spfilgene' number of cells (default = 5)
# spfilcell: keep cells with at least 'spfilcell' counts (default = 50)

# scdir: file directory for single cell data 
# -> each single-cell data should be separated in different folders 
# -> each file formats should be among 10x format or others (.mtx.gz, .h5ad, h5, .csv, .tsv, or .txt)
# -> and metadata with corresponding barcode name as index should be included in metadata folder of each single-cell data
# -> metadata should be csv format
## Example directory (reference)
# -> one single cell dataset (10x mtx format) with metadata
# ./Mouse_sc/first/barcodes.tsv, ./Mouse_sc/first/genes.tsv, ./Mouse_sc/first/matrix.mtx, ./Mouse_sc/first/metadata/metadata.csv

# sc10x_mtx: check whether single-cell data is 10x genomics formatted mtx directory (T) or not (F)
# sc10x_h5: check whether single-cell data is 10x genomics formatted hdf5 file (T) or not (F)
# sctranspose: if sc10x_mtx and sc10x_h5 is F, check whether loaded matrix should be transposed (T) or not (F)

# celltype: the metadata column name for the cell type clusters (default = cell_type)
# num_markers: number of selected marker genes in each cell-type (default = 20)
# nmix: sampling number of cells in pseudospot (default = 10)
# npseudo: a total number of pseudospots (default = 20000)

# alpha: loss weights of domain classifier to source classifier
# alpha_lr: learning rate for domain classifier (alpha_lr*0.001, default = 5)
# batchsize: minibatch size during the training (default = 512)
# emb_dim: output size of dimensions for feature extractor (default = 64)

# n_iterations: iteration number for the adversarial learning (default = 3000)
# init_train_epoch: iteration number of pre-train (default = 10)

# outdir: directory to save the output file (predicted cell fraction & trained models)
