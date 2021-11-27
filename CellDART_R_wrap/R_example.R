### 0. Install required packages
# Install Seurat
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")

# Install dplyr
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

## Install sceasy
# https://github.com/cellgeni/sceasy
# 1. Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2. Install LoomExperiment
if (!requireNamespace("LoomExperiment", quietly = TRUE))
  BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))

# 3. Install sceasy
if (!requireNamespace("sceasy", quietly = TRUE))
  devtools::install_github("cellgeni/sceasy")



### 1. Example of using function "pred_cellf_celldart"
library(Seurat)

# Find the directory for active script file (file_path)
file_path <- rstudioapi::getSourceEditorContext()$path
file_path <- strsplit(file_path, split='/')
file_path <- paste(file_path[[1]][-length(file_path[[1]])],
                    collapse=.Platform$file.sep)

# Set working directory
setwd(file_path)

# Import the CellDART function
source('cellf_pred.R')

# Make output folder
output_folder_name <- 'CellDART_output'
if (!file.exists(output_folder_name)){
  dir.create(output_folder_name)
}


## Load single-cell and spatial datasets
# Load single-cell dataset (RDS file with Seurat object): GSE115746 (mouse single cell: ALS and VISp)
sc_data <- readRDS('sc_data.rds')

# Load spatial dataset (RDS file with Seurat object): 10X genomics data repository
# V1_Mouse_Brain_Sagittal_Anterior
# V1_Mouse_Brain_Sagittal_Posterior
sp_data <- readRDS('sp_data.rds')


## Check the size of spatial dataset
dim(sp_data)

## Set the number of pseudospots: 5 times the number of spatial spots
npseudo <- 5*dim(sp_data)[2]



## Perform CellDART analysis
## Explanation of the variables
# celldart.dir: R wrapper python file for the celldart analysis
# outdir: file directory to save the model and result
# sp_data: spatial data (Seurat object) to be used in predicting cell fraction (default: NULL)
# -> If NULL, then use the spatial data already saved in spdir
# sc_data: single-cell data (Seurat object) to be used in making pseudospots (default: NULL)
# -> If NULL, then use the spatial data already saved in scdir

# sp_subset: whether to subset spatial data and calculate for specific spot cluster (default = FALSE)
# spot.cluster.name: group name of the cluster used for subsetting spatial data (default: 'seurat_clusters')
# spot.cluster.of.interest: name of each spot clusters to be used (default: NULL)
# metadata_celltype: column name for single-cell annotation data in metadata (default: 'celltype')

# virtual.env.name: name of the conda environment to use for CellDART analysis (default: 'CellDART')

# gpu: check whether to use gpu (TRUE) or not (FALSE) (default = TRUE)

# sp10x: whether the spatial data is 10x Visium format (TRUE) or not (FALSE) (default: TRUE)
# spfilter: check whether to filter the number of cells and genes in spatial data (TRUE: run filter)
# spfilgene: keep genes that are expressed in at least 'spfilgene' number of cells (default = 5)
# spfilspot: keep spots with at least 'spfilcell' counts (default = 50)

# sc10x_mtx: check whether single-cell data is 10x genomics formatted mtx directory (TRUE) or not (FALSE)
# sc10x_h5: check whether single-cell data is 10x genomics formatted hdf5 file (TRYE) or not (FALSE)
# sctranspose: if sc10x_mtx and sc10x_h5 is F, check whether loaded matrix should be transposed (TRUE) or not (FALSE)

# num_markers: number of selected marker genes in each cell-type (default = 20)

# seed_num: seed to be used in random sampling (default = 0)

# nmix: sampling number of cells in pseudospot (default = 8)
# npseudo: a total number of pseudospots (default = 20000)

# alpha: loss weights of domain classifier to source classifier
# alpha_lr: learning rate for domain classifier (alpha_lr*0.001, default = 5)
# batch_size: minibatch size during the training (default = 512)
# emb_dim: output size of dimensions for feature extractor (default = 64)

# n_iterations: iteration number for the adversarial learning (default = 3000)
# init_train_epoch: iteration number of pre-train (default = 10)

sp_data_cellf <- pred_cellf_celldart(celldart.dir='CellDART_R_wrap.py',
                                     outdir=file.path(output_folder_name),
                                     sp_data=sp_data, sc_data=sc_data,
                                     spdir=NULL,scdir=NULL,
                                     sp_subset=FALSE,spot.cluster.name='seurat_clusters',
                                     spot.cluster.of.interest=NULL,
                                     metadata_celltype='cell_subclass',
                                     virtual.env.name='spatial',gpu=TRUE,
                                     sp10x=FALSE,spfilter=FALSE,spfilgene=0,spfilspot=0,
                                     sc10x_mtx=FALSE,sc10x_h5=FALSE,sctranspose=FALSE,
                                     seed_num=0, 
                                     num_markers=20, nmix=8, npseudo=npseudo,
                                     alpha=0.6, alpha_lr=5,
                                     emb_dim=64,batch_size=512,
                                     n_iterations=3000, init_train_epoch=10)

## Save seurat object with cell fraction
saveRDS(sp_data_cellf, file.path('CellDART_output', 'sp_data_cellf.rds'))



### 2. Visualization of spatial cell fraction
# Remove '_cellf' from the column names cell fraction metadata
cellf.data <- sp_data_cellf@meta.data
cellf.data.colname <- sapply(colnames(cellf.data), function(x){
  if (grepl('_cellf',x)){return(strsplit(x, split='_cellf')[[1]][1])}
  else {return(x)}
})
sp_data_cellf.mod <- sp_data_cellf
colnames(sp_data_cellf.mod@meta.data) <- cellf.data.colname

# Visualize the layer-specific excitatory neuons
cell_types <- c("L2.3.IT","L4","L5.IT","L5.PT","L6b","L6.CT","L6.IT")
p <- SpatialFeaturePlot(sp_data_cellf.mod, features = cell_types, 
                        ncol = 4, alpha=0.6, combine = FALSE)

patchwork::wrap_plots(p, ncol = 8)
