# CellDART: Cell type inference by domain adaptation of single-cell and spatial transcriptomic data
CellDART is a tool to estimate cell fraction of spatial transcriptomic spots using domain adaptation of deep neural networks.

![figure1png](https://user-images.githubusercontent.com/14209383/114880774-528b8100-9e3d-11eb-9b60-41c9d0acd5fd.png)


## Optimal parameter choices (for brain)
  Number of total marker genes = 200 ~ 400 (or number of markers per cluster: 10 ~ 20)  
  Number of pseudospots = 5 to 10 times the number of real spots (20,000~40,000 per Visium slide)  
  Number of sampled cells in a pseudospot (virtual mixture of single-cell data) = 8  
  Iteration number = 3,000  
  Mini-batch size = 512  
  Loss weights between source and domain classifier (alpha) = 0.6  
  Learning rate = 0.001 * alpha_lr = 0.005  

## Code Example  
python: CellDART_example_mousebrain_markers.ipynb  
R wrap: Please refer to the '/vignettes/introduction.Rmd' file  

## Python function for CellDART (pred_cellf_celldart)  
### Install conda environment and add jupyter kernel  
```Plain Text
  conda create -n CellDART python=3.7.12 -c conda-forge  
  conda activate CellDART  
  pip install git+https://github.com/mexchy1000/CellDART.git  
  python -m ipykernel install --user --name CellDART --display-name CellDART  
```
### Dependency (python)  
```Plain Text
python 3.7  
numpy 1.21.6  
pandas 1.3.5  
tensorflow 1.14.0  
scanpy 1.5.1  
seaborn 0.11.2  
keras 2.3.1  
h5py 2.10.0 
```
### Function and parameters
```Plain Text
from CellDART.pred_cellf_celldart import pred_cellf_celldart  
adata_sp = pred_cellf_celldart(adata_sp=adata_sp, adata_sc=adata_sc, count_from_raw = False,  
        　　　　　　　　　　　   gpu=True, celltype='celltype', num_markers=20,  
        　　　　　　　　　　　　 nmix=8, npseudo=20000, alpha=0.6, alpha_lr=5, batch_size=512,  
        　　　　　　　　　　　　 emb_dim=64, n_iterations=3000, init_train_epoch=10,  
        　　　　　　　　　　　　 outdir='./CellDART_output', return_anndata=True)
```        
**(1) adata_sp:** spatial data (AnnData object) with raw count matrix to be used in predicting cell fraction (default: None)  
**(2) adata_sc:** single-cell data (AnnData object) with raw count matrix to be used in making pseudospots (default: None)  
**(3) count_from_raw:** whether to extract count matrix frow .raw of AnnData (default: False)  
-> non-normalized raw count matrix should be contained in the AnnData .raw file  
-> if False, then utilize the count matrices saved in adata_sp and adata_sc directly  
**(4) gpu:** check whether to use gpu (True) or not (False) (default = True)  
**(5) celltype:** column name for single-cell annotation data in .obs (default: 'celltype')  
**(6) num_markers:** number of selected marker genes in each celltype (default = 20)   
**(7) nmix:** sampling number of cells in pseudospot (default = 10)  
**(8) npseudo:** a total number of pseudospots (default = 20,000)  
**(9) alpha:** loss weights of the domain classifier to the source classifier (default = 0.6)  
**(10) alpha_lr:** learning rate for the domain classifier (alpha_lr*0.001, default = 5)  
**(11) batch_size:** minibatch size for pseudospots and spatial data during the training (default = 512)  
**(12) n_iterations:** iteration number for the adversarial learning (default = 3,000)  
**(13) init_train_epoch:** iteration number for the pre-training process (default = 10)  
**(14) outdir:** the directory to save output files (models and results)  
**(15) return_anndata:** whether to return spatial AnnData file with predicted cell fraction in .obs (default: False)  

## R wrap function for CellDART using reticulate  
  ```Plain Text
  devtools::install_github("mexchy1000/CellDART", build_vignettes = T, force = T)  
  library(CellDART)  
  help(pred_cellf_celldart)  # Explanation for the parameters and short examples  
  browseVignettes("CellDART")  # Browse for the vignettes (/vignettes/introduction.Rmd)
  ```
  ### Function and additional parameters
  ```Plain Text
  # Using conda environment (environment will be automatically installed in Linux distributions)
  # If using Windows, then install conda environment first and then run the function below with python.install = F
  sp_data_cellf <- pred_cellf_celldart(sp_data, sc_data, outdir = '.',
                                       sp_subset=F, spot.cluster.name='seurat_clusters',
                                       spot.cluster.of.interest=NULL,
                                       env.select='conda',python.install=T,
                                       python_path=NULL, env.name='CellDART',
                                       gpu=TRUE, metadata_celltype='celltype',
                                       num_markers=20, seed_num=0,
                                       nmix=8, npseudo=20000, alpha=0.6,alpha_lr=5,
                                       emb_dim=64,batch_size=512,n_iterations=3000, init_train_epoch=10)
  ```
  ```Plain Text
  # Using virtual environment (environment will be automatically installed in Linux distributions)
  # Not recommended for Windows
  sp_data_cellf <- pred_cellf_celldart(sp_data, sc_data, outdir = '.',
                                       sp_subset=F, spot.cluster.name='seurat_clusters',
                                       spot.cluster.of.interest=NULL,
                                       env.select='virtual',python.install=T,
                                       python_path=NULL, env.name='CellDART',
                                       gpu=TRUE, metadata_celltype='celltype',
                                       num_markers=20, seed_num=0,
                                       nmix=8, npseudo=20000, alpha=0.6,alpha_lr=5,
                                       emb_dim=64,batch_size=512,n_iterations=3000, init_train_epoch=10)
  ```
  **(1) outdir:** the directory to save output files (models and results) (default = '.')  
  **(2) sp_subset:** whether to subset spatial data and calculate for specific spot cluster (default = FALSE)  
  **(3) spot.cluster.name:** group name of the cluster used for subsetting spatial data (default = 'seurat_clusters')  
  **(4) spot.cluster.of.interest:** name of each spot clusters to be used (default = NULL)  
  **(5) env.select:** select between using reticulate virtual environment or conda environment (default = 'conda')  
  -> either of the selection will search the already installed environment  
  -> if environment is not found, then it will automatically install the new environment  
  **(6) python.install:** whether to automatically install python version 3.7.12 (default = F)  
  -> For Windows, set python.install = F  
  **(7) python_path:** path for the python 3.7.12 (default = NULL)  
  **(8) env.name:** name of the virtual or conda environment to use for the analysis (default = 'CellDART')  
  **(9) metadata_celltype:** column name for single-cell annotation data in metadata (default = 'celltype')  

  ### Dependency (R wrapper)
  ```Plain Text
  Seurat 4.0.5  
  dplyr 1.0.7  
  sceasy 0.0.6  
  reticulate 1.22  
  ```
  ### Installation in Linux distributions  
  Virtual environment (env.select="virtual") or conda environment (env.select="conda") will be automatically installed while running function 'pred_cellf_celldart'  
  Detailed explanation is in '/R/Read_R_wrap.md' file.  
  ### Installation in Windows  
  Install conda environment first and then run the function with env.select='conda' and python.install=F   
  
  
## R shiny application for CellDART (under preparation)
Shiny application for preprocessing and CellDART analysis. (inside 'shiny')  
Application is based on Seurat, sceasy, dplyr, stringr, vroom, and CellDART.  
The R shiny application panel consists of main, upload, and analysis sections.  

### 1. Upload files:  
### A. Token is proivded for each Shiny sesion.  
(1) First you should click 'Set dir with current token' to start.  
(2) If you want to return to the previous working directory, then write down the previous session token to 'Enter session token' and then click 'Find & Set dir'. You will return to the previous working directory.  
(3) Write down the name of the output folder for saving the results and click 'Create'. You can return to the previous folder by entering the according name. Also, you can create multiple folders and save different results in each folder.  
(4) When the session ends (30 min), all the results will disappear.  

### B. To upload RDS or 10X format single-cell or spatial count matrix (*.rds/*.h5). (<500 MB)  
To upload rds matrix file for example, 'sparse_matrix.rds' or 10X matrix file for example, 'filtered_feature_bc_matrix.h5' from your computer, select the corresponding data format (single-cell/spatial) and click 'Browse' button. Use this button only when the data format is rds or 10X formatted hdf5 file.  

### C. To upload tsv/csv/txt format count matrix or metadata file. (<500 MB)  
(1) Click 'Convert' button in 'Convert to matrix or metadata'. A window will appear.  
(2) Select data format (single-cell, spatial, cell metadata) and click 'Browse' to search for the file to upload.  
(3) For the cell metadata, the rownames should be barcodes names and it should correspond with single-cell data barcode names. Unless, an error will occur and CellDART may abruptly stop.  
(4) For the single-cell or spatial count matrix, rownames should be gene names and column names should be barcode names.  
(5) First to check the data composition, click 'Check' to explore the contents (~10rows) of the file uploaded.  
(6) Change the delimeter if it is not correct.  
(7) Check if the matrix should be transposed.  
(8) If the name of the first column is not empty or wrongly assigned (normally, it should be empty or non-meaningful character such as .X) and the column should be shifted.  
(9) If you finished checking, then click 'Convert' button to generate data file for CellDART analysis.  

### D. Download predicted cell fraction.  
'Download' button in 'Download predicted cell fraction' should only be used after CellDART analysis is finished and the zipped file for cell fraction data is ready to be donwloaded.  

### 2. CellDART:  
A. Select the column name of metadata representing cell types in 'Group for classifying celltypes'.  

B. Select the number of marker genes per cell type (default: 20), number of cells in a pseudospot (default: 8), and number of pseudospots (generally optimal in the range of 5~10 times the number of real spots).  

C. Training parameters may be changed, but it is recommended to run the CellDART with default parameters.  

D. Click the 'Start' button to start the analysis. Return to 'Upload files - Download' to download predicted cell fraction data.  
