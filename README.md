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
Example File: CellDART_example_mousebrain_markers.ipynb  

## Dependency (python)
  json 2.0.9  
  numpy 1.17.0  
  pandas 0.25.0  
  tensorflow 1.14.0  
  ipywidgets 7.5.1  
  scanpy 1.5.1  
  pandas 0.25.0  
  numpy 1.17.0  
  seaborn 0.10.1  
  keras 2.2.4  

## R wrapper for CellDART using reticulate
  Please refer to the '/CellDART_R_wrap/R_example.R' file.  
  All the files in 'CellDART_R_wrap' should be in the same folder.  
  Reticulate virtual environment or conda environment can be used while running function 'pred_cellf_celldart'  
  Detailed explanation is in '/CellDART_R_wrap/Read_R_wrap.md' file.  
  
## Dependency (R wrapper)
  Seurat 4.0.5  
  dplyr 1.0.7  
  sceasy 0.0.6  
  reticulate 1.22  
  
## R shiny application for CellDART
Shiny application for preprocessing and CellDART analysis.  
Application is based on Seurat, sceasy, dplyr, stringr, vroom, and CellDART.  
The web application panel consists of main, upload, and analysis sections.  
Please try among the four URLs for the application **(ID: CellDART, PW: user) ** 
URL 1: http://34.64.214.85:3838/CellDART1  
URL 2: http://34.64.214.85:3838/CellDART2/  
URL 3: http://34.64.214.85:3838/CellDART3/  
URL 4: http://34.64.214.85:3838/CellDART4/  

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
