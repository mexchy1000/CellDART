## R wrap function for CellDART
  ```Plain Text
  devtools::install_github("mexchy1000/CellDART", build_vignettes = T, force = T)  
  library(CellDART)  
  help(pred_cellf_celldart)  # Explanation for the parameters and short examples  
  browseVignettes("CellDART")  # Browse for the vignettes (/doc/introduction.html)
  ```
  ### Installation in Linux distributions  
  Virtual environment (env.select="virtual") or conda environment (env.select="conda") will be automatically installed while running function 'pred_cellf_celldart'  
  ### Installation in Windows  
  Install conda environment first and then run the function with env.select='conda' and python.install=F   
  Example: Please refer to the '/vignettes/introduction.Rmd' file.  

### Datasets
  #### 1. Description  
  Example single-cell Seurat object file: sc_data.rds (GSE115746: mouse from ALS and VISp)  
  Example spatial Seurat object file: sp_data.rds  
  (10X Genomics Data Repository: V1_Mouse_Brain_Sagittal_Anterior, V1_Mouse_Brain_Sagittal_Posterior)  
  
  #### 2. Download  
  sc_data.rds and sp_data.rds can be downloaded from:  
  https://drive.google.com/drive/folders/1zvn3rJ6vH-LYPNnuWIjqb8USSfp7YfIy?usp=sharing
  
### Potential error in reticulate::install_python
  "ModuleNotFoundError: No module named '_ctypes'"  
  Then try on the below command (suggested from https://stackoverflow.com/questions/27022373)  
  sudo apt-get -y update  
  sudo apt-get -y upgrade  
  sudo apt-get -y dist-upgrade  
  sudo apt-get -y install build-essential python-dev python-setuptools python-pip python-smbus  
  sudo apt-get -y install libncursesw5-dev libgdbm-dev libc6-dev  
  sudo apt-get -y install zlib1g-dev libsqlite3-dev tk-dev  
  sudo apt-get -y install libssl-dev openssl  
  sudo apt-get -y install libffi-dev  
