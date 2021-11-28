## R wrap function for CellDART
  Please refer to the 'R_example.R' file.  
  All the files in 'CellDART_R_wrap' should be in the same folder.  
  When 'pred_cellf_celldart' function defined in 'cellf_pred.R' is running, python or virtual environment will be automatically installed via reticulate.  
  Python will be automatically installed via reticulate::python_install on Linux if python_path = NULL (Tested on Ubuntu 18.04.5 LTS)  
  Python path can be used to install virtual environment by python_path = [python version 3.7: python.exe path]  
  Python installed for r-reticulate can be used by python_path = "current"; however, version should be 3.7 for virtual environement should be successfully installed.  

  ** Datasets **
  1. Description  
  Example single-cell Seurat object file: sc_data.rds (GSE115746: mouse from ALS and VISp)  
  Example spatial Seurat object file: sp_data.rds  
  (10X Genomics Data Repository: V1_Mouse_Brain_Sagittal_Anterior, V1_Mouse_Brain_Sagittal_Posterior)  
  
  2. Download  
  sc_data.rds and sp_data.rds can be downloaded from:  
  https://drive.google.com/drive/folders/1cBCeFWvSjxIP1naHBpNZ7nmJV6Ql8OdI?usp=sharing
