## R wrap function for CellDART
  Please refer to the 'R_example.R' file.  
  All the files in 'CellDART_R_wrap' should be in the same folder.  
  
  1. Using reticulate virtual environment (env.select == 'virtual')  
  When 'pred_cellf_celldart' function defined in 'cellf_pred.R' is running, python or virtual environment will be automatically installed via reticulate.  
  Python will be automatically installed via reticulate::install_python on Linux if python_path = NULL (Tested on Ubuntu 18.04.5 LTS)  
  Python path can be given to install virtual environment by setting python_path = [python version 3.7: python.exe path]  
  (example: python_path = ~\\Anaconda3\\envs\\spatial\\python.exe)  
  Python installed in r-reticulate can be used by setting python_path = "current"; however, the version should be 3.7 for the virtual environement to be successfully installed.  
  
  2. Using conda environment (env.select == 'conda'): recommended for Windows  
  Conda environment should already be installed in Anaconda and it will be loaded for CellDART analysis.  
  Python dependencies: "scanpy==1.5.1","numba==0.52.0","pandas","numpy",
                       "keras==2.3.1","tensorflow==1.14.0","tensorflow-gpu"  

  ** Datasets **
  1. Description  
  Example single-cell Seurat object file: sc_data.rds (GSE115746: mouse from ALS and VISp)  
  Example spatial Seurat object file: sp_data.rds  
  (10X Genomics Data Repository: V1_Mouse_Brain_Sagittal_Anterior, V1_Mouse_Brain_Sagittal_Posterior)  
  
  2. Download  
  sc_data.rds and sp_data.rds can be downloaded from:  
  https://drive.google.com/drive/folders/1cBCeFWvSjxIP1naHBpNZ7nmJV6Ql8OdI?usp=sharing
  
## Potential error in reticulate::install_python
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
