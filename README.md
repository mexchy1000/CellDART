# CellDART: Cell type inference by domain adaptation of single-cell and spatial transcriptomic data
CellDART is a tool to estimate cell fraction of spatial transcriptomic spots using domain adaptation of deep neural networks.

![figure1png](https://user-images.githubusercontent.com/14209383/114880774-528b8100-9e3d-11eb-9b60-41c9d0acd5fd.png)


## Optimal parameter choices (for brain)
  Number of total marker genes = 200 ~ 400  
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

## R Wrapper for CellDART using reticulate 
  Please refer to the '/R_wrap/R_example.R' file
  All the files in 'R_wrap' folder
  
## Dependency (R wrapper)
  Seurat 4.0.5  
  dplyr 1.0.7  
  sceasy 0.0.6  
  reticulate 1.22  
  
## R shiny application for CellDART
  ..
  
