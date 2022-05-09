
library(SeuratObject)
library(dplyr)


### CellDART implementation

## Explanation of the variables
# celldart.dir: R wrapper python file for the celldart analysis
# outdir: file directory to save the model and result
# sp_data: spatial data (Seurat object) to be used in predicting cell fraction (default: NULL)
# -> If NULL, then use the spatial data already saved in spdir
# sc_data: single-cell data (Seurat object) to be used in making pseudospots (default: NULL)
# -> If NULL, then use the spatial data already saved in scdir

# spdir: file directory to find or save spatial data
# -> In case of utilizing already saved spatial data, otherwise, put NULL
# -> Visium data should be separated in different folders
## Example directory (spatial)
# -> two spatial datasets (10x visium format)
# ./Mouse_sp/first/filtered_feature_bc_matrix.h5, ./Mouse_sp/first/spatial/tissue_hires_image.png, ./Mouse_sp/first/spatial/tissue_lowres_image.png,
# ./Mouse_sp/first/spatial/scalefactors_json.json, ./Mouse_sp/first/spatial/tissue_positions_list.csv
# second dataset directory starts with ./Mouse_sp/second/.., others are same as above.

# scdir: file directory to find or save single-cell data
# -> In case of utilizing already saved sc data, otherwise, put NULL
# -> each single-cell data should be separated in different folders 
# -> each file formats should be among 10x format or others (.mtx.gz, .h5ad, h5, .csv, .tsv, or .txt)
# -> and metadata with corresponding barcode name as index should be included in metadata folder of each single-cell data
# -> metadata should be csv format
## Example directory (single-cell)
# -> two single cell dataset (10x mtx format) with metadata
# ./Mouse_sc/first/barcodes.tsv, ./Mouse_sc/first/genes.tsv, ./Mouse_sc/first/matrix.mtx, ./Mouse_sc/first/metadata/metadata.csv
# ./Mouse_sc/second/barcodes.tsv, ./Mouse_sc/second/genes.tsv, ./Mouse_sc/second/matrix.mtx, ./Mouse_sc/first/second/metadata.csv

# sp_subset: whether to subset spatial data and calculate for specific spot cluster (default = FALSE)
# spot.cluster.name: group name of the cluster used for subsetting spatial data (default: 'seurat_clusters')
# spot.cluster.of.interest: name of each spot clusters to be used (default: NULL)
# metadata_celltype: column name for single-cell annotation data in metadata (default: 'celltype')

# env.select: select between using reticulate virtual environment or conda environment ("virtual" or "conda")
# 1. In case of env.select == 'virtual'
# python_path: path for the python 3.7. (default: NULL)
# If NULL, python version 3.7.9 will be installed (valid for Linux) 
# If "current", python interpreter associated with current virtual env (ex: r-reticulate) will be used. (version should be 3.7)
# virtual.env.name: name of the virtual environment to use for CellDART analysis (default: 'CellDART')

# 2. In case of env.select == 'conda'
# conda.evn.name: name of the conda environment (already installed with Anaconda) to use for CellDART analysis (default: 'CellDART')

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

# alpha: loss weights of domain classifier to source classifier (default = 0.6)
# alpha_lr: learning rate for domain classifier (alpha_lr*0.001, default = 5)
# batch_size: minibatch size during the training (default = 512)
# emb_dim: output size of dimensions for feature extractor (default = 64)

# n_iterations: iteration number for the adversarial learning (default = 3000)
# init_train_epoch: iteration number of pre-train (default = 10)

pred_cellf_celldart <- function(celldart.dir, outdir,sp_data=NULL,sc_data=NULL,
                                spdir=NULL,scdir=NULL,
                                sp_subset=T,spot.cluster.name='seurat_clusters',
                                spot.cluster.of.interest=NULL,metadata_celltype='celltype',
                                env.select='virtual',
                                python_path=NULL,virtual.env.name='CellDART',
                                conda.env.name='CellDART',
                                gpu=T,sp10x=T,spfilter=F,spfilgene=0,
                                spfilspot=0,sc10x_mtx=F,sc10x_h5=F,sctranspose=F,
                                num_markers=20,seed_num=0,
                                nmix=8,npseudo=20000,alpha=0.6,alpha_lr=5,
                                emb_dim=64,batch_size=512,n_iterations=3000,init_train_epoch=10
                                ){
  # Suppress warnings
  defaultW <- getOption("warn") 
  options(warn = -1)
  
  # Select between using reticulate virtual environment or conda environment ("virtual" or "conda")
  if (env.select=="virtual"){
    # Setting virtual environment with reticulate
    # Check if python 3.7.9 is installed
    if (is.null(python_path)){
      reticulate::install_python(version = '3.7.12')
    }
  
    if (!(virtual.env.name %in% reticulate::virtualenv_list())){
      ## Python dependencies use python version 3.7
      if (is.null(python_path)){
        reticulate::virtualenv_create(envname = virtual.env.name, version = '3.7.12')
      } else if (python_path=="current") {
        reticulate::virtualenv_create(envname = virtual.env.name, python = NULL)
      } else {
        reticulate::virtualenv_create(envname = virtual.env.name, python = python_path)
      }
      # python_depend = c("scanpy==1.5.1","numba==0.52.0","pandas","numpy",
      #                   "keras==2.3.1","tensorflow==1.14.0","tensorflow-gpu")

      # Create virtual env and install dependencies
      reticulate::virtualenv_install(virtual.env.name, packages = 'CellDART', ignore_installed=T, 
                                     pip_options = "git+https://github.com/mexchy1000/CellDART.git")
      reticulate::use_virtualenv(virtual.env.name, required = T)
    }
    # Apply virtual environment
    reticulate::use_virtualenv(virtual.env.name, required = T)
  } else if (env.select=="conda"){
    # Apply conda environment
    reticulate::use_condaenv(conda.env.name, required = T)
  } else {
    stop("'env.select' should be either 'virtual' or 'conda'")
  }
  
  
  if (!is.null(sc_data)){
    ## Create file directory for single-cell data
    sc_dir <- file.path(outdir,'sc_data')
    if (!file.exists(sc_dir)){
      dir.create(sc_dir)
      dir.create(file.path(sc_dir,'a'))
      dir.create(file.path(sc_dir,'a','metadata'))
    } else if (!file.exists(file.path(sc_dir,'a'))){
      dir.create(file.path(file.path(sc_dir,'a')))
      dir.create(file.path(sc_dir,'a','metadata'))
    } else if (!file.exists(file.path(sc_dir,'a','metadata'))){
      dir.create(file.path(sc_dir,'a','metadata'))
    }
    
    ## 1. Save single-cell cell-type annotation data
    meta.data <- eval(parse(text=paste0('sc_data@meta.data %>% dplyr::select(',metadata_celltype,')')))
    write.csv(meta.data, file.path(sc_dir,'a','metadata','sc_metadata.csv'),
              row.names = T)
  
    ## 2. Saving single-cell data in anndata format
    sceasy::convertFormat(sc_data, from="seurat", to="anndata",
                          main_layer = 'counts',
                          outFile=file.path(sc_dir,'a','sc_data.h5ad'))
  } else {
    if (!is.null(scdir)){
      sc_dir <- scdir
    } else {
      stop("Directory for single-cell data should be provided in 'scdir'")
    }
  }
  
  if (!is.null(sp_data)){
    ## Create file directory for spatial data
    sp_dir = file.path(outdir,'sp_data')
    if (!file.exists(sp_dir)){
      dir.create(sp_dir)
      dir.create(file.path(sp_dir,'a'))
    } else if (!file.exists(sp_dir,'a')){
      dir.create(file.path(sp_dir,'a'))
    }
    
    ## 3. Subsetting spatial data and save
    cluster_info <- eval(parse(text=paste0('sp_data$',spot.cluster.name)))
    Idents(sp_data) <- spot.cluster.name
    
    if (is.null(spot.cluster.of.interest)){
      sp_data_sub <- sp_data
      sceasy::convertFormat(sp_data_sub, from="seurat", to="anndata",
                            assay = 'Spatial', main_layer = 'counts',
                            outFile = file.path(sp_dir,'a','sp_data.h5ad'))
    } else if (sum(spot.cluster.of.interest%in%levels(cluster_info))==length(spot.cluster.of.interest)){
      sp_data_sub <- subset(sp_data, idents=spot.cluster.of.interest)
      sceasy::convertFormat(sp_data_sub, from="seurat", to="anndata",
                            assay = 'Spatial', main_layer = 'counts',
                            outFile = file.path(sp_dir,'a','sp_data.h5ad'))
    } else {
      stop("'spot.cluster.of.interest' should be among the levels of 'spot.cluster.name' provided")
    }
  } else {
    if (!is.null(spdir)){
      sp_dir <- spdir
    } else {
      stop("Directory for spatial data should be provided in 'spdir'")
    }
  }
  
  ## RUN CellDART and save the result in Seurat object
  # Assign the output directory for the models generated
  out_dir <- file.path(outdir, 'results')
  
  # Change the type of boolean input
  gpu <- ifelse(gpu, 'T', 'F')
  sp10x <- ifelse(sp10x, 'T', 'F')
  spfilter <- ifelse(spfilter, 'T', 'F')
  sc10x_mtx <- ifelse(sc10x_mtx, 'T', 'F')
  sc10x_h5 <- ifelse(sc10x_h5, 'T', 'F')
  sctranspose <- ifelse(sctranspose, 'T', 'F')
  
  path_os <- paste0("python ",celldart.dir," --gpu ",gpu," --sp10x ",sp10x,
                    " --spdir ",sp_dir," --spfilter ",spfilter,
                    " --spfilgene ",spfilgene," --spfilspot ",spfilspot,
                    " --scdir ",sc_dir," --sc10x_mtx ",sc10x_mtx,
                    " --sc10x_h5 ",sc10x_h5," --sctranspose ",sctranspose,
                    " --celltype ",metadata_celltype," --num_markers ",num_markers,
                    " --seed_num ",seed_num,
                    " --nmix ",nmix,
                    " --npseudo ",npseudo," --alpha ",alpha," --alpha_lr ",alpha_lr,
                    " --batch_size ",batch_size," --emb_dim ",emb_dim,
                    " --n_iterations ",n_iterations,
                    " --init_train_epoch ",init_train_epoch,
                    " --outdir ",out_dir)
  
  # RUN CellDART
  system(path_os)
  
  # Loading cell fraction data
  results <- list.files(file.path(out_dir))
  for (i in 1:length(results)){
    if (grepl('csv',results[i])&grepl('cellf',results[i])){
      cellf <- read.csv(file.path(out_dir,results[i]),header=TRUE,row.names='X')  
    }
  }
  if (!is.null(sp_data)){
    # Saving cell fraction data into the file
    sp_data_sub <- Seurat::AddMetaData(sp_data_sub, cellf)
    return(sp_data_sub)
  }
  options(warn = defaultW)  
}
