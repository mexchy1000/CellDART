#' R wrap function to implement CellDART
#' @description Cell type inference by domain adaptation of single-cell and spatial transcriptomic data
#'
#' @param sp_data spatial data (AnnData object) to be used in predicting cell fraction (default: None): count matrix should be non-normalied raw data.
#' @param sc_data single-cell data (AnnData object) to be used in making pseudospots (default: None): count matrix should be non-normalied raw data.
#'
#' @param outdir the directory to save output files (models and results) (default = '.')
#'
#' @param sp_subset whether to subset spatial data and calculate for specific spot cluster (default = FALSE)
#' @param spot.cluster.name group name of the cluster used for subsetting spatial data (default: 'seurat_clusters')
#' @param spot.cluster.of.interest name of each spot clusters to be used (default: NULL)
#' @param metadata_celltype column name for single-cell annotation data in metadata (default: 'celltype')
#'
#' @param env.select select between using reticulate virtual environment or conda environment (default: "conda")
#' @param python.install whether to automatically install python version 3.7.12
#'
#' @param python_path path for the python 3.7. (default: NULL)
#' \itemize{
#'   \item If NULL, python version 3.7.12 will be installed (valid for Linux)
#'   \item If "current", python interpreter associated with current virtual env (ex: r-reticulate) will be used. (version should be 3.7)
#' }
#'
#' @param env.name name of the virtual or conda environment to use for CellDART analysis (default: 'spSeudoMap')
#'
#' @param gpu check whether to use gpu (True) or not (False) (default = True)
#' @param metadata_celltype column name for single-cell annotation data in metadata (default: 'celltype')
#' @param num_markers number of selected marker genes in each cell-type (default = 20)
#' @param seed_num seed to be used in random sampling (default = 0)
#' @param nmix the number of cells sampled from single-cell data when making a pseudospot (default = 10)
#' @param npseudo a total number of pseudospots (default = 20000)
#'
#' @param alpha loss weights of domain classifier to the source classifier (default = 0.6)
#' @param alpha_lr learning rate for the domain classifier (alpha_lr*0.001, default = 5)
#' @param emb_dim output size of dimensions for feature extractor (default = 64)
#'
#' @param batch_size minibatch size for pseudospots and spatial data during the training (default = 512)
#' @param n_iterations iteration number for the adversarial learning (default = 3000)
#' @param init_train_epoch iteration number for the pre-training process (default = 10)
#'
#' @return spatial data (Seurat object) with predicted cell fraction in metadata (meta.data)
#' @examples
#' Using conda environment (environment will be automatically installed in Linux distributions)
#' If using Windows, then install conda environment first and then run the function below with python.install = F
#' sp_data_cellf <- pred_cellf_celldart(sp_data, sc_data, outdir = '.',
#'                                      sp_subset=F, spot.cluster.name='seurat_clusters',
#'                                      spot.cluster.of.interest=NULL,
#'                                      env.select='conda',python.install=T,
#'                                      python_path=NULL, env.name='CellDART',
#'                                      gpu=TRUE, metadata_celltype='celltype',
#'                                      num_markers=20, seed_num=0,
#'                                      nmix=8, npseudo=20000, alpha=0.6,alpha_lr=5,
#'                                      emb_dim=64,batch_size=512,n_iterations=3000, init_train_epoch=10)
#' 
#' Using virtual environment (environment will be automatically installed in Linux distributions)
#' Not recommended for Windows
#' sp_data_cellf <- pred_cellf_celldart(sp_data, sc_data, outdir = '.',
#'                                      sp_subset=F, spot.cluster.name='seurat_clusters',
#'                                      spot.cluster.of.interest=NULL,
#'                                      env.select='virtual',python.install=T,
#'                                      python_path=NULL, env.name='CellDART',
#'                                      gpu=TRUE, metadata_celltype='celltype',
#'                                      num_markers=20, seed_num=0,
#'                                      nmix=8, npseudo=20000, alpha=0.6,alpha_lr=5,
#'                                      emb_dim=64,batch_size=512,n_iterations=3000, init_train_epoch=10)
#' @export 
pred_cellf_celldart <- function(sp_data, sc_data, outdir='.',
                                sp_subset=FALSE, spot.cluster.name='seurat_clusters',
                                spot.cluster.of.interest=NULL,
                                env.select='virtual', python.install=F,
                                python_path=NULL, env.name='CellDART',
                                gpu=TRUE, metadata_celltype='celltype',
                                num_markers=20, seed_num=0,
                                nmix=8, npseudo=20000, alpha=0.6, alpha_lr=5,
                                emb_dim=64, batch_size=512, n_iterations=3000, 
                                init_train_epoch=10){
  # Suppress warnings
  defaultW <- getOption("warn") 
  options(warn = -1)

  if (python.install){
    reticulate::install_python(version = '3.7.12')
  }
  
  # Select between using reticulate virtual environment or conda environment ("virtual" or "conda")
  if (env.select=="virtual"){
    # Setting virtual environment with reticulate  
    if (!(virtual.env.name %in% reticulate::virtualenv_list())){
      ## Python dependencies use python version 3.7
      if (is.null(python_path)){
        reticulate::virtualenv_create(envname = env.name, version = '3.7.12')
      } else if (python_path=="current") {
        reticulate::virtualenv_create(envname = env.name, python = NULL)
      } else {
        reticulate::virtualenv_create(envname = env.name, python = python_path)
      }
      # python_depend = c("scanpy==1.5.1","numba==0.52.0","pandas","numpy",
      #                   "keras==2.3.1","tensorflow==1.14.0","tensorflow-gpu")

      # Create virtual env and install dependencies
      reticulate::virtualenv_install(env.name, packages = 'CellDART', ignore_installed=T, 
                                     pip_options = "git+https://github.com/mexchy1000/CellDART.git")
      reticulate::use_virtualenv(env.name, required = T)
    }
    # Apply virtual environment
    reticulate::use_virtualenv(env.name, required = T)
  } else if (env.select=="conda"){
    if (!(env.name %in% reticulate::conda_list()[['name']])){
      ## Python dependencies use python version 3.7
      if (is.null(python_path)){
        reticulate::conda_create(envname = env.name, python_version = '3.7.12')
      } else if (python_path=="current") {
        reticulate::conda_create(envname = env.name, python = NULL)
      } else {
        reticulate::conda_create(envname = env.name, python = python_path)
      }

      # Create conda env and install dependencies
      reticulate::conda_install(env.name, packages='CellDART', ignore_installed=T,
                                pip = TRUE, "git+https://github.com/mexchy1000/CellDART.git")
    }
    # Apply conda environment
    reticulate::use_condaenv(env.name, required = T)
  } else {
    stop("'env.select' should be either 'virtual' or 'conda'")
  }

  ## Import anndata
  ann <- reticulate::import('anndata', convert = FALSE)

  ## Import python function
  CellDART <- reticulate::import('CellDART', convert = FALSE)
  
  ## 1. Saving single-cell data in anndata format
  # Define count matrix
  sparse_mtx <- Seurat::GetAssayData(sc_data, slot = "counts", assay = "RNA")

  # Define obs and var (reference from sceasy library: https://github.com/cellgeni/sceasy)
  obs <- sc_data@meta.data
  if (!metadata_celltype %in% colnames(obs)){
    stop("Column name for the cell annotation should be provided.")
  } else {
    obs <- obs[metadata_celltype]
    obs[[metadata_celltype]] <- factor(obs[[metadata_celltype]])
  }
  var <- data.frame(matrix(nrow=dim(sc_data)[1],ncol=0,
                           dimnames = list(rownames(sc_data),NULL)))
  var[['name']] <- rownames(var)

  adata_sc <- ann$AnnData(
    X = Matrix::t(sparse_mtx),
    obs = obs,
    var = var
  )

  ## 2. Subsetting spatial data and save in anndata format
  if (sp_subset){
    cluster_info <- sp_data[[spot.cluster.name]]
    Seurat::Idents(sp_data) <- spot.cluster.name
  }

  if (is.null(spot.cluster.of.interest)){
    sp_data_sub <- sp_data
  } else if (sum(spot.cluster.of.interest%in%levels(cluster_info))==length(spot.cluster.of.interest)){
    sp_data_sub <- subset(sp_data, idents=spot.cluster.of.interest)
  } else {
    stop("'spot.cluster.of.interest' should be among the levels of 'spot.cluster.name' provided")
  }

  # Define count matrix
  sparse_mtx <- Seurat::GetAssayData(sp_data_sub, slot = "counts", assay = "Spatial")

  # Define obs and var (reference from sceasy library)
  obs <- sp_data_sub@meta.data
  var <- data.frame(matrix(nrow=dim(sp_data_sub)[1],ncol=0,
                           dimnames = list(rownames(sp_data_sub),NULL)))
  var[['name']] <- rownames(var)

  adata_sp <- scanpy_data$AnnData(
    X = Matrix::t(sparse_mtx),
    obs = obs,
    var = var
  )

  # Assign the output directory for the models generated
  if (!file.exists(outdir)){
    dir.create(file.path(outdir, 'results'))
  }
  out_dir <- file.path(getwd(), outdir, 'results')
  
  # Run CellDART
  try({
    df <- CellDART$pred_cellf_celldart$pred_cellf_celldart(adata_sp=adata_sp, adata_sc=adata_sc, count_from_raw=False, 
                                                           gpu=gpu, celltype=metadata_celltype, num_markers=num_markers,
                                                           nmix=nmix, npseudo=npseudo, alpha=alpha, alpha_lr=alpha_lr, 
                                                           batch_size=batch_size, emb_dim=emb_dim, n_iterations=n_iterations,
                                                           init_train_epoch=init_train_epoch, 
                                                           outdir=out_dir, return_anndata=False)

    # Saving cell fraction data into the metadata of spatial Seurat object
    sp_data_sub <- Seurat::AddMetaData(sp_data_sub, reticulate::py_to_r(df))
  })

  options(warn = defaultW)
  return(sp_data_sub)
}
