## Function to implement CellDART in python

# adata_sp: spatial data (AnnData object) to be used in predicting cell fraction (default: None)
# -> If None, then provide spdir where spatial datasets are saved (formats are explained below)
# adata_sc: single-cell data (AnnData object) to be used in making pseudospots (default: None)
# -> If None, then provide scdir where single-cell datasets are saved (formats are explained below)
# count_from_raw: whether to extract count matrix frow .raw of AnnData
# -> non-normalized count matrix should be contained in the AnnData .raw file
# -> if False, then utilize the count matrices saved in adata_sp and adata_sc directly

# gpu: check whether to use gpu (True) or not (False) (default = True)

# spdir: file directory to find or save spatial data
# -> In case of utilizing already saved spatial data, otherwise, put None
# -> Visium data should be separated in different folders
## Example directory (spatial)
# -> two spatial datasets (10x visium format)
# ./Mouse_sp/first/filtered_feature_bc_matrix.h5, ./Mouse_sp/first/spatial/tissue_hires_image.png, ./Mouse_sp/first/spatial/tissue_lowres_image.png,
# ./Mouse_sp/first/spatial/scalefactors_json.json, ./Mouse_sp/first/spatial/tissue_positions_list.csv
# second dataset directory starts with ./Mouse_sp/second/.., others are same as above.

# sp10x: whether the spatial data is 10x Visium format (True) or not (False) (default: True)
# spfilter: check whether to filter the number of cells and genes in spatial data (True: run filter)
# spfilgene: keep genes that are expressed in at least 'spfilgene' number of cells (default = 5)
# spfilspot: keep spots with at least 'spfilcell' counts (default = 50)

# scdir: file directory to find or save single-cell data
# -> In case of utilizing already saved sc data, otherwise, put None
# -> each single-cell data should be separated in different folders 
# -> each file formats should be among 10x format or others (.mtx.gz, .h5ad, h5, .csv, .tsv, or .txt)
# -> and metadata with corresponding barcode name as index should be included in metadata folder of each single-cell data
# -> metadata should be csv format
## Example directory (single-cell)
# -> two single cell dataset (10x mtx format) with metadata
# ./Mouse_sc/first/barcodes.tsv, ./Mouse_sc/first/genes.tsv, ./Mouse_sc/first/matrix.mtx, ./Mouse_sc/first/metadata/metadata.csv
# ./Mouse_sc/second/barcodes.tsv, ./Mouse_sc/second/genes.tsv, ./Mouse_sc/second/matrix.mtx, ./Mouse_sc/first/second/metadata.csv

# sc10x_mtx: check whether single-cell data is 10x genomics formatted mtx directory (True) or not (False)
# sc10x_h5: check whether single-cell data is 10x genomics formatted hdf5 file (True) or not (False)
# sctranspose: if sc10x_mtx and sc10x_h5 is F, check whether loaded matrix should be transposed (True) or not (False)

# celltype: column name for single-cell annotation data in .obs (default: 'cluster')
# num_markers: number of selected marker genes in each cell-type (default = 20)

# seed_num: seed to be used in random sampling (default = 0)

# nmix: sampling number of cells in pseudospot (default = 10)
# npseudo: a total number of pseudospots (default = 20000)

# alpha: loss weights of domain classifier to source classifier (default = 0.6)
# alpha_lr: learning rate for domain classifier (alpha_lr*0.001, default = 5)
# batch_size: minibatch size during the training (default = 512)
# emb_dim: output size of dimensions for feature extractor (default = 64)

# n_iterations: iteration number for the adversarial learning (default = 3000)
# init_train_epoch: iteration number of pre-train (default = 10)

# outdir: the directory to save output files (models and results)
# return_anndata: return spatial AnnData file with predicted cell fraction in .obs (default = True)

def pred_cellf_celldart(adata_sp=None, adata_sc=None, count_from_raw=False, 
                        gpu=True, spdir=None, sp10x=True, spfilter=False, spfilgene=5, spfilspot=50, 
                        scdir=None, sc_10x_mtx=True, sc10x_h5=False, sctranspose=False, 
                        celltype='cluster', num_markers=20, seed_num=0, 
                        nmix=10, npseudo=20000, alpha=0.6, alpha_lr=5, batch_size=512, emb_dim=64, n_iterations=3000, init_train_epoch=10, 
                        outdir='./CellDART_output', return_anndata=True):

    import os
    if gpu:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"]= "0" # Use only gpu-0
        print('GPU is available and will be used')
    else:
        os.environ['CUDA_VISIBLE_DEVICES'] = "-1" # Use CPU
        print('CPU will be used')
    
    from warnings import simplefilter 
    simplefilter(action='ignore', category=Warning)

    import scanpy as sc
    import pandas as pd
    import numpy as np

    from CellDART import utils
    from CellDART import da_cellfraction

    ## Change float variables into integer (during conversion from R to python)
    num_markers, seed_num, \
    nmix, npseudo, batch_size, emb_dim, n_iterations, init_train_epoch = \
        int(num_markers), int(seed_num), \
        int(nmix), int(npseudo), int(batch_size), int(emb_dim), \
        int(n_iterations), int(init_train_epoch)
    
    ## Create directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    ## Load and preprocess spatial dataset
    if adata_sp is not None:
        if spdir is not None: 
            raise ValueError("'spdir' should be None when 'adata_sp' is provided.")
        if count_from_raw: spatial_all = adata_sp.raw.to_adata()
        else: spatial_all = adata_sp.copy()     
        sc.pp.normalize_total(spatial_all, target_sum=1e4, inplace=True)
        print('Shape of the provided spatial data is',spatial_all.shape)
    else:
        if spdir is None:
            raise ValueError("'spdir' should be provided when 'adata_sp' is None")
        # Load and normalize spatial data
        sp_list = os.listdir(spdir)
        adata_sp = []
        for i, sp_data in enumerate(sp_list):
            if sp10x:
                adata = sc.read_visium(os.path.join(spdir,sp_data))
            else:
                sp = os.listdir(os.path.join(spdir,sp_data))[0]
                adata = sc.read(os.path.join(spdir,sp_data,sp))
            adata.var_names_make_unique()
            if spfilter:
                sc.pp.filter_genes(adata, min_cells=spfilgene)
                sc.pp.filter_cells(adata, min_counts=spfilspot)
            sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
            adata_sp.append(adata)
            print('Shape of spatial data',i,'is',adata.shape)
            
        # Merge spatial data
        if len(adata_sp)==1:
            spatial_all = adata_sp[0]
        else:
            spatial_all = adata_sp[0].concatenate(adata_sp[1:], join='inner',
                                                    uns_merge='unique')
        print('Shape of the merged spatial data is',spatial_all.shape)

    
    ## Load and preprocess single-cell dataset
    if adata_sc is not None:
        if scdir is not None: 
            raise ValueError("'scdir' should be None when 'adata_sc' is provided.")
        if count_from_raw: single_all = adata_sc.raw.to_adata()
        else: single_all = adata_sc.copy()

        # Check if the column for the cell type is included in .obs
        if celltype not in list(single_all.obs):
            raise ValueError('Column for cell type is not found')
        
        sc.pp.normalize_total(single_all, target_sum=1e4, inplace=True)
        print('Shape of the provided single-cell data is',single_all.shape)
    else:
        if scdir is None:
            raise ValueError("'scdir' should be provided when 'adata_sc' is None")
        # Load single cell data
        sc_list = os.listdir(scdir)
        if sc_10x_mtx:
            adata_sc = [sc.read_10x_mtx(os.path.join(scdir,y), cache=True) for y in sc_list]
        elif sc10x_h5:
            adata_sc = [sc.read_10x_h5(os.path.join(scdir,y)) for y in sc_list]
        else:
            if sctranspose:
                adata_sc = [sc.read(os.path.join(scdir,y,z), cache=True).T for y in sc_list \
                        for z in [i for i in os.listdir(os.path.join(scdir,y)) \
                        if i.endswith('mtx.gz') or i.endswith('h5ad') or i.endswith('h5') or \
                            i.endswith('csv') or i.endswith('tsv') or i.endswith('txt')]]
            else: 
                adata_sc = [sc.read(os.path.join(scdir,y,z), cache=True) for y in sc_list \
                        for z in [i for i in os.listdir(os.path.join(scdir,y)) \
                        if i.endswith('mtx.gz') or i.endswith('h5ad') or i.endswith('h5') or \
                            i.endswith('csv') or i.endswith('tsv') or i.endswith('txt')]]
                
        # preprocess each of the dataset
        for i, adata in enumerate(adata_sc):
            adata.var_names_make_unique()
            sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
            sc_meta_list = os.listdir(os.path.join(scdir,sc_list[i],'metadata'))
            sc_meta_list = [i for i in sc_meta_list if i.endswith('.csv')]
            if len(sc_meta_list)==0:
                raise NotImplementedError('No csv format metadata in the folder')
            tmp = pd.read_csv(os.path.join(scdir,sc_list[i],'metadata',sc_meta_list[0]),index_col=0)
            
            if (set(tmp.index)<=set(adata.obs.index)): 
                print('All barcode names in metadata are found')
            else:
                raise ValueError('Unidentified barcode names in metadata of '+sc_list[i])
            if celltype not in list(tmp):
                raise ValueError('Column for cell type is not found')
            
            # subset the data to include only the barcodes in metadata
            adata = adata[adata.obs.index.isin(tmp.index)].copy()
            # rearrange the metadata index according to adata
            tmp = tmp.reindex(adata.obs.index)
            adata.obs = tmp
            adata_sc[i] = adata
            print('Shape of single cell data',i,'is',adata.shape)
        
        if len(adata_sc)==1:
            single_all = adata_sc[0]
        else:
            single_all = adata_sc[0].concatenate(adata_sc[1:], join='inner')

        print('Shape of merged single cell data is',single_all.shape)
    
    # save the normalized data in raw
    single_all.raw = single_all
    
    # log-transform the count matrix
    sc.pp.log1p(single_all)
    
    # Find marker genes for single cell data
    single_all.obs[celltype] = single_all.obs[celltype].astype('category', copy=False)
    sc.tl.rank_genes_groups(single_all, celltype, method='wilcoxon')
    genelists=single_all.uns['rank_genes_groups']['names']
    df_genelists = pd.DataFrame.from_records(genelists)
    
    # Combining top marker genes representing each cell type
    res_genes = []
    for column in df_genelists.head(num_markers): 
        res_genes.extend(df_genelists.head(num_markers)[column].tolist())
    res_genes_ = list(set(res_genes))
   
    # Calculate intersecting genes
    inter_genes_comb = [val for val in res_genes_ if val in spatial_all.var.index]
    print('Total number of marker genes: ',len(inter_genes_comb))

    
    # Generation of an array representing cell type number
    df_sc = single_all.obs
    lab_sc_sub = df_sc[celltype]
    sc_sub_dict = dict(zip(range(len(set(lab_sc_sub))), set(lab_sc_sub)))
    sc_sub_dict2 = dict((y,x) for x,y in sc_sub_dict.items())
    lab_sc_num = [sc_sub_dict2[ii] for ii in lab_sc_sub]
    # Make an array for cell type numbers following the sequence of single cell barcodes
    lab_sc_num = np.asarray(lab_sc_num, dtype='int')
    
    # Call original normalized count (not log-normalized count)
    adata_final = single_all.raw.to_adata()

    # Generate count matrix for single-cell data (mat_sc)
    adata_final = adata_final[:,inter_genes_comb].copy()
    if isinstance(adata_final.X, np.ndarray):
        mat_sc = adata_final.X
    else:
        mat_sc = adata_final.X.toarray()

    # Raw file for merged spatial data
    spatial_raw = spatial_all
    
    # Generate count matrix for spatial data (mat_sp)
    spatial_all = spatial_all[:,inter_genes_comb].copy()
    if isinstance(spatial_all.X, np.ndarray):
        mat_sp = spatial_all.X
    else: 
        mat_sp = spatial_all.X.toarray()
    
    # Generate pseudospot: random mixture of cells
    sc_mix, lab_mix = utils.random_mix(mat_sc, lab_sc_num, nmix=nmix, n_samples=npseudo, seed=seed_num)
    
    # Log-normalize and scale the data 
    def log_minmaxscale(arr):
        arrd = len(arr)
        arr = np.log1p(arr)
        e = 1e-8 # modified by adding e
        return (arr-np.reshape(np.min(arr,axis=1),(arrd,1)))/np.reshape((np.max(arr,axis=1)-np.min(arr,axis=1))+e,(arrd,1))

    sc_mix_s = log_minmaxscale(sc_mix)
    mat_sp_s = log_minmaxscale(mat_sp)
    mat_sc_s = log_minmaxscale(mat_sc)

    print('Size of spatial, single-cell, pseudospot, and cell fraction data:',
        mat_sp_s.shape, mat_sc_s.shape, sc_mix_s.shape, lab_mix.shape)

 
    # Train the CellDART model
    embs, clssmodel = da_cellfraction.train(sc_mix_s, lab_mix, mat_sp_s, enable_dann = True,
                                            alpha=alpha, alpha_lr=alpha_lr,
                                            emb_dim = emb_dim, 
                                            batch_size = batch_size,
                                            n_iterations = n_iterations,
                                            initial_train=True,
                                            initial_train_epochs=init_train_epoch)
    # Prediction of cell fraction in each spot
    pred_sp = pd.DataFrame(clssmodel.predict(mat_sp_s))
    pred_sp.index = spatial_all.obs.index


    # Make directory for the model save
    if not os.path.exists(os.path.join(outdir,'model')):
        os.makedirs(os.path.join(outdir,'model'))

    # Save the cell fraction in .obs of spatial_raw file
    for visnum in range(len(sc_sub_dict)):
        spatial_raw.obs[str(sc_sub_dict[visnum])+'_cellf'] = pred_sp.iloc[pred_sp.index.isin(spatial_raw.obs.index),visnum]

    # Save cell fraction data
    df = spatial_raw.obs.filter(regex='_cellf', axis=1)
    df.to_csv(os.path.join(outdir,'cellfraction.csv'),header=True,index=True)
    print('Cell fraction data for was saved')

    # Save model files
    embs.save_weights(os.path.join(outdir,'model/embedder.h5'))
    clssmodel.save_weights(os.path.join(outdir,'model/classifier.h5'))

    # Save spatial anndata
    spatial_raw.write_h5ad(os.path.join(outdir,'model/sp_data.h5ad'))
    print('Spatial anndata was saved')

    print('Model and python data files were saved')

    if return_anndata: return(spatial_raw)
    else: return(df)
