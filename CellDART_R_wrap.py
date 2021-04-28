## Implement CellDART in R
import os
import argparse

import scanpy as sc
import pandas as pd
import numpy as np

from utils import random_mix

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--gpu', type=str, default = 'T')
    parser.add_argument('--spdir', type=str)
    parser.add_argument('--spfilter', type=str, default = 'F')
    parser.add_argument('--spfilgene', type=int, default = 5)
    parser.add_argument('--spfilcell', type=int, default = 50)
    parser.add_argument('--scdir', type=str)
    parser.add_argument('--sc10x_mtx', type=str, default='T')
    parser.add_argument('--sc10x_h5', type=str, default='F')
    parser.add_argument('--sctranspose', type=str, default='F')
    parser.add_argument('--celltype', type=str, default = 'cell_type')
    parser.add_argument('--num_markers', type=int, default = 20)
    parser.add_argument('--nmix', type=int, default = 10)
    parser.add_argument('--npseudo', type=int, default = 20000)
    parser.add_argument('--alpha', type=int, default = 1)
    parser.add_argument('--alpha_lr', type=int, default = 5)
    parser.add_argument('--emb_dim', type=int, default = 64)
    parser.add_argument('--batch_size', type=int, default = 512)
    parser.add_argument('--n_iterations', type=int, default = 3000)
    parser.add_argument('--init_train_epoch', type=int, default = 10)
    parser.add_argument('--outdir', type=str, default='./CellDART_output/')
    args = parser.parse_args()
    
    if args.gpu == 'T':
        import tensorflow as tf
        if tf.test.is_gpu_available: 
            print('GPU is available and will be used')
            os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
            os.environ["CUDA_VISIBLE_DEVICES"]= "0" # Use only gpu-0
        else:
            print('GPU is unavailable: CPU will be used')
    else:
        os.environ['CUDA_VISIBLE_DEVICES'] = "-1" # Use CPU
        print('CPU will be used')
    import da_cellfraction
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    # Load 10x data
    sp_list = os.listdir(args.spdir)
    adata_sp = [sc.read_visium(args.spdir+y) for y in sp_list]
    
    # Normalize spatial data
    for i, adata in enumerate(adata_sp):
        adata.var_names_make_unique()
        if args.spfilter == 'T':
            sc.pp.filter_genes(adata, min_cells=args.spfilgene)
            sc.pp.filter_cells(adata, min_counts=args.spfilcell)
        sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
        adata_sp[i] = adata
        print('Shape of spatial data',i,'is',adata.shape)
        
    # Merge spatial data
    if len(adata_sp)==1:
        spatial_all = adata_sp[0]
    else:
        spatial_all = adata_sp[0].concatenate(adata_sp[1:], join='inner',
                                                   uns_merge='unique')
        print('Shape of merged spatial data is',spatial_all.shape)
    
    # Load single cell data
    sc_list = os.listdir(args.scdir)
    if args.sc10x_mtx == 'T':
        adata_sc = [sc.read_10x_mtx(args.scdir+y, cache=True) for y in sc_list]
    elif args.sc10x_h5 == 'T':
        adata_sc = [sc.read_10x_h5(args.scdir+y) for y in sc_list]
    else:
        if args.sctranspose == 'T':
            adata_sc = [sc.read(args.scdir+y+'/'+z, cache=True).T for y in sc_list \
                       for z in [i for i in os.listdir(args.scdir+y) \
                       if i.endswith('mtx.gz') or i.endswith('h5ad') or i.endswith('h5') or \
                        i.endswith('csv') or i.endswith('tsv') or i.endswith('txt')]]
        else: 
            adata_sc = [sc.read(args.scdir+y+'/'+z, cache=True) for y in sc_list \
                       for z in [i for i in os.listdir(args.scdir+y) \
                       if i.endswith('mtx.gz') or i.endswith('h5ad') or i.endswith('h5') or \
                        i.endswith('csv') or i.endswith('tsv') or i.endswith('txt')]]
            
    # preprocess each of the dataset
    for i, adata in enumerate(adata_sc):
        adata.var_names_make_unique()
        sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
        sc_meta_list = os.listdir(args.scdir+sc_list[i]+'/metadata')
        sc_meta_list = [i for i in sc_meta_list if i.endswith('.csv')]
        if len(sc_meta_list)==0:
            print('No csv format metadata in the folder')
            raise NotImplementedError
        tmp = pd.read_csv(args.scdir+sc_list[i]+'/metadata/'+sc_meta_list[0],index_col=0)
        
        if (set(tmp.index)<set(adata.obs.index)): 
            print('All barcode names in metadata are found')
        else:
            print('Unidentified barcode names in metadata of '+sc_list[i])
            raise NotImplementedError
        if args.celltype not in list(tmp):
            print('Column for cell type is not found')
            raise NotImplementedError
        
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
    
    # standard preprocessing steps
    sc.pp.log1p(single_all)
    sc.pp.scale(single_all, max_value=10)
    sc.tl.pca(single_all, svd_solver='arpack')

    if len(adata_sc)==1:
        sc.pp.neighbors(single_all, n_neighbors=10, n_pcs=40)
        sc.tl.umap(single_all)
        sc.tl.leiden(single_all, resolution = 0.5)
        
    else:
        import bbknn
        sc.external.pp.bbknn(single_all, batch_key='batch')
        sc.tl.umap(single_all)
        sc.tl.leiden(single_all, resolution = 0.5)
    
    # Find marker genes for single cell data
    adata_final = single_all.raw.to_adata()
    sc.tl.rank_genes_groups(adata_final, args.celltype, method='wilcoxon')
    genelists=adata_final.uns['rank_genes_groups']['names']
    df_genelists = pd.DataFrame.from_records(genelists)
    
    # Combining top marker genes representing each cell type
    res_genes = []
    for column in df_genelists.head(args.num_markers): 
        res_genes.extend(df_genelists.head(args.num_markers)[column].tolist())
    res_genes_ = list(set(res_genes))
    print('Total number of marker genes is '+str(len(res_genes_)))
    
    # Find overlapping genes
    inter_genes = [val for val in res_genes_ if val in spatial_all.var.index]
    print('Intersecting gene numbers:',len(inter_genes))
    
    adata_final = adata_final[:,inter_genes]
    spatial_all = spatial_all[:,inter_genes]
    
    # Generation of count matrices
    mat_sc = adata_final.X.todense()
    mat_sp = spatial_all.X.todense()
    
    # Generation of an array representing cell type number
    df_sc = adata_final.obs
    lab_sc_sub = df_sc.cluster
    sc_sub_dict = dict(zip(range(len(set(lab_sc_sub))), set(lab_sc_sub)))
    sc_sub_dict2 = dict((y,x) for x,y in sc_sub_dict.items())
    lab_sc_num = [sc_sub_dict2[ii] for ii in lab_sc_sub]
    # Make an array for cell type numbers following the sequence of single cell barcodes
    lab_sc_num = np.asarray(lab_sc_num, dtype='int')
    
    # Generation of pseudospots (select nmix cells to make n_samples pseudospots)
    sc_mix, lab_mix = random_mix(mat_sc, lab_sc_num, nmix=args.nmix, n_samples=args.npseudo)
    
    # Log-normalize and scale the data 
    def log_minmaxscale(arr):
        arrd = len(arr)
        arr = np.log1p(arr)
        return (arr-np.reshape(np.min(arr,axis=1),(arrd,1)))/np.reshape((np.max(arr,axis=1)-np.min(arr,axis=1)),(arrd,1))

    sc_mix_s = log_minmaxscale(sc_mix)
    mat_sp_s = log_minmaxscale(mat_sp)
    mat_sc_s = log_minmaxscale(mat_sc)
    
    print('Size of spatial, single-cell, and pseudospot data:',
          mat_sp_s.shape, mat_sc_s.shape, sc_mix_s.shape)
    
    
    # Train the CellDART model
    embs, clssmodel = da_cellfraction.train(sc_mix_s, lab_mix, mat_sp_s, enable_dann = True,
                                            alpha=args.alpha, alpha_lr=args.alpha_lr, 
                                            emb_dim = args.emb_dim, 
                                            batch_size = args.batch_size,
                                            n_iterations = args.n_iterations,
                                            initial_train=True, 
                                            initial_train_epochs=args.init_train_epoch)
    
    
    # Prediction of cell fraction in each spot
    pred_sp = pd.DataFrame(clssmodel.predict(mat_sp_s))
    pred_sp.index = spatial_all.obs.index
    
    # Save the cell fraction in observation
    for i in range(len(adata_sp)):
        if len(adata_sc) == 1:
            adata = spatial_all
        else:
            adata = spatial_all[spatial_all.obs.batch == str(i), :].copy()
        for visnum in range(len(sc_sub_dict)):
            adata.obs[sc_sub_dict[visnum]] = pred_sp.iloc[pred_sp.index.isin(adata.obs.index),visnum]
        df = adata.obs
        df.to_csv(args.outdir+'cellfraction_'+sp_list[i]+'.csv',header=True,index=True)
        print('Cell fraction data for',sp_list[i],'is saved')
        
    # Save all the files    
    if not os.path.exists(args.outdir+'Model'):
        os.makedirs(args.outdir+'Model')
    embs.save_weights(args.outdir+'Model/embedder.h5')
    clssmodel.save_weights(args.outdir+'Model/classifier.h5')

    # AnnData and dataframe save
    single_all.write_h5ad(args.outdir+'Model/single_cell_proc.h5')
    adata_final.write_h5ad(args.outdir+'Model/single_cell_final.h5')
    spatial_all.write_h5ad(args.outdir+'Model/spatial_final.h5')
    
    # save numpy array as csv file
    np.savetxt(args.outdir+'Model/lab_mix.csv', lab_mix, delimiter=',')
    np.savetxt(args.outdir+'Model/sc_mix.csv', sc_mix, delimiter=',')
    
    # save prediction file
    np.savetxt(args.outdir+'Model/pred_sp.csv', pred_sp, delimiter=',')
    
    dat = pd.DataFrame(list(sc_sub_dict.values()))
    dat.to_csv(args.outdir+'Model/sc_sub_dict.csv', sep = ',',
               header=True, index=False)
    
    print('Model and python data files are saved')
    
