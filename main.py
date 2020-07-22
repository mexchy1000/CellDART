import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import transfer_model
from sklearn.manifold import TSNE

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3

adata_spatial_anterior = sc.datasets.visium_sge(
    sample_id="V1_Mouse_Brain_Sagittal_Anterior"
)
adata_spatial_posterior = sc.datasets.visium_sge(
    sample_id="V1_Mouse_Brain_Sagittal_Posterior"
)

#Normalize and log1P
for adata in [
    adata_spatial_anterior,
    adata_spatial_posterior,
]:
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace=True)
    
#scale factor
adata_spatial_anterior.uns['spatial']['V1_Mouse_Brain_Sagittal_Anterior']['scalefactors']

##################
#Sc data GSE115746

adata_cortex = sc.read_csv('../data/GSE115746_cells_exon_counts.csv').T
adata_cortex_meta = pd.read_csv('../data/GSE115746_complete_metadata_28706-cells.csv', index_col=0)
adata_cortex_meta_ = adata_cortex_meta.loc[adata_cortex.obs.index,]

adata_cortex.obs = adata_cortex_meta_

adata_cortex.var_names_make_unique()  
sc.pl.highest_expr_genes(adata_cortex, n_top=20, )

adata_cortex.var['mt'] = adata_cortex.var_names.str.startswith('Mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_cortex, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata_cortex, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata_cortex, x='total_counts', y='n_genes_by_counts')

sc.pp.normalize_total(adata_cortex)
sc.pp.log1p(adata_cortex)
sc.pp.highly_variable_genes(adata_cortex,flavor="seurat", n_top_genes=2000, inplace=True)
sc.pl.highly_variable_genes(adata_cortex)

#PCA and clustering
sc.tl.pca(adata_cortex, svd_solver='arpack')
sc.pp.neighbors(adata_cortex, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_cortex)
sc.tl.leiden(adata_cortex, resolution = 0.5)
sc.pl.umap(adata_cortex, color=['leiden','cell_subclass'])


#Int Genes
adata_spatial_anterior.var_names_make_unique() 
inter_genes = [val for val in adata_cortex.var.index if val in adata_spatial_anterior.var.index]
adata_cortex = adata_cortex[:,inter_genes]

adata_spatial_anterior = adata_spatial_anterior[:,inter_genes]

#####
#To arrays#
###########
mat_sc = adata_cortex.X
mat_sp = adata_spatial_anterior.X.todense()

df_sc = adata_cortex.obs
lab_sc = np.asarray(df_sc.leiden, dtype='int')

lab_sc_sub = df_sc.cell_subclass
sc_sub_dict = dict(zip(range(len(set(lab_sc_sub))), set(lab_sc_sub)))
sc_sub_dict2 = dict((y,x) for x,y in sc_sub_dict.items())
lab_sc_num = [sc_sub_dict2[ii] for ii in lab_sc_sub]
lab_sc_num = np.asarray(lab_sc_num, dtype='int')


embs, clssmodel = transfer_model.train(mat_sc, lab_sc_num, mat_sp, enable_dann = True,
                                 alpha=1, alpha_lr=100, emb_dim = 256, batch_size = 1024,
                                 n_iterations = 2000)

#Predicted Embedding
z_sc = embs.predict(mat_sc)
z_sp = embs.predict(mat_sp)

z_scsp= np.concatenate([z_sc,z_sp], axis=0)
z_tsne = TSNE(n_components=2).fit_transform(z_scsp)

pred_sc = clssmodel.predict(mat_sc)
pred_sc_ = np.argmax(pred_sc, axis=1)

pred_sp = clssmodel.predict(mat_sp)
pred_sp_ = np.argmax(pred_sp, axis=1)

sns.scatterplot(x=z_tsne[:,0], y= z_tsne[:,1],
                hue = [0]*z_sc.shape[0]+[1]*z_sp.shape[0], alpha=0.1, size=1,
                linewidth=0)
sns.scatterplot(x=z_tsne[:z_sc.shape[0],0], 
                y= z_tsne[:z_sc.shape[0],1],
                hue = [str(i) for i in lab_sc_sub.tolist()],
                palette = 'Set1',
                alpha=0.1, size=1,
                linewidth=0)

sns.scatterplot(x=z_tsne[:z_sc.shape[0],0], 
                y= z_tsne[:z_sc.shape[0],1],
                hue = [str(i) for i in pred_sc_.tolist()],
                palette = 'Set1',
                alpha=0.1, size=1,
                linewidth=0)

sns.scatterplot(x=z_tsne[z_sc.shape[0]:,0], 
                y= z_tsne[z_sc.shape[0]:,1],
                hue = [str(i) for i in pred_sp_.tolist()],
                palette = 'Set1',
                alpha=0.5, size=1,
                linewidth=0)

#Score for specific cell types.. 
adata_spatial_anterior.obs['Pred_label'] = pred_sp[:,1]
sc.pl.spatial(
        adata_spatial_anterior,
        img_key="hires",
        color='Pred_label',
        palette='Set1',
        size=1.5,
        legend_loc=None)


#Or
#adata = AnnData(
#    pd.read_csv('output/X.csv').asarray(),
#    pd.read_csv('output/obs.csv'),
#    pd.read_csv('output/var.csv'),
#    { # adata.uns
#        'some_thing': pd.read_csv('output/some_thing.csv'),
#    },
#    pd.read_csv('output/obsm.csv'),
#    pd.read_csv('output/varm.csv'),
#)