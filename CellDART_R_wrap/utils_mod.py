#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from keras.utils import to_categorical


## Make random mix for pseudospot
# Return subsetted array for single-cell data, composite gene expression array for pseudospot (Xs_new), cell fraction (ys_new)
# Also, return gene expression by cell type (zs) if feat_by_celltype = True

# adata_sc: original single cell data that contains all the features
# ys: array for cell type numbers following the sequence of single cell barcodes
# genes_of_interest: genes of interest to generate pseudospot mixture
# nmix: sampling number of cells to make pseudospots (default: nmix=8)
# n_samples: number of pseudospots to generate (default: 20000)
# seed: set seed number during the module score generation (default: seed=0)

def random_mix(adata_sc, ys, genes_of_interest, nmix=8, n_samples=20000, seed=0):

    Xs_new, ys_new, zs = [], [], []
    
    # Convert single-cell anndata to array
    adata_sc_sub = adata_sc[:,genes_of_interest].copy()
    if isinstance(adata_sc.X, np.ndarray):
        Xs = adata_sc_sub.X
    else:
        Xs = adata_sc_sub.X.toarray()
        
    ys_ = to_categorical(ys)
    
    rstate = np.random.RandomState(seed)
    fraction_all = rstate.rand(n_samples,nmix)
    randindex_all = rstate.randint(len(Xs),size=(n_samples,nmix))
    
    for i in range(n_samples):
        # Random fraction of the cell-type & pseudo cell-type
        fraction = fraction_all[i]
        fraction = fraction/np.sum(fraction)
        fraction = np.reshape(fraction, (nmix,1))
        
        # Random selection of the single cell data
        randindex = randindex_all[i]
        ymix = ys_[randindex]
        # Calculate weighted cell fraction
        yy = np.sum(ymix*np.reshape(fraction,(nmix,1)), axis=0)
        # Calculate weighted gene expression
        XX = np.asarray(Xs[randindex])*fraction
        XX_ = np.sum(XX, axis=0)
        
        # Add cell fraction & gene expression
        ys_new.append(yy)
        Xs_new.append(XX_)

    Xs_new = np.asarray(Xs_new)
    ys_new = np.asarray(ys_new)
    
    return Xs, Xs_new, ys_new
